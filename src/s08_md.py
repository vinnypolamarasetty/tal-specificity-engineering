"""
Step 08 - MD relaxation and pocket-stability check.

Docking scores a rigid receptor. That is the single biggest approximation in
this pipeline, and it is exactly the approximation that fails for the variants
we most want to trust: a substitution can score well simply because the
grafted rotamer happened to sit out of the way, while in reality it collapses
the pocket or is too mobile to hold the substrate in a productive pose.

This step runs short implicit-solvent MD on the wild type and on shortlisted
variants, and reports three things per system:

  * backbone RMSD        - did the fold stay put, or is the mutant unstable?
  * active-site RMSF     - is the pocket rigid, or is the new side chain
                           flapping and unlikely to give clean selectivity?
  * pocket width         - did the substitution narrow or open the channel?

Implicit solvent (GBn2) is used so this is tractable on a CPU. That trades
solvent structure for speed, which is acceptable for a stability triage but
would not be acceptable for binding free energies - those would need explicit
solvent and an alchemical method, which is out of scope here and is called
out as such.

Usage:
    python src/s08_md.py --structure data/structures/2O6Y_clean.pdb --ns 0.5
    python src/s08_md.py --all-mutants --ns 0.25 --site 89
"""

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
from openmm import LangevinMiddleIntegrator, Platform, unit
from openmm.app import (CutoffNonPeriodic, DCDReporter, ForceField, HBonds,
                        NoCutoff, PDBFile, Simulation, StateDataReporter)

from config import MD_DIR, MD_TEMPERATURE, STRUCT_DIR

IMPLICIT = ("amber14-all.xml", "implicit/gbn2.xml")


def run_md(pdb_path, ns, site_resid, temperature=MD_TEMPERATURE,
           timestep_fs=2.0, equil_ps=20.0, report_ps=2.0, seed=42,
           cutoff_nm=2.0):
    """Run one implicit-solvent trajectory. Returns a metrics dict.

    `cutoff_nm` matters more than it looks. Generalised Born with no cutoff is
    an O(N^2) calculation: fine for a small monomer, but an ammonia-lyase is a
    homotetramer of ~500-residue subunits, and at that size an uncut run is
    hours per nanosecond on a CPU. A 2 nm cutoff makes the same system
    tractable at a modest cost in long-range accuracy. Pass cutoff_nm=None to
    disable it.
    """
    pdb_path = Path(pdb_path)
    name = pdb_path.stem.replace("_clean", "")
    print(f"\n=== MD: {name} ===")

    pdb = PDBFile(str(pdb_path))
    n_atoms = pdb.topology.getNumAtoms()
    forcefield = ForceField(*IMPLICIT)

    if cutoff_nm:
        system = forcefield.createSystem(
            pdb.topology, nonbondedMethod=CutoffNonPeriodic,
            nonbondedCutoff=cutoff_nm * unit.nanometer, constraints=HBonds)
        print(f"  {n_atoms:,} atoms, GBn2 implicit solvent, "
              f"{cutoff_nm} nm cutoff")
    else:
        system = forcefield.createSystem(
            pdb.topology, nonbondedMethod=NoCutoff, constraints=HBonds)
        print(f"  {n_atoms:,} atoms, GBn2 implicit solvent, no cutoff "
              "(O(N^2) - slow for large systems)")

    integrator = LangevinMiddleIntegrator(
        temperature * unit.kelvin, 1.0 / unit.picosecond,
        timestep_fs * unit.femtoseconds)
    integrator.setRandomNumberSeed(seed)

    try:
        platform = Platform.getPlatformByName("CPU")
        sim = Simulation(pdb.topology, system, integrator, platform)
    except Exception:
        sim = Simulation(pdb.topology, system, integrator)

    sim.context.setPositions(pdb.positions)
    print("  minimising...")
    sim.minimizeEnergy(maxIterations=1000)

    steps_per_ps = int(1000 / timestep_fs)
    equil_steps = int(equil_ps * steps_per_ps)
    prod_steps = int(ns * 1000 * steps_per_ps)
    report_every = max(1, int(report_ps * steps_per_ps))

    sim.context.setVelocitiesToTemperature(temperature * unit.kelvin, seed)
    print(f"  equilibrating {equil_ps:.0f} ps...")
    sim.step(equil_steps)

    traj = MD_DIR / f"{name}.dcd"
    log = MD_DIR / f"{name}_log.csv"
    sim.reporters.append(DCDReporter(str(traj), report_every))
    sim.reporters.append(StateDataReporter(
        str(log), report_every, step=True, time=True,
        potentialEnergy=True, temperature=True, speed=True))

    n_frames = prod_steps // report_every
    print(f"  production {ns*1000:.0f} ps ({prod_steps:,} steps, "
          f"~{n_frames} frames)...")
    t0 = time.time()
    sim.step(prod_steps)
    wall = time.time() - t0
    print(f"  done in {wall/60:.1f} min "
          f"({ns*1000/max(wall,1e-9)*60:.1f} ps/min)")

    final = MD_DIR / f"{name}_final.pdb"
    state = sim.context.getState(getPositions=True)
    with open(final, "w") as fh:
        PDBFile.writeFile(pdb.topology, state.getPositions(), fh, keepIds=True)

    return analyse(pdb_path, traj, site_resid, name, ns, wall)


def analyse(pdb_path, traj_path, site_resid, name, ns, wall):
    """Backbone RMSD, active-site RMSF and pocket width from the trajectory."""
    try:
        import MDAnalysis as mda
        from MDAnalysis.analysis import rms
    except ImportError:
        print("  (MDAnalysis not installed - skipping trajectory analysis)",
              file=sys.stderr)
        return {"system": name, "ns": ns, "wall_min": round(wall / 60, 2)}

    u = mda.Universe(str(pdb_path), str(traj_path))
    n_frames = len(u.trajectory)

    # Backbone RMSD against the minimised starting structure.
    R = rms.RMSD(u, u, select="backbone", ref_frame=0)
    R.run()
    rmsd = R.results.rmsd[:, 2]

    metrics = {
        "system": name,
        "ns": ns,
        "frames": n_frames,
        "wall_min": round(wall / 60, 2),
        "rmsd_mean": round(float(rmsd.mean()), 3),
        "rmsd_final": round(float(rmsd[-1]), 3),
        "rmsd_max": round(float(rmsd.max()), 3),
    }

    # Active-site flexibility: RMSF of residues near the specificity position.
    if site_resid is not None:
        try:
            site = u.select_atoms(f"protein and around 8 (resid {site_resid})")
            resids = sorted(set(site.residues.resids))
            sel = u.select_atoms(
                "name CA and resid " + " ".join(str(r) for r in resids))
            coords = np.array([sel.positions.copy() for _ in u.trajectory])
            mean = coords.mean(axis=0)
            rmsf = np.sqrt(((coords - mean) ** 2).sum(axis=2).mean(axis=0))
            metrics["site_residues"] = len(resids)
            metrics["site_rmsf_mean"] = round(float(rmsf.mean()), 3)
            metrics["site_rmsf_max"] = round(float(rmsf.max()), 3)

            # Pocket width proxy: mean pairwise CA spread across the site.
            widths = []
            for ts in u.trajectory:
                p = sel.positions
                widths.append(float(np.linalg.norm(
                    p - p.mean(axis=0), axis=1).mean()))
            metrics["pocket_radius_mean"] = round(float(np.mean(widths)), 3)
            metrics["pocket_radius_sd"] = round(float(np.std(widths)), 3)
        except Exception as exc:
            print(f"  (site analysis failed: {exc})", file=sys.stderr)

    print(f"  backbone RMSD {metrics['rmsd_mean']:.2f} A (max "
          f"{metrics['rmsd_max']:.2f})")
    if "site_rmsf_mean" in metrics:
        print(f"  active-site RMSF {metrics['site_rmsf_mean']:.2f} A, "
              f"pocket radius {metrics.get('pocket_radius_mean', float('nan')):.2f} A")
    return metrics


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--structure", action="append", default=None)
    ap.add_argument("--all-mutants", action="store_true")
    ap.add_argument("--ns", type=float, default=0.5,
                    help="production length in nanoseconds")
    ap.add_argument("--site", type=int, default=None,
                    help="residue number defining the active site")
    ap.add_argument("--equil-ps", type=float, default=20.0)
    ap.add_argument("--cutoff", type=float, default=2.0,
                    help="nonbonded cutoff in nm (0 = no cutoff)")
    args = ap.parse_args()

    if args.all_mutants:
        targets = sorted(STRUCT_DIR.glob("*_clean.pdb"))
    elif args.structure:
        targets = [Path(s) for s in args.structure]
    else:
        raise SystemExit("give --structure or --all-mutants")
    targets = [t for t in targets if t.exists()]
    if not targets:
        raise SystemExit("no structures found")

    print(f"{len(targets)} system(s), {args.ns} ns each, implicit solvent (GBn2)")
    results = []
    for t in targets:
        try:
            results.append(run_md(t, args.ns, args.site, equil_ps=args.equil_ps,
                                  cutoff_nm=args.cutoff or None))
        except Exception as exc:
            print(f"  {t.name}: FAILED - {type(exc).__name__}: {exc}",
                  file=sys.stderr)

    if results:
        out = MD_DIR / "md_summary.json"
        out.write_text(json.dumps(results, indent=2))
        print(f"\n-> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
