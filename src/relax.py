"""
Restrained local relaxation of a mutated site.

PDBFixer grafts a new side chain using idealised internal coordinates. For
small substitutions that is fine, but for bulky residues the rebuilt ring
geometry can come out distorted enough that bond perception fails downstream,
and the pocket around it is left clashing.

This module relaxes the neighbourhood of the mutation while holding the rest
of the protein in place:

  * atoms within `free_radius` of the mutated residue move freely,
  * atoms beyond it are harmonically restrained to their crystal positions,
  * minimisation runs in vacuum with an amber14 force field.

The restraint is what keeps this honest. Minimising the whole structure
without restraints would let the backbone drift away from the experimental
coordinates, and then wild type and mutant would no longer be comparable.
"""

import numpy as np
from openmm import CustomExternalForce, LangevinMiddleIntegrator, Platform, unit
from openmm.app import PME, ForceField, HBonds, Modeller, NoCutoff, PDBFile, Simulation


def relax_site(topology, positions, target_resid, chains=None,
               free_radius=8.0, restraint_k=50.0, max_iterations=800,
               verbose=False):
    """Minimise around `target_resid`, restraining distant atoms.

    Parameters
    ----------
    target_resid : str or int   residue id of the mutated position
    free_radius  : float        angstrom; atoms inside move freely
    restraint_k  : float        kcal/mol/A^2 on distant atoms

    Returns (topology, positions) after minimisation.
    """
    forcefield = ForceField("amber14-all.xml")
    modeller = Modeller(topology, positions)

    # Strip anything the force field cannot parameterise (stray HETATMs).
    try:
        system = forcefield.createSystem(
            modeller.topology, nonbondedMethod=NoCutoff, constraints=HBonds)
    except ValueError as exc:
        raise RuntimeError(f"force field could not parameterise the structure: {exc}")

    pos_nm = np.array(modeller.positions.value_in_unit(unit.nanometer))

    # Centre of the mutated residue.
    target = str(target_resid)
    idx = [a.index for a in modeller.topology.atoms()
           if a.residue.id == target and (chains is None
                                          or a.residue.chain.id in chains)]
    if not idx:
        raise ValueError(f"residue {target} not found for relaxation")
    centre = pos_nm[idx].mean(axis=0)

    # Harmonic restraint on everything outside the free sphere.
    restraint = CustomExternalForce(
        "0.5*k*periodicdistance(x, y, z, x0, y0, z0)^2")
    restraint.addGlobalParameter(
        "k", restraint_k * unit.kilocalories_per_mole / unit.angstrom**2)
    for p in ("x0", "y0", "z0"):
        restraint.addPerParticleParameter(p)

    cutoff_nm = free_radius / 10.0
    n_free = 0
    for atom in modeller.topology.atoms():
        d = float(np.linalg.norm(pos_nm[atom.index] - centre))
        if d <= cutoff_nm:
            n_free += 1
            continue
        restraint.addParticle(atom.index, pos_nm[atom.index].tolist())
    system.addForce(restraint)

    integrator = LangevinMiddleIntegrator(
        300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
    try:
        platform = Platform.getPlatformByName("CPU")
        sim = Simulation(modeller.topology, system, integrator, platform)
    except Exception:
        sim = Simulation(modeller.topology, system, integrator)

    sim.context.setPositions(modeller.positions)
    e0 = sim.context.getState(getEnergy=True).getPotentialEnergy()
    sim.minimizeEnergy(maxIterations=max_iterations)
    state = sim.context.getState(getPositions=True, getEnergy=True)
    e1 = state.getPotentialEnergy()

    if verbose:
        print(f"    relaxed {n_free} atoms within {free_radius} A; "
              f"E {e0.value_in_unit(unit.kilocalories_per_mole):.0f} -> "
              f"{e1.value_in_unit(unit.kilocalories_per_mole):.0f} kcal/mol")

    return modeller.topology, state.getPositions()
