from sys import stdout

import simtk.openmm as mm
import simtk.unit as u
from simtk.openmm.app import PDBFile, ForceField, Simulation, DCDReporter, StateDataReporter
from md_utils import plot_data
from md_utils import gen_sin_array
import numpy as np


STATE_FNAME = 'state.csv'
STEPS = 10000
LE_FORCE_SCALE = 3 * u.kilocalories_per_mole / u.angstroms ** 2
STEPS_PER_CYCLE = 200

#Macierz z parametrami sił wiązań
#Dodano funkcje generacji macierzy o wartościach sinusoidalnych. Funkcja ta przyjmuje dwa argumenty. Pierwszy oznacza liczbę kroków które ma posiadać macierz a drugi
#stanowi regulacje maksymalnej siły (tzn jeśli wstawimy 3 to maksymalna siła bedzie tyle wynosić)
LE_FORCE_MATRIX = gen_sin_array(200,3)

pdb = PDBFile('initial_structure.pdb')
forcefield = ForceField('polymer_ff.xml')
system = forcefield.createSystem(pdb.topology, nonbondedCutoff=1 * u.nanometer)
integrator = mm.LangevinIntegrator(100 * u.kelvin, 0.2, 1 * u.femtoseconds)

# Distance constraint
for i in range(system.getNumParticles() - 1):
    system.addConstraint(i, i + 1, 0.1 * u.nanometer)

# Pinning ends with rubber
pin_force = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
pin_force.addGlobalParameter("k", 50 * u.kilocalories_per_mole / u.angstroms ** 2)
pin_force.addPerParticleParameter("x0")
pin_force.addPerParticleParameter("y0")
pin_force.addPerParticleParameter("z0")
pin_force.addParticle(0, [15 * u.angstrom, 0 * u.angstrom, 0 * u.angstrom])
pin_force.addParticle(system.getNumParticles() - 1, [-15 * u.angstrom, 0 * u.angstrom, 0 * u.angstrom])
system.addForce(pin_force)

# Loop extrusion force
le_force = mm.HarmonicBondForce()
le_force.addBond(48, 50, 1 * u.angstrom, LE_FORCE_MATRIX[2][0] * u.kilocalories_per_mole / u.angstroms ** 2)
for i in range(2, 35):
    p1, p2 = 49 - i, 49 + i
    le_force.addBond(p1, p2, 1 * u.angstrom, LE_FORCE_MATRIX[1][0] * u.kilocalories_per_mole / u.angstroms ** 2)
system.addForce(le_force)

simulation = Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)
simulation.minimizeEnergy()
simulation.reporters.append(DCDReporter('trj.dcd', 1))
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True))
simulation.reporters.append(StateDataReporter(STATE_FNAME, 10, step=True, potentialEnergy=True))

simulation.step(STEPS_PER_CYCLE)

for i in range(2, 35):
    p1, p2 = 49 - i, 49 + i
    for j in range(0, STEPS_PER_CYCLE):
        le_force.setBondParameters(i - 2, p1 + 1, p2 - 1, 1 * u.angstrom,
                                  LE_FORCE_MATRIX[1][j] * u.kilocalories_per_mole / u.angstroms ** 2)
        le_force.setBondParameters(i - 1, p1, p2, 1 * u.angstrom, LE_FORCE_MATRIX[2][j])
        le_force.updateParametersInContext(simulation.context)
        simulation.step(1)
    # le_force.setBondParameters(i - 2, p1 + 1, p2 - 1, 1 * u.angstrom,
    #                           LE_FORCE_MATRIX[1][0] * u.kilocalories_per_mole / u.angstroms ** 2)
    # le_force.setBondParameters(i - 1, p1, p2, 1 * u.angstrom, LE_FORCE_MATRIX[2][0])
    # simulation.minimizeEnergy()
    # le_force.updateParametersInContext(simulation.context)
    # simulation.step(STEPS_PER_CYCLE)

    plot_data(STATE_FNAME, 'energy.pdf')

print('#1: repr stick; color white; color red :1,100; repr sphere :1,100; vdwdefine 0.5')
print('#1: color green :49,51; repr sphere :49,51; color #ffffa2e8a2e8 :50;')
for i in range(1, 35):
    p1, p2 = 50 - i - 1, 50 + i + 1
    print(
        f'#{i*STEPS_PER_CYCLE+1}: color green :{p1},{p2}; repr sphere :{p1},{p2}; repr stick :{p1+1},{p2-1}; color #ffffa2e8a2e8 :{p1+1}-{p2-1};')
