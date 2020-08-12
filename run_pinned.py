from sys import stdout

import simtk.openmm as mm
import simtk.unit as u
from simtk.openmm.app import PDBFile, ForceField, Simulation, DCDReporter, StateDataReporter
from md_utils import plot_data
import numpy as np
from md_utils import gen_sin_array
from md_utils import gen_line_array

STATE_FNAME = 'state.csv'
STEPS = 10000
LE_FORCE_SCALE = 3 * u.kilocalories_per_mole / u.angstroms ** 2
STEPS_PER_CYCLE = 200

#Macierz z parametrami sił wiązań
LE_FORCE_MATRIX = np.array([[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201],
[0.000001,0.015,0.030,0.045,0.060,0.075,0.090,0.105,0.120,0.135,0.150,0.165,0.180,0.195,0.210,0.225,0.240,0.255,0.270,0.285,0.300,0.315,0.330,0.345,0.360,0.375,0.390,0.405,0.420,0.435,0.450,0.465,0.480,0.495,0.510,0.525,0.540,0.555,0.570,0.585,0.600,0.615,0.630,0.645,0.660,0.675,0.690,0.705,0.720,0.735,0.750,0.765,0.780,0.795,0.810,0.825,0.840,0.855,0.870,0.885,0.900,0.915,0.930,0.945,0.960,0.975,0.990,1.005,1.020,1.035,1.050,1.065,1.080,1.095,1.110,1.125,1.140,1.155,1.170,1.185,1.200,1.215,1.230,1.245,1.260,1.275,1.290,1.305,1.320,1.335,1.350,1.365,1.380,1.395,1.410,1.425,1.440,1.455,1.470,1.485,1.500,1.515,1.530,1.545,1.560,1.575,1.590,1.605,1.620,1.635,1.650,1.665,1.680,1.695,1.710,1.725,1.740,1.755,1.770,1.785,1.800,1.815,1.830,1.845,1.860,1.875,1.890,1.905,1.920,1.935,1.950,1.965,1.980,1.995,2.010,2.025,2.040,2.055,2.070,2.085,2.100,2.115,2.130,2.145,2.160,2.175,2.190,2.205,2.220,2.235,2.250,2.265,2.280,2.295,2.310,2.325,2.340,2.355,2.370,2.385,2.400,2.415,2.430,2.445,2.460,2.475,2.490,2.505,2.520,2.535,2.550,2.565,2.580,2.595,2.610,2.625,2.640,2.655,2.670,2.685,2.700,2.715,2.730,2.745,2.760,2.775,2.790,2.805,2.820,2.835,2.850,2.865,2.880,2.895,2.910,2.925,2.940,2.955,2.970,2.985,3.000],
[3.000,2.985,2.970,2.955,2.940,2.925,2.910,2.895,2.880,2.865,2.850,2.835,2.820,2.805,2.790,2.775,2.760,2.745,2.730,2.715,2.700,2.685,2.670,2.655,2.640,2.625,2.610,2.595,2.580,2.565,2.550,2.535,2.520,2.505,2.490,2.475,2.460,2.445,2.430,2.415,2.400,2.385,2.370,2.355,2.340,2.325,2.310,2.295,2.280,2.265,2.250,2.235,2.220,2.205,2.190,2.175,2.160,2.145,2.130,2.115,2.100,2.085,2.070,2.055,2.040,2.025,2.010,1.995,1.980,1.965,1.950,1.935,1.920,1.905,1.890,1.875,1.860,1.845,1.830,1.815,1.800,1.785,1.770,1.755,1.740,1.725,1.710,1.695,1.680,1.665,1.650,1.635,1.620,1.605,1.590,1.575,1.560,1.545,1.530,1.515,1.500,1.485,1.470,1.455,1.440,1.425,1.410,1.395,1.380,1.365,1.350,1.335,1.320,1.305,1.290,1.275,1.260,1.245,1.230,1.215,1.200,1.185,1.170,1.155,1.140,1.125,1.110,1.095,1.080,1.065,1.050,1.035,1.020,1.005,0.990,0.975,0.960,0.945,0.930,0.915,0.900,0.885,0.870,0.855,0.840,0.825,0.810,0.795,0.780,0.765,0.750,0.735,0.720,0.705,0.690,0.675,0.660,0.645,0.630,0.615,0.600,0.585,0.570,0.555,0.540,0.525,0.510,0.495,0.480,0.465,0.450,0.435,0.420,0.405,0.390,0.375,0.360,0.345,0.330,0.315,0.300,0.285,0.270,0.255,0.240,0.225,0.210,0.195,0.180,0.165,0.150,0.135,0.120,0.105,0.090,0.075,0.060,0.045,0.030,0.015,0.000001]])

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
