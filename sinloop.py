from sys import stdout

import simtk.openmm as mm
import simtk.unit as u
from simtk.openmm.app import PDBFile, ForceField, Simulation, DCDReporter, StateDataReporter
from md_utils import plot_data
import numpy as np


STATE_FNAME = 'state.csv'
STEPS = 10000
LE_FORCE_SCALE = 3 * u.kilocalories_per_mole / u.angstroms ** 2
STEPS_PER_CYCLE = 200

#Macierz z parametrami sił wiązań
LE_FORCE_MATRIX =np.array([ [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200],[0,0.023561703,0.047121952,0.070679295,0.094232277,0.117779447,0.141319352,0.16485054,0.188371559,0.211880958,0.235377287,0.258859097,0.28232494,0.305773367,0.329202933,0.352612192,0.375999701,0.399364016,0.422703696,0.446017302,0.469303395,0.49256054,0.515787301,0.538982246,0.562143944,0.585270966,0.608361886,0.631415279,0.654429724,0.677403801,0.700336092,0.723225182,0.746069661,0.768868119,0.79161915,0.81432135,0.836973318,0.859573658,0.882120976,0.90461388,0.927050983,0.949430902,0.971752255,0.994013665,1.016213761,1.038351171,1.060424531,1.082432479,1.104373658,1.126246714,1.148050297,1.169783063,1.191443672,1.213030787,1.234543076,1.255979213,1.277337875,1.298617745,1.31981751,1.340935862,1.361971499,1.382923124,1.403789443,1.42456917,1.445261022,1.465863724,1.486376005,1.506796599,1.527124247,1.547357695,1.567495694,1.587537003,1.607480385,1.62732461,1.647068454,1.666710699,1.686250134,1.705685552,1.725015756,1.744239553,1.763355757,1.782363189,1.801260676,1.820047053,1.838721161,1.857281848,1.875727969,1.894058386,1.912271969,1.930367594,1.948344145,1.966200513,1.983935596,2.001548301,2.019037541,2.036402237,2.053641318,2.070753721,2.08773839,2.104594277,2.121320344,2.137915557,2.154378893,2.170709338,2.186905882,2.202967528,2.218893285,2.23468217,2.250333209,2.265845437,2.281217897,2.29644964,2.311539728,2.32648723,2.341291222,2.355950793,2.370465037,2.38483306,2.399053975,2.413126906,2.427050983,2.440825349,2.454449152,2.467921554,2.481241723,2.494408837,2.507422084,2.520280662,2.532983777,2.545530645,2.557920493,2.570152557,2.582226081,2.594140322,2.605894543,2.617488021,2.62892004,2.640189895,2.65129689,2.662240341,2.673019573,2.68363392,2.694082727,2.704365351,2.714481157,2.724429521,2.73420983,2.743821479,2.753263877,2.76253644,2.771638598,2.780569787,2.789329458,2.797917069,2.806332092,2.814574008,2.822642307,2.830536492,2.838256076,2.845800584,2.853169549,2.860362517,2.867379044,2.874218698,2.880881057,2.887365709,2.893672255,2.899800306,2.905749483,2.91151942,2.917109761,2.922520161,2.927750286,2.932799813,2.937668432,2.942355841,2.946861752,2.951185887,2.955327978,2.959287772,2.963065022,2.966659496,2.970070973,2.973299242,2.976344104,2.979205371,2.981882866,2.984376425,2.986685894,2.988811129,2.990752001,2.992508389,2.994080185,2.995467292,2.996669625,2.997687109,2.998519681,2.999167291,2.999629897,2.999907473,3],
[3,2.976438297,2.952878048,2.929320705,2.905767723,2.882220553,2.858680648,2.83514946,2.811628441,2.788119042,2.764622713,2.741140903,2.71767506,2.694226633,2.670797067,2.647387808,2.624000299,2.600635984,2.577296304,2.553982698,2.530696605,2.50743946,2.484212699,2.461017754,2.437856056,2.414729034,2.391638114,2.368584721,2.345570276,2.322596199,2.299663908,2.276774818,2.253930339,2.231131881,2.20838085,2.18567865,2.163026682,2.140426342,2.117879024,2.09538612,2.072949017,2.050569098,2.028247745,2.005986335,1.983786239,1.961648829,1.939575469,1.917567521,1.895626342,1.873753286,1.851949703,1.830216937,1.808556328,1.786969213,1.765456924,1.744020787,1.722662125,1.701382255,1.68018249,1.659064138,1.638028501,1.617076876,1.596210557,1.57543083,1.554738978,1.534136276,1.513623995,1.493203401,1.472875753,1.452642305,1.432504306,1.412462997,1.392519615,1.37267539,1.352931546,1.333289301,1.313749866,1.294314448,1.274984244,1.255760447,1.236644243,1.217636811,1.198739324,1.179952947,1.161278839,1.142718152,1.124272031,1.105941614,1.087728031,1.069632406,1.051655855,1.033799487,1.016064404,0.998451699,0.980962459,0.963597763,0.946358682,0.929246279,0.91226161,0.895405723,0.878679656,0.862084443,0.845621107,0.829290662,0.813094118,0.797032472,0.781106715,0.76531783,0.749666791,0.734154563,0.718782103,0.70355036,0.688460272,0.67351277,0.658708778,0.644049207,0.629534963,0.61516694,0.600946025,0.586873094,0.572949017,0.559174651,0.545550848,0.532078446,0.518758277,0.505591163,0.492577916,0.479719338,0.467016223,0.454469355,0.442079507,0.429847443,0.417773919,0.405859678,0.394105457,0.382511979,0.37107996,0.359810105,0.34870311,0.337759659,0.326980427,0.31636608,0.305917273,0.295634649,0.285518843,0.275570479,0.26579017,0.256178521,0.246736123,0.23746356,0.228361402,0.219430213,0.210670542,0.202082931,0.193667908,0.185425992,0.177357693,0.169463508,0.161743924,0.154199416,0.146830451,0.139637483,0.132620956,0.125781302,0.119118943,0.112634291,0.106327745,0.100199694,0.094250517,0.08848058,0.082890239,0.077479839,0.072249714,0.067200187,0.062331568,0.057644159,0.053138248,0.048814113,0.044672022,0.040712228,0.036934978,0.033340504,0.029929027,0.026700758,0.023655896,0.020794629,0.018117134,0.015623575,0.013314106,0.011188871,0.009247999,0.007491611,0.005919815,0.004532708,0.003330375,0.002312891,0.001480319,0.000832709,0.000370103,0.00009,0]])




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
