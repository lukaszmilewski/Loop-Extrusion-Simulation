from sinloop import loop_extrusion as twosided
from sinloop2 import loop_extrusion as onesided
from kontakty import heatmap

STEPS = 10000
LE_FORCE_SCALE = 3
MATRIX_LENGTH = 1000
STEPS_PER_CYCLE = 10
STEPS_PER_IT = 1

dcdfile2 = r'C:\Users\lukas\PycharmProjects\Loop-Extrusion-Simulation\2sided-trj.dcd'
dcdfile1 = r'C:\Users\lukas\PycharmProjects\Loop-Extrusion-Simulation\1sided-trj.dcd'
pdbfile = r'C:\Users\lukas\PycharmProjects\Loop-Extrusion-Simulation\initial_structure.pdb'
outputmap2 = "2sided-heatmap.png"
outputmap1 = "1sided-heatmap.png"

print("Simulating two-sided model...")
twosided(STEPS, LE_FORCE_SCALE, MATRIX_LENGTH, STEPS_PER_CYCLE, STEPS_PER_IT)
print("Making heatmap...")
heatmap(dcdfile2, pdbfile, outputmap2)
print("Done")

print("Simulating one-sided model...")
onesided(STEPS, LE_FORCE_SCALE, MATRIX_LENGTH, STEPS_PER_CYCLE, STEPS_PER_IT)
print("Making heatmap...")
heatmap(dcdfile1, pdbfile, outputmap1)
print("Done")