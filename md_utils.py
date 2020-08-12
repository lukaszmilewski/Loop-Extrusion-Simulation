import csv

import numpy as np
from matplotlib import pyplot as plt

def gen_line_array (N,F):
    a=np.arange(1,N+1)
    b=F*a/N
    c=b[::-1]
def gen_sin_array (N,F):
    a=np.arange(1,N+1)
    b=F*np.sin((a/(N/(np.pi/2))))
    c=b[::-1]
   
    return np.vstack((a,b,c))
def plot_data(state_fname, plot_fname):
    with open(state_fname) as csv_file:
        csv_reader = csv.DictReader(csv_file, fieldnames=['step', 'energy'], delimiter=",")
        next(csv_reader, None)  # skip headers
        steps = []
        energies = []
        n = 0
        for row in csv_reader:
            steps.append(int(row['step']))
            energies.append(float(row['energy']))
            n += 1
        steps = np.array(steps)
        energies = np.array(energies)
        fig, ax1 = plt.subplots()

        ax1.plot(steps, energies, lw=1, label="Potential energy")
        ax1.set_xlabel("Step")
        ax1.set_ylabel("Potential energy [kJ/mole]")

        plt.savefig(plot_fname)
        print("Plot saved in {} file".format(plot_fname))
