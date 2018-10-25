#!/usr/bin/python3

import csv
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sys

ylabels = [r'Energy $E$', r'Magnetization $m$', r'Heat capacity $C_v$', r'Magnetic susceptibility $\chi$', r"Binder's 4th order cumulant $U_L$", r'Vortex density $\rho_v$']

labels = [16, 20, 32]

for i in range(1, len(sys.argv)):
    ls = []
    with open(sys.argv[i]) as csv_file:
        reader = csv.reader(csv_file, delimiter=',')
        for row in reader:
            ls.append(list(map(float, row)))
        a = np.array(ls, dtype=float)
        N = a.shape[1]
        xs = a[:,0]

        for j in range(1, N):
            plt.subplot(N // 3, 3, j)
            ys = a[:,j]
            plt.xlabel(r'Temperature $T$')
            plt.ylabel(ylabels[j - 1])
            if j == 3:
                ys = [y * (labels[i - 1] ** 2) for y in ys] # Hint to fix Cv wrong scale
            plt.plot(xs, ys, '.-', label=str(labels[i - 1]))
            if j == 1:
                plt.legend()
            plt.grid()

plt.show()
