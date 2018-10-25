#!/usr/bin/python3

import csv
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sys

for i in range(1, len(sys.argv)):
    ls = []
    with open(sys.argv[i]) as csv_file:
        reader = csv.reader(csv_file, delimiter=',')
        for row in reader:
            ls.append(list(map(float, row)))
        a = np.array(ls, dtype=float)
        N = a.shape[1]
        xs = a[:,0]

        for i in range(1, N):
            plt.subplot(N, 1, i)
            ys = a[:,i]
            plt.xlabel(r"$N$ steps")
            plt.ylabel(r"Energy")
            plt.loglog(xs, -ys, '.-')
            plt.grid()

plt.show()
