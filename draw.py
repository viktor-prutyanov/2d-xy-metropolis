#!/usr/bin/python3

from PIL import Image
import numpy as np
import sys
import csv
import colorsys
import math

ls = []

with open(sys.argv[1]) as f:
    reader = csv.reader(f, delimiter=' ')
    for row in reader:
        ls.append(list(map(lambda x: float(x) / (2 * math.pi), row)))

a = np.array(ls, dtype=np.float)
r = np.vectorize(lambda x: np.uint8(colorsys.hsv_to_rgb(x, 1.0, 1.0)[0] * 255))(a)
g = np.vectorize(lambda x: np.uint8(colorsys.hsv_to_rgb(x, 1.0, 1.0)[1] * 255))(a)
b = np.vectorize(lambda x: np.uint8(colorsys.hsv_to_rgb(x, 1.0, 1.0)[2] * 255))(a)
pic = np.dstack((r, g, b))

v = np.zeros_like(a)

def d(x, y):
    return min((x - y) % 256, (y - x) % 256)

w = np.empty((2, 2), dtype=np.uint8)

for i in range(a.shape[0]):
    for j in range(a.shape[1]):
        for k in range(2):
            for l in range(2):
                w[k, l] = a[(i + k) % a.shape[0], (j + l) % a.shape[1]] * 255
        s = d(w[1, 0], w[0,0]) + d(w[1, 1], w[1, 0]) + d(w[0, 1], w[1, 1]) + d(w[0, 1], w[0, 0])
        if abs(s - 256) == 0:
            #print("w =\n", w)
            v[i, j] = 1

for i in range(a.shape[0]):
    for j in range(a.shape[1]):
        if v[i, j] == 1:
            for k in range(2):
                for l in range(2):
                    pic[(i + k) % a.shape[0], (j + l) % a.shape[1], :] = pic[(i + k) % a.shape[0], (j + l) % a.shape[1], :] / 2

im = Image.fromarray(pic)
im.save(sys.argv[1] + '.bmp')

nr_vortices = len(v[v == 1])
print(nr_vortices, nr_vortices / (a.shape[0] * a.shape[1]))
