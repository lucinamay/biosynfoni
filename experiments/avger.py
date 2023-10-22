from sys import argv
import numpy as np

lines = np.loadtxt(argv[1], dtype=float, delimiter=",")

print(np.sum(lines) / len(lines))
