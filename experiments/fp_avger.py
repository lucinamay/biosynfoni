from sys import argv
import numpy as np
import matplotlib.pyplot as plt

lines = []
with open(argv[1], "r") as f:
    for line in f:
        lines.append(line.strip().split(","))

arr = np.array(lines, dtype=int)
mean_fp = arr.mean(axis=1)
print(mean_fp)

with open("fp_avg.csv", "w") as f:
    for i in range(len(mean_fp)):
        f.write(f"{mean_fp[i]}\t")

plt.ioff()
dataset = np.transpose(arr)
plt.violinplot(dataset=dataset, showmeans=True)
plt.savefig("fp_avg.svg")
