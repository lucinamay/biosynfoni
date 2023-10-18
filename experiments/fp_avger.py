from sys import argv
import numpy as np
import matplotlib.pyplot as plt

lines = []
with open(argv[1], "r") as f:
    for line in f:
        lines.append(line.strip().split(","))

arr = np.array(lines, dtype=int)
print(arr.mean(axis=1))

with open("fp_avg.csv", "w") as f:
    for i in range(len(arr)):
        f.write(f"{i},{arr[i].mean()}\t")

plt.ioff()
dataset = np.transpose(arr)
plt.violinplot(dataset=dataset, showmeans=True)
plt.savefig("fp_avg.svg")
