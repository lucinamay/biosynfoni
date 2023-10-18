from sys import argv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

matplotlib.use("Agg")
plt.ioff()

arr = np.loadtxt(argv[1], dtype=int, delimiter=",")
bsf_name = argv[1].split("/")[-1]

mean_fp = arr.mean(axis=1)
print(mean_fp)

with open("fp_avg.csv", "w") as f:
    for i in range(len(mean_fp)):
        f.write(f"{mean_fp[i]}\t")

plt.ioff()
plt.figure().set_figwidth(15)
dataset = np.transpose(arr)
print("making plot")
plt.violinplot(dataset=arr, showmeans=True)
print("saving plot")
plt.savefig(f"fp_avg_{bsf_name}.svg")
plt.close()
