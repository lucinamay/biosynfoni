import sys, os
from sys import argv

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.ioff()

sys.path.append(os.path.abspath(os.path.join(sys.path[0], os.pardir, "src")))
# for intra-biosynfoni-code running
sys.path.append(
    os.path.abspath(os.path.join(sys.path[0], os.pardir, "src", "biosynfoni"))
)
from biosynfoni.inoutput import outfile_namer


def fp_stats(fp_arr: np.array, set_name: str) -> np.array:
    mean_fp = fp_arr.mean(axis=0)
    std_fp = fp_arr.std(axis=0)
    outfile = f"{outfile_namer('fp_avg', set_name)}.csv"
    with open(outfile, "w") as f:
        f.write("\t".join([str(x) for x in mean_fp]))
        f.write("\n")
        f.write("\t".join([str(x) for x in std_fp]))
    return mean_fp


def fp_means_plots(fp_mean_arr1, fp_mean_arr2, filename):
    plt.figure()
    print("making plot")

    joined = np.array([fp_mean_arr1, fp_mean_arr2])
    both = np.transpose(np.array([fp_mean_arr1, fp_mean_arr2]))

    plt.stairs(fp_mean_arr1, color="green", label="natural products")
    plt.stairs(fp_mean_arr2, color="grey", label="synthetic compounds")

    plt.xlabel("substructure no.")  # , fontsize=20)
    plt.ylabel("mean fp count")  # , fontsize=20)
    plt.title(
        "mean fp count per substructure for natural products and synthetic compounds"
    )  # , fontsize=20)
    print("saving plot")
    plt.legend()
    plt.savefig(filename)
    plt.close()
    return None


def fp_plots(fp_arr, bsf_name):
    plt.ioff()
    plt.figure().set_figwidth(15)
    print("making plot")
    plt.violinplot(dataset=fp_arr, showmeans=True)
    print("saving plot")
    plt.savefig(f"fp_avg_{bsf_name}.svg")
    plt.close()
    return None


def main():
    print("hello")
    fingerprintfile_coco = argv[1]  # natural products
    fingerprintfile_zinc = argv[2]  # synthetic compounds
    coco = np.loadtxt(fingerprintfile_coco, dtype=int, delimiter=",")
    zinc = np.loadtxt(fingerprintfile_zinc, dtype=int, delimiter=",")
    coco_name = fingerprintfile_coco.split("/")[-1].split(".")[0]
    zinc_name = fingerprintfile_zinc.split("/")[-1].split(".")[0]

    coco_mean = fp_stats(coco, coco_name)
    zinc_mean = fp_stats(zinc, zinc_name)
    fp_plots(coco, coco_name)
    fp_plots(zinc, zinc_name)

    fp_means_plots(coco_mean, zinc_mean, outfile_namer(f"{coco_name}_{zinc_name}.svg"))


if __name__ == "__main__":
    main()
