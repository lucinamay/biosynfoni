from sys import argv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

matplotlib.use("Agg")
plt.ioff()


def fp_stats(fp_arr, set_name):
    mean_fp = fp_arr.mean(axis=0)
    outfile = f"outfile_namer('fp_avg', set_name).csv"
    with open(outfile, "w") as f:
        f.write("\t".join([str(x) for x in mean_fp]))
        f.write("\n")
        f.write("\t".join([str(x) for x in std_fp]))


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
    plt.savefig(f"mean_fp_both.svg")
    plt.close()
    return None


def fp_plots(fp_arr, set_name):
    plt.ioff()
    plt.figure().set_figwidth(15)
    print("making plot")
    plt.violinplot(dataset=arr, showmeans=True)
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

    fp_stats(coco, coco_name)
    fp_stats(zinc, zinc_name)
    fp_plots(coco, coco_name)
    fp_plots(zinc, zinc_name)
