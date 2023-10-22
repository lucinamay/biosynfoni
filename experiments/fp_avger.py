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
from biosynfoni.def_biosynfoni import FP_VERSIONS
from utils.figuremaking import heatmap, annotate_heatmap


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


def heatmap_array(fp_arr: np.array, max_height: int = 30):
    heat_array = np.zeros((max_height, fp_arr.shape[1]))
    for i in range(max_height):
        countrow = np.count_nonzero(fp_arr > i, axis=0)
        heat_array[max_height - 1 - i] = countrow
    return heat_array.astype(int)


def fp_heatmap(
    fp_arr: np.array,
    filename: str,
    subsnames: list = [],
    bsfname: str = "",
):
    height = fp_arr.shape[0]
    fig, ax = plt.subplots()
    if not subsnames:
        subsnames = [f"subs{i}" for i in range(1, fp_arr.shape[1] + 1)]
    im, cbar = heatmap(
        fp_arr,
        [i for i in range(1, height + 1)],
        subsnames,
        ax=ax,
        cmap="YlGn",
        cbarlabel="number of compounds",
    )
    # texts = annotate_heatmap(im, valfmt="{x:.1f}")
    texts = annotate_heatmap(im, valfmt="{x:.0e)")
    # plt.figure(figsize=(10,6))
    fig.tight_layout()
    # plt.show()
    # for filename, using collection name is nice
    plt.savefig(f"{outfile_namer(filename, bsfname)}_heatmap.sgv")
    return None


def main():
    print("hello")
    fingerprintfile_coco = argv[1]  # natural products
    fingerprintfile_zinc = argv[2]  # synthetic compounds
    coco = np.loadtxt(fingerprintfile_coco, dtype=int, delimiter=",")
    zinc = np.loadtxt(fingerprintfile_zinc, dtype=int, delimiter=",")
    coco_name = fingerprintfile_coco.split("/")[-1].split(".")[0]
    zinc_name = fingerprintfile_zinc.split("/")[-1].split(".")[0]
    # bsf_name_coco = coco_name.replace("_noblock", '').replace(split("_")[-1]
    if len(argv) > 3:
        bsf_name = argv[3]
    if FP_VERSIONS[bsf_name]:
        substructure_names = FP_VERSIONS[bsf_name]

    coco_mean = fp_stats(coco, coco_name)
    zinc_mean = fp_stats(zinc, zinc_name)
    fp_plots(coco, coco_name)
    fp_plots(zinc, zinc_name)
    fp_heatmap(
        heatmap_array(coco, max_height=30),
        coco_name,
        subsnames=substructure_names,
        bsfname=bsf_name,
    )
    fp_heatmap(
        heatmap_array(zinc, max_height=30),
        zinc_name,
        subsnames=substructure_names,
        bsfname=bsf_name,
    )

    fp_means_plots(coco_mean, zinc_mean, outfile_namer(f"{coco_name}_{zinc_name}.svg"))


if __name__ == "__main__":
    main()
