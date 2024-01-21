import sys, os
import argparse

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
from biosynfoni.subkeys import fpVersions, defaultVersion, get_names
from utils.figuremaking import heatmap, annotate_heatmap
from utils import set_style


def cli():
    parser = argparse.ArgumentParser()
    # parser.add_argument("fingerprintfile_coco", type=str)
    # parser.add_argument("fingerprintfile_zinc", type=str)
    # parser.add_argument("bsf_name", type=str, default=defaultVersion)
    parser.add_argument(
        "fingerprints", type=str, help="path to tsv or csv file of fingerprints"
    )
    parser.add_argument(
        "-n",
        "--name",
        type=str,
        help="name of compound collection for title and filename",
    )
    args = parser.parse_args()
    return args


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


def fp_plots(fps, bsf_name):
    plt.ioff()
    plt.figure().set_figwidth(15)
    print("making plot")
    plt.violinplot(dataset=fps, showmeans=True)
    print("saving plot")
    plt.savefig(f"fp_avg_{bsf_name}.svg")
    plt.close()
    return None


def heatmap_array(
    fps: np.array,
    max_height: int = 30,
    percentages=False,
    accumulative=True,
    end_accumulative=False,
):
    heat_array = np.zeros((max_height, fps.shape[1]))
    for i in range(max_height):
        if accumulative:
            countrow = np.count_nonzero(fps > i, axis=0)
        else:
            if end_accumulative and i == max_height - 1:
                # for last height, count all remaining values
                countrow = np.count_nonzero(fps > i, axis=0)
            else:
                countrow = np.count_nonzero((fps == i + 1), axis=0)
        heat_array[max_height - 1 - i] = countrow

    if percentages:
        heat_array = heat_array / fps.shape[0] * 100
    return heat_array.astype(int)


def fp_heatmap(
    fp_hm_array: np.array,
    subslabels: list = [],
    size: tuple[int] = (10, 6),
    percentages: bool = False,
    annotate: bool = False,
    color_scheme: str = "Purples",
    title: str = "Representative substructure count for compound collection",
    top_acc_array=None,
):
    cbarlab = "number of compounds"
    if percentages:
        cbarlab = "percentage of compounds"

    print("saving heatmap")
    height = fp_hm_array.shape[0]
    fig, ax = plt.subplots(figsize=size, dpi=500)
    if not subslabels:
        subslabels = [f"subs{i}" for i in range(1, fp_hm_array.shape[1] + 1)]

    yaxlabels = [(height + 1 - i) for i in range(1, height + 1)]
    if top_acc_array is not None:
        yaxlabels[0] = f"≥{height}"
        im2, cbar2 = heatmap(
            top_acc_array,
            # ['>11']+[(height+1-i) for i in range(1, height + 1)],
            yaxlabels,
            subslabels,
            ax=ax,
            cmap="Greys",
            # cmap = "PiYG",
            cbarlabel=f"{cbarlab}, accumulative",
            cbar_kw={"drawedges": False, "shrink": 0.5, "pad": 0.05, "aspect": 10},
        )

    im, cbar = heatmap(
        fp_hm_array,
        # [(height+1-i) for i in range(1, height + 1)],
        yaxlabels,
        subslabels,
        ax=ax,
        cmap=color_scheme,
        # cmap = "PiYG",
        cbarlabel=cbarlab,
        cbar_kw={"drawedges": False, "shrink": 0.5, "pad": 0.05, "aspect": 10},
    )

    # texts = annotate_heatmap(im, valfmt="{x:.1f}")
    if annotate:
        texts = annotate_heatmap(im, valfmt="{x:.0f}", size=7)
    # plt.figure(figsize=(10,6))
    ax.set_xlabel("substructure", labelpad=10)
    ax.set_ylabel("counts", labelpad=10)
    ax.set_title(title, loc="center", pad=20)
    fig.tight_layout()
    # plt.show()
    # for filename, using collection name is nice
    # plt.savefig(f"{outfile_namer(name)}_heatmap.png")
    return None


def over_under_divide(fps: np.array, limit: int = 10, percentages: bool = True):
    full = heatmap_array(
        fps,
        max_height=limit + 1,
        percentages=percentages,
        accumulative=False,
        end_accumulative=True,
    )
    under, over = full.astype(float).copy(), full.astype(float).copy()
    under[0] = np.nan
    over[1:] = np.nan
    return under, over


def fp_heatmap_accumulative(fp_arr: np.array, limit: int = 10, *args, **kwargs):
    under, over = over_under_divide(fp_arr, limit, percentages=True)
    fp_heatmap(
        under,
        *args,
        percentages=True,
        top_acc_array=over,
        **kwargs,
    )
    return None


def main():
    print("hello")
    set_style()
    args = cli()
    fps = np.loadtxt(args.fingerprints, dtype=int, delimiter=",")
    if args.name:
        fp_name = args.name
    else:
        fp_name = args.fingerprints.split("/")[-1].split(".")[0]
    # coco_name = fingerprintfile_coco.split("/")[-1].split(".")[0]
    # bsf_name_coco = coco_name.replace("_noblock", '').replace(split("_")[-1]

    substructure_names = get_names(version=defaultVersion)
    if fps.shape[1] != len(substructure_names):
        substructure_names = [f"{i}" for i in range(1, fps.shape[1] + 1)]

    # if all values in fps are 0 or 1, then it's a binary fingerprint
    if np.all(np.isin(fps, [0, 1])):
        fp_heatmap(
            heatmap_array(fps, max_height=1, percentages=True, accumulative=False),
            subslabels=substructure_names,
            title=f"Distribution of {fp_name} substructure counts",
            color_scheme="Greys",
            percentages=True,
            size=(15, 1),
        )
        plt.savefig(f"{outfile_namer(fp_name)}_heatmap.png")
        return None

    # coco_mean = fp_stats(coco, coco_name)
    # zinc_mean = fp_stats(zinc, zinc_name)
    # fp_plots(coco, coco_name)
    # fp_plots(zinc, zinc_name)

    fp_heatmap_accumulative(
        fps,
        limit=10,
        title=f"Distribution of {fp_name} substructure counts",
        subslabels=substructure_names,
        color_scheme="GnBu",
    )
    plt.savefig(f"{outfile_namer(fp_name)}_heatmap.png")

    # fp_means_plots(coco_mean, zinc_mean, outfile_namer(f"{coco_name}_{zinc_name}.svg"))


if __name__ == "__main__":
    main()
