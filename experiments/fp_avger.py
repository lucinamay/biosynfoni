import sys, os, logging
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
from biosynfoni.subkeys import fpVersions, defaultVersion, get_names, get_pathway
from utils.figuremaking import (
    heatmap,
    annotate_heatmap,
    savefig,
    set_label_colors_from_categories,
    custom_cmap,
)
from utils import set_style
from utils.colours import colourDict


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
    parser.add_argument(
        "-c",
        "--classes",
        type=str,
        required=False,
        help="path to tsv or csv file of classes, uses first column as class label. pass if you want to also plot classes separately",
    )
    args = parser.parse_args()
    args.fingerprints = os.path.abspath(args.fingerprints)
    args.classes = os.path.abspath(args.classes) if args.classes else None
    return args


def write_stats(fp_arr: np.array) -> np.array:
    mean_fp = fp_arr.mean(axis=0)
    std_fp = fp_arr.std(axis=0)
    median_fp = np.median(fp_arr, axis=0)
    with open("stats.tsv", "w") as f:
        f.write("mean\n")
        f.write("\t".join([str(x) for x in mean_fp]))
        f.write("\nstdev\n")
        f.write("\t".join([str(x) for x in std_fp]))
        f.write("\nmedian\n")
        f.write("\t".join([str(x) for x in median_fp]))
    return mean_fp


def fp_means_plots(fp_mean_arr1, fp_mean_arr2, filename):
    plt.figure()
    logging.info("making plot")

    joined = np.array([fp_mean_arr1, fp_mean_arr2])
    both = np.transpose(np.array([fp_mean_arr1, fp_mean_arr2]))

    plt.stairs(fp_mean_arr1, color="green", label="natural products")
    plt.stairs(fp_mean_arr2, color="grey", label="synthetic compounds")

    plt.xlabel("substructure no.")  # , fontsize=20)
    plt.ylabel("mean fp count")  # , fontsize=20)
    plt.title(
        "mean fp count per substructure for natural products and synthetic compounds"
    )  # , fontsize=20)
    logging.info("saving plot")
    plt.legend()
    plt.savefig(filename)
    plt.close()
    return None


def fp_plots(fps, bsf_name):
    plt.ioff()
    plt.figure().set_figwidth(15)
    logging.info("making plot")
    plt.violinplot(dataset=fps, showmeans=True)
    logging.info("saving plot")
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
    standard_colour: bool = False,
):
    cbarlab = "number of compounds"
    if percentages:
        cbarlab = "% of compounds"

    logging.info("saving heatmap")
    height = fp_hm_array.shape[0]
    fig, ax = plt.subplots(figsize=size, dpi=500)
    if not subslabels:
        subslabels = [f"subs{i}" for i in range(1, fp_hm_array.shape[1] + 1)]
    subslabels = [x.replace("_", " ") for x in subslabels]

    yaxlabels = [(height + 1 - i) for i in range(1, height + 1)]
    if top_acc_array is not None:
        yaxlabels[0] = f"≥{height}"
        maxtop = top_acc_array[~np.isnan(top_acc_array)].max()
        maxfp = fp_hm_array[~np.isnan(fp_hm_array)].max()
        maxval = max(maxtop, maxfp)
        im2, cbar2 = heatmap(
            top_acc_array,
            # ['>11']+[(height+1-i) for i in range(1, height + 1)],
            yaxlabels,
            subslabels,
            ax=ax,
            cmap=custom_cmap("Greys", first_color="#ffffff00"),
            # cmap = "PiYG",
            cbar_kw={
                "drawedges": False,
                "shrink": 0.3,
                "pad": -0.05,
                "aspect": 10,
            },
            vmin=0,
            vmax=maxval,
        )
        # rotate cbar labels -90
        cbar2.set_label(f"{cbarlab} ≥{height}", rotation=90, va="bottom", labelpad=10)

    im, cbar = heatmap(
        fp_hm_array,
        # [(height+1-i) for i in range(1, height + 1)],
        yaxlabels,
        subslabels,
        ax=ax,
        cmap=custom_cmap(color_scheme, first_color="#ffffff00"),
        # cmap = "PiYG",
        cbarlabel=cbarlab,
        vmin=0,
        cbar_kw={"drawedges": False, "shrink": 0.3, "pad": 0.02, "aspect": 10},
    )
    cbar.set_label(f"{cbarlab}", rotation=90, va="bottom", labelpad=10)

    # texts = annotate_heatmap(im, valfmt="{x:.1f}")
    if annotate:
        texts = annotate_heatmap(im, valfmt="{x:.0f}", size=7)
    if standard_colour:
        set_label_colors_from_categories(
            ax.get_xticklabels(),
            get_pathway(version=defaultVersion),
            colourDict["pathways"],
        )
    # plt.figure(figsize=(10,6))
    ax.set_xlabel("substructure", labelpad=10)
    ax.set_ylabel("counts", labelpad=10)
    ax.set_title(title, loc="center", pad=20)
    fig.tight_layout()
    return fig


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
    hm = fp_heatmap(
        under,
        *args,
        percentages=True,
        top_acc_array=over,
        **kwargs,
    )
    return hm


def main():
    print("hello")
    set_style()
    ft = "svg"
    args = cli()
    fps = np.loadtxt(args.fingerprints, dtype=int, delimiter=",")

    fp_name = "_".join(args.fingerprints.split("/")[-1].split(".")[0].split("_")[1:])
    set_name = args.fingerprints.split("/")[-1].split(".")[0].split("_")[0]

    iwd = os.getcwd()
    os.makedirs(f"{iwd}/heatmaps/{set_name}/{fp_name}", exist_ok=True)
    os.chdir(f"{iwd}/heatmaps/{set_name}/{fp_name}")
    fp_name = fp_name.replace("_", " ")

    # save mean, std, median, etc of fps
    write_stats(fps)

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
            standard_colour=True,
        )
        plt.savefig(f"{fp_name}_heatmap.{ft}")
        return None

    hm = fp_heatmap_accumulative(
        fps,
        limit=10,
        title=f"Distribution of {fp_name} substructure counts",
        subslabels=substructure_names,
        color_scheme="GnBu",
        standard_colour=True,
    )
    savefig(hm, f"heatmap.{ft}")

    if args.classes:
        classes = np.loadtxt(args.classes, dtype="str", delimiter="\t", usecols=0)
        classes[classes == "fatty_acid,isoprenoid"] = "isoprenoid"
        classes = np.where(
            np.core.defchararray.find(classes, ",") != -1, "multiple", classes
        )
        classes[classes == ""] = "None"

        if len(classes) != fps.shape[0]:
            logging.warning(
                "classes file not same length as fingerprints file; will lead to errors"
            )
        for classif in np.unique(classes):
            idx = np.where(classes == classif)
            focus = fps[idx]
            if not classif:
                classif = "None"

            hm = fp_heatmap_accumulative(
                focus,
                limit=10,
                title=f"Distribution of {fp_name} substructure counts for {len(focus)} {classif.replace('_', ' ')} compounds",
                subslabels=substructure_names,
                color_scheme="GnBu",
                standard_colour=True,
            )
            savefig(hm, f"{classif}_heatmap.{ft}")

    os.chdir(iwd)
    # fp_means_plots(coco_mean, zinc_mean, outfile_namer(f"{coco_name}_{zinc_name}.svg"))


if __name__ == "__main__":
    main()
