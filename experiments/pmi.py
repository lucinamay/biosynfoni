#!/usr/bin/env python3
import argparse, logging
import math
from collections import Counter

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import svds
from tqdm import tqdm

from utils import set_style
from utils.figuremaking import set_label_colors, annotate_heatmap


def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", required=True, type=str)
    parser.add_argument("-o", required=True, type=str)
    return parser.parse_args()


def parse_fingerprints(path: str) -> np.array:
    fps = []

    with open(path, "r") as f:
        for line in tqdm(f):
            fp = line.strip().split(",")
            fp = [int(x) for x in fp]
            fps.append(np.array(fp))

    return np.array(fps)


def add_minuses(heatmap, array):
    for i in range(array.shape[0]):
        for j in range(array.shape[1]):
            text = ""
            # if array[i, j] > 0:
            #     text = "+"
            if array[i, j] < 0:
                text = "-"
            heatmap.text(
                j + 0.5,
                i + 0.5,
                text,
                horizontalalignment="center",
                verticalalignment="center",
                color="white",
            )
    return None


def get_labels() -> list[str]:
    keys = [
        "coa",
        "nadh",
        "nadph",
        "all standard aminos",
        "non-standard aminos",
        "open pyranose",
        "open furanose",
        "pyranose",
        "furanose",
        "indoleC2N",
        "phenylC2N",
        "c5n",
        "c4n",
        "phenylC3",
        "phenylC2",
        "phenylC1",
        "isoprene",
        "acetyl",
        "methylmalonyl",
        "ethyl",
        "methyl",
        "phosphate",
        "sulfonate",
        "fluorine",
        "chlorine",
        "bromine",
        "iodine",
        "nitrate",
        "epoxy",
        "ether",
        "hydroxyl",
        "c3 ring",
        "c4 ring",
        "c5 ring",
        "c6 ring",
        "c7 ring",
        "c8 ring",
        "c9 ring",
        "c10 ring",
    ]
    return keys


def get_colours() -> list[str]:
    colours = [
        "grey",
        "grey",
        "grey",
        "#FFEAA0",
        "#FFEAA0",
        "#FFC4CE",  # pink
        "#FFC4CE",  # pink
        "#FFC4CE",  # pink
        "#FFC4CE",  # pink
        "#A783B6",
        "#A783B6",
        "#FF8B61",
        "#FF8B61",
        "#A783B6",
        "#A783B6",
        "#A783B6",
        "#B9C311",  # green
        "#FF8B61",
        "#FF8B61",
        "#FF8B61",
        "#FF8B61",
        "grey",
        "grey",
        "grey",
        "grey",
        "grey",
        "grey",
        "grey",
        "grey",
        "grey",
        "grey",
        "grey",
        "grey",
        "grey",
        "grey",
        "grey",
        "grey",
        "grey",
        "grey",
    ]
    return colours


def _pearson_correlation_matrix(fps: np.array) -> np.array:
    # randomly subsample 10000

    # fps = fps[np.random.choice(fps.shape[0], 10000, replace=False), 3:]
    # fps = fps[:, 3:]
    mat = np.corrcoef(fps, rowvar=False, dtype=np.float16)
    # print(fps.shape)
    # print(mat.shape)
    return mat


def correlation_heatmap(fps: np.array) -> None:
    keys = get_labels()
    correlations = _pearson_correlation_matrix(fps)
    np.savetxt("correlations.tsv", correlations, delimiter="\t", fmt="%.2f")
    logging.warning(correlations.shape)
    hm = sns.heatmap(
        correlations,
        xticklabels=keys,
        yticklabels=keys,
        # cmap="coolwarm",
        cmap="PiYG",  # more colorblind-friendly diverging colormap
        vmin=np.min(fps),
        vmax=np.max(fps),
        center=0,  # center of the colormap
    )
    # hm.text(0.5, 0.5, "test", horizontalalignment="center", verticalalignment="center")
    # for all negative values, add a minus sign
    add_minuses(hm, fps)
    colours = get_colours()
    xt, yt = hm.get_xticklabels(), hm.get_yticklabels()
    set_label_colors(hm.get_xticklabels(), colours)
    set_label_colors(hm.get_yticklabels(), colours)
    return hm


def main() -> None:
    args = cli()
    fps = parse_fingerprints(args.i)
    title_text = args.o.replace(".png", "").replace("_", " ")

    set_style()

    # subsample 1000
    # fps = fps[:1000]

    # Count the number of times each bit is set.
    cx = Counter()
    cxy = Counter()

    for idx in tqdm(range(fps.shape[0])):
        for bit_idx, bit in enumerate(fps[idx]):
            if bit > 0:
                cx[bit_idx] += bit

            for bit_idx2, bit2 in enumerate(fps[idx]):
                # if bit_idx == bit_idx2:
                #     continue
                if bit > 0 and bit2 > 0:
                    # cxy[(bit_idx, bit_idx2)] += 1
                    cxy[(bit_idx, bit_idx2)] += min(bit, bit2)

    # Create lookup between key and fingerprint index.
    x2i, i2x = {}, {}
    keys = get_labels()
    for i, x in enumerate(keys):
        x2i[x] = i
        i2x[i] = x

    # Build sparse PMI matrix.
    sx = sum(cx.values())
    sxy = sum(cxy.values())
    data, rows, cols = [], [], []
    for (x, y), n in cxy.items():
        rows.append(x)
        cols.append(y)
        data.append(math.log((n / sxy) / (cx[x] / sx) / (cx[y] / sx)))

    PMI = csc_matrix((data, (rows, cols)))
    mat = PMI.toarray()

    # Visualize matrix.
    fig = plt.figure(figsize=(8, 6))
    # sns.set(font_scale=0.5)
    hm = sns.heatmap(
        mat,
        xticklabels=keys,
        yticklabels=keys,
        # cmap="coolwarm",
        cmap="PiYG",  # more colorblind-friendly diverging colormap
        vmin=np.min(mat),
        vmax=np.max(mat),
        center=0,  # center of the colormap
    )
    # hm.text(0.5, 0.5, "test", horizontalalignment="center", verticalalignment="center")
    # for all negative values, add a minus sign
    add_minuses(hm, mat)

    # set_label_colors(hm.get_xticklabels(), colours)
    # set_label_colors(hm.get_yticklabels(), colours)
    colours = get_colours()
    set_label_colors(hm.get_xticklabels(), colours)
    set_label_colors(hm.get_yticklabels(), colours)

    # plt.xticks(rotation=90)
    # plt.yticks(rotation=0)

    plt.ylabel("Bit 1")
    plt.xlabel("Bit 2")
    plt.title(
        f"Pointwise Mutual Information - non-overlap Biosynfoni on COCONUT", size=16
    )
    # plt.title(f"Pointwise Mutual Information - {title_text}", size=22)

    plt.savefig(args.o, bbox_inches="tight")
    plt.close()

    # get a correlation heatmap as well
    corr_hm = correlation_heatmap(fps)
    plt.title("Pearson correlation - non-overlap Biosynfoni on COCONUT", size=16)
    plt.savefig(args.o.replace(".png", "_correlation.png"), bbox_inches="tight")
    plt.close()

    exit(0)


if __name__ == "__main__":
    main()
