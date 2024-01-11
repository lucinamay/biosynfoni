#!/usr/bin/env python3
import argparse
import math
from collections import Counter

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import svds
from tqdm import tqdm

from utils import set_style
from utils.figuremaking import label_colourcode


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


def main() -> None:
    args = cli()
    fps = parse_fingerprints(args.i)

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
    sns.set(font_scale=0.5)
    hm = sns.heatmap(
        mat,
        xticklabels=keys,
        yticklabels=keys,
        # cmap="coolwarm",
        cmap="PiYG",  # more colorblind-friendly diverging colormap
        vmin=np.min(mat),
        vmax=np.max(mat),
        center=0,
    )

    colours = [
        "grey",
        "grey",
        "grey",
        "yellow",
        "yellow",
        "pink",
        "pink",
        "pink",
        "pink",
        "purple",
        "purple",
        "orange",
        "orange",
        "purple",
        "purple",
        "purple",
        "teal",
        "orange",
        "orange",
        "orange",
        "orange",
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

    label_colourcode(hm.get_xticklabels(), colours)

    plt.ylabel("Bit 1")
    plt.xlabel("Bit 2")
    plt.title("Pointwise Mutual Information")

    plt.savefig(args.o, bbox_inches="tight")

    exit(0)


if __name__ == "__main__":
    main()
