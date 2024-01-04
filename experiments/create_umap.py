import argparse

import numpy as np
import umap
import matplotlib.pyplot as plt

from utils import set_style


def cli():
    parser = argparse.ArgumentParser()
    set_style()
    parser.add_argument("fingerprint", help="Fingerprint file")
    parser.add_argument("labels", help="Labels file")
    parser.add_argument(
        "-s",
        "--synthetic",
        required=False,
        type=str,
        help="fingerprints of synthetic compounds",
    )
    return parser.parse_args()


def main():
    args = cli()
    fp = np.loadtxt(args.fingerprint, delimiter=",", dtype=int)
    labels = np.loadtxt(args.labels, delimiter="\t", dtype=str, usecols=(0,))

    if args.synthetic:
        synthetic_fp = np.loadtxt(args.synthetic, delimiter=",", dtype=int)
        synthetic_fp = np.random.choice(
            synthetic_fp.shape[0], fp.shape[0], replace=False
        )
        fp = np.concatenate((fp, synthetic_fp))
        synthetic_labels = np.array(["synthetic" for _ in range(synthetic_fp.shape[0])])
        labels = np.concatenate((labels, synthetic_labels))

    # remove any where there is a , in the label
    idx = np.where(["," not in label for label in labels])
    fp = fp[idx]
    labels = labels[idx]

    # # random subsample
    # np.random.seed(42)
    # idx = np.random.choice(fp.shape[0], 10000, replace=False)
    # fp = fp[idx]
    # labels = labels[idx]

    # get unique labels
    unique_labels = list(np.unique(labels))
    print(type(unique_labels))
    # map labels to indexes
    label_to_idx = {str(label): int(idx) for idx, label in enumerate(unique_labels)}
    print(label_to_idx)
    # convert labels to indexes
    labels_i = np.array([label_to_idx[label] for label in list(labels)])

    print(fp.shape)
    print(labels.shape)
    reducer = umap.UMAP()
    embedding = reducer.fit_transform(fp)
    # plot umap with different colours for each label, and a legend on the right side
    s = plt.scatter(embedding[:, 0], embedding[:, 1], c=labels_i, cmap="Spectral")
    plt.legend(s.legend_elements()[0], list(set(labels)), loc="lower right")
    plt.xlabel("UMAP 1")
    plt.ylabel("UMAP 2")
    plt.title(f"UMAP of {args.fingerprint.split('/')[-1]}")
    plt.savefig("umap.png")

    np.savetxt("embedding.txt", embedding, fmt="%s")
    np.savetxt("labels.txt", labels, fmt="%s")


if __name__ == "__main__":
    main()
