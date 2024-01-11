import argparse, logging

import numpy as np
import umap
import matplotlib as mpl
import matplotlib.pyplot as plt

from utils import set_style
from utils.colours import COLOUR_DICT as col_dicts


def cli():
    parser = argparse.ArgumentParser()

    parser.add_argument("fingerprint", help="Fingerprint file")
    parser.add_argument("labels", help="Labels file")
    parser.add_argument("--smiles", required=False, help="Smiles file")
    parser.add_argument(
        "-s",
        "--synthetic",
        required=False,
        type=str,
        help="fingerprints of synthetic compounds",
    )
    return parser.parse_args()


def clean_labels(fp, labels):
    """
    removes any where there is a ',' in the label
    """
    idx = np.where(["," not in label for label in labels])
    fp = fp[idx]
    labels = labels[idx]
    return fp, labels


def umap_2d(embedding, mols_info):
    annotations = []

    fp = mols_info["fp"]
    labels = mols_info["labels"]
    smiles = mols_info["smiles"]
    colors = mols_info["colors"]
    labels_cl = mols_info["labels_cl"]

    def onpick(event):
        print("onpick scatter")
        ind = event.ind
        print(
            "onpick scatter:",
            ind,
            fp[ind],
            labels[ind],
            smiles[ind],
            # np.take(embedding[:, 0], ind),
            # np.take(embedding[:, 1], ind),
            # np.take(embedding[:, 2], ind),
        )
        for i in ind:
            annot = ax.text(
                (
                    embedding[i, 0] + 0.1
                ),  # x coordinate + 2 to the right to avoid overlap
                (
                    embedding[i, 1] + 0.05
                ),  # y coordinate + 2 to the right to avoid overlap
                f"{i} {smiles[i]}",  # text
                size=2,
                zorder=1,
                color="k",
            )
            annotations.append(annot)
        # force redraw
        fig.canvas.draw_idle()
        return annotations

    # plot umap with different colours for each label, and a legend on the right side
    fig = plt.figure(figsize=(3, 3))
    ax = fig.add_subplot(111)
    s = plt.scatter(
        embedding[:, 0],
        embedding[:, 1],
        # c=labels_i,
        # cmap="Spectral",
        c=colors,
        # alpha=0.5,
        edgecolors="none",
        picker=True,
    )
    s.set_alpha(0.5)  # set afterwards

    # set background color
    ax.set_facecolor("white")
    fig.canvas.mpl_connect("pick_event", onpick)
    # press 'e' to erase all annotations
    fig.canvas.mpl_connect(
        "key_press_event",
        lambda event: [
            (annotation.remove() for annotation in annotations),
            fig.canvas.draw_idle() if event.key == "e" else None,
        ],
    )
    ax.legend(
        # cannot get to
        handles=s.legend_elements()[0],
        # labels=label_to_idx.keys(),
        labels=labels_cl,
        loc="lower right",
        title="Classes",
        # draw far outside the plot
        bbox_to_anchor=(1.1, 0.0),
        borderaxespad=0,
    )
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    # ax.set_zlabel("UMAP 3")
    plt.title(f"UMAP of {args.fingerprint.split('/')[-1]}")
    plt.show()
    plt.savefig("umap.png")

    # save embedding and labels
    np.savetxt("embedding.txt", embedding, fmt="%s")
    np.savetxt("labels.txt", labels, fmt="%s")


def main():
    set_style()
    args = cli()

    colour_dict = col_dicts["class"]

    fp = np.loadtxt(args.fingerprint, delimiter=",", dtype=int)
    labels = np.loadtxt(args.labels, delimiter="\t", dtype=str, usecols=(0,))
    smiles = np.zeros(labels.shape[0]).astype(str)
    if args.smiles:
        smiles = np.loadtxt(args.smiles, delimiter="\t", dtype=str, usecols=(0,))

    if args.synthetic:
        synthetic_fp = np.loadtxt(args.synthetic, delimiter=",", dtype=int)
        ind = np.random.choice(synthetic_fp.shape[0], fp.shape[0], replace=False)
        synthetic_fp = synthetic_fp[ind]
        synthetic_labels = np.array(["synthetic" for _ in range(synthetic_fp.shape[0])])
        synthetic_smiles = np.array(["synthetic" for _ in range(synthetic_fp.shape[0])])
        fp = np.concatenate((fp, synthetic_fp))
        labels = np.concatenate((labels, synthetic_labels))
        smiles = np.concatenate((smiles, synthetic_smiles))

    # remove any where there is a , in the label
    idx = np.where(["," not in label for label in labels])
    fp = fp[idx]
    labels = labels[idx]
    smiles = smiles[idx]

    # remove any where there is a * in the smiles
    idx = np.where(["*" not in smile for smile in smiles])
    fp = fp[idx]
    labels = labels[idx]
    smiles = smiles[idx]

    # remove any where the smiles is the same as others
    idx = np.where([smiles[i] not in smiles[:i] for i in range(smiles.shape[0])])
    fp = fp[idx]
    labels = labels[idx]
    smiles = smiles[idx]

    # # random subsample
    # np.random.seed(42)
    # idx = np.random.choice(fp.shape[0], 10000, replace=False)
    # fp = fp[idx]
    # labels = labels[idx]

    # get unique labels
    unique_labels = list(np.unique(labels))

    # map labels to indexes
    label_to_idx = {str(label): int(idx) for idx, label in enumerate(unique_labels)}
    idx_to_label = {val: key for key, val in label_to_idx.items()}
    # print(label_to_idx)
    # print(idx_to_label)

    # convert labels to indexes
    labels_i = np.array([label_to_idx[label] for label in list(labels)])
    colors = [colour_dict[label] for label in labels]
    # make colors into matplotlib colors
    colors = [mpl.colors.to_rgba(color) for color in colors]

    logging.debug(fp.shape, labels_i.shape, smiles.shape)
    assert fp.shape[0] == labels_i.shape[0], "fp and labels must have same length"

    # fit an embedding to data
    reducer = umap.UMAP(n_components=2)

    if args.synthetic:
        s_label = label_to_idx["synthetic"]
        # fit an embedding to data except class 6
        embedding = reducer.fit_transform(fp[labels_i != s_label])
        # show class 6 on that embedding
        embedding = reducer.transform(fp)
    else:
        embedding = reducer.fit_transform(fp)

    # # show loadings of the embedding
    # print(reducer.embedding_)

    mols_info = [
        {
            "fp": fp[i],
            "labels": labels[i],
            "smiles": smiles[i],
            "colors": colors[i],
            "labels_cl": idx_to_label[labels_i[i]],
        }
        for i in range(fp.shape[0])
    ]
    # # plot umap with different colours for each label, and a legend on the right side
    # 2d_umap = umap_2d(embedding, mols_info)

    s.set_alpha(0.5)  # set afterwards

    # set background color
    ax.set_facecolor("white")

    # def onpick(event):
    #     print("onpick scatter")
    #     ind = event.ind
    #     print(
    #         "onpick scatter:",
    #         ind,
    #         fp[ind],
    #         labels[ind],
    #         smiles[ind],
    #         # np.take(embedding[:, 0], ind),
    #         # np.take(embedding[:, 1], ind),
    #         # np.take(embedding[:, 2], ind),
    #     )
    #     for i in ind:
    #         ax.text(
    #             embedding[i, 0] + 2,  # x coordinate + 2 to the right to avoid overlap
    #             embedding[i, 1] + 2,  # y coordinate + 2 to the right to avoid overlap
    #             embedding[i, 2] + 2,  # z coordinate + 2 to the right to avoid overlap
    #             f"{i} {smiles[i]}",  # text
    #             size=5,
    #             zorder=1,
    #             color="k",
    #         )
    #     fig.canvas.draw_idle()  # redraw the canvas

    # # plot umap in three dimensions with adjustable viewing angle ================================
    # fig = plt.figure(figsize=(3, 3))
    # ax = fig.add_subplot(111, projection="3d")
    # s = ax.scatter(
    #     embedding[:, 0],
    #     embedding[:, 1],
    #     embedding[:, 2],
    #     # c=labels_i,
    #     # cmap="Spectral",
    #     c=colors,
    #     alpha=0.5,
    #     s=3,
    #     edgecolors="none",
    #     picker=True,
    # )
    # ax.view_init(30, 45) #==================== change viewing angle here ====================

    # fig.canvas.mpl_connect("pick_event", onpick3)

    # add hover information within the scatterplot, at (10,10,10) from the data point that is hovered over

    # # also show the label in the left corner
    # ax.text(
    #     0.05,
    #     0.95,
    #     f"{labels[ind[0]]}",
    #     transform=ax.transAxes,
    #     size=10,
    #     zorder=1,
    #     color="k",
    # )
    # fig.canvas.draw_idle() # redraw the canvas

    fig.canvas.mpl_connect("pick_event", onpick)
    # press 'e' to erase all annotations
    fig.canvas.mpl_connect(
        "key_press_event",
        lambda event: [
            (annotation.remove() for annotation in annotations),
            fig.canvas.draw_idle() if event.key == "e" else None,
        ],
    )

    # for i, txt in enumerate(labels):
    #     ax.annotate(txt, (embedding[i, 0], embedding[i, 1], embedding[i, 2]))

    # # Add annotations and legends

    # plt.xlabel("UMAP 1")
    # plt.ylabel("UMAP 2")
    # ax.set_zlabel("UMAP 3")
    # plt.title(f"UMAP of {args.fingerprint.split('/')[-1]}")
    # plt.savefig("umap.png")
    # add annotations and legends, making the legend max 1/10th of the plot
    ax.legend(
        # cannot get to
        handles=s.legend_elements()[0],
        # labels=label_to_idx.keys(),
        labels=labels_cl,
        loc="lower right",
        title="Classes",
        # draw far outside the plot
        bbox_to_anchor=(1.1, 0.0),
        borderaxespad=0,
    )
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    # ax.set_zlabel("UMAP 3")
    plt.title(f"UMAP of {args.fingerprint.split('/')[-1]}")
    plt.show()
    plt.savefig("umap.png")

    # save embedding and labels
    np.savetxt("embedding.txt", embedding, fmt="%s")
    np.savetxt("labels.txt", labels, fmt="%s")


if __name__ == "__main__":
    main()
