import argparse, logging, os, tracemalloc
from time import perf_counter

import numpy as np
import pandas as pd
import umap
import matplotlib as mpl
import matplotlib.pyplot as plt

from utils import set_style
from utils.colours import colourDict as col_dicts
from utils.figuremaking import scatter, savefig


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
    parser.add_argument(
        "-o", "--output", required=False, default="umap.svg", help="Output file"
    )
    parser.add_argument(
        "-a", "--additional", required=False, help="Additional labels file"
    )

    args = parser.parse_args()
    if not "." in args.output:
        args.output = f"{args.output}.svg"

    args.fingerprint = os.path.abspath(args.fingerprint)
    args.labels = os.path.abspath(args.labels)
    if args.smiles:
        args.smiles = os.path.abspath(args.smiles)
    if args.synthetic:
        args.synthetic = os.path.abspath(args.synthetic)
    if args.additional:
        args.additional = os.path.abspath(args.additional)
    return args


def clean_labels(fp, labels):
    """
    removes any where there is a ',' in the label
    """
    idx = np.where(["," not in label for label in labels])
    fp = fp[idx]
    labels = labels[idx]
    return fp, labels


# def onpick(event):
#     fp = mols_info["fp"]
#     labels = mols_info["labels"]
#     smiles = mols_info["smiles"]
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
#         annot = ax.text(
#             (embedding[i, 0] + 0.1),  # x coordinate + 2 to the right to avoid overlap
#             (embedding[i, 1] + 0.05),  # y coordinate + 2 to the right to avoid overlap
#             f"{i} {smiles[i]}",  # text
#             size=2,
#             zorder=1,
#             color="k",
#         )
#         annotations.append(annot)
#     # force redraw
#     fig.canvas.draw_idle()
#     return annotations


# def get_cmap_norm(labels, colour_dict):
#     c = np.random.randint(1, 5, size=N)
#     colors = [f"C{i}" for i in np.arange(1, c.max() + 1)]
#     cmap, norm = mpl.colors.from_levels_and_colors(np.arange(1, c.max() + 2), colors)


# def umap_2d(embedding, mols_info):
#     annotations = []

#     colors = mols_info["colors"]
#     labels_cl = mols_info["labels_cl"]

#     # plot umap with different colours for each label, and a legend on the right side
#     fig = plt.figure(figsize=(3, 3))
#     ax = fig.add_subplot(111)
#     s = plt.scatter(
#         embedding[:, 0],
#         embedding[:, 1],
#         # c=labels_i,
#         # cmap="Spectral",
#         c=colors,
#         # alpha=0.5,
#         edgecolors="none",
#         # picker=True,
#     )
#     s.set_alpha(0.5)  # set afterwards

#     # set background color
#     ax.set_facecolor("white")
#     # fig.canvas.mpl_connect("pick_event", onpick)

#     # # press 'e' to erase all annotations
#     # fig.canvas.mpl_connect(
#     #     "key_press_event",
#     #     lambda event: [
#     #         (annotation.remove() for annotation in annotations),
#     #         fig.canvas.draw_idle() if event.key == "e" else None,
#     #     ],
#     # )

#     # legend
#     print(ax.get_legend_handles_labels())
#     ax.legend(
#         handles=s.legend_elements()[0],
#         labels=labels_cl,
#         loc="lower right",
#         title="Classes",
#         # draw far outside the plot
#         bbox_to_anchor=(1.1, 0.0),
#         borderaxespad=0,
#     )

#     # ax.set_zlabel("UMAP 3")
#     return fig, ax


# def umap_3d(embedding, mols_info):
#     fp = mols_info["fp"]
#     labels = mols_info["labels"]
#     smiles = mols_info["smiles"]
#     colors = mols_info["colors"]
#     labels_cl = mols_info["labels_cl"]
#     #  def onpick(event):
#     #     print("onpick scatter")
#     #     ind = event.ind
#     #     print(
#     #         "onpick scatter:",
#     #         ind,
#     #         fp[ind],
#     #         labels[ind],
#     #         smiles[ind],
#     #         # np.take(embedding[:, 0], ind),
#     #         # np.take(embedding[:, 1], ind),
#     #         # np.take(embedding[:, 2], ind),
#     #     )
#     #     for i in ind:
#     #         ax.text(
#     #             embedding[i, 0] + 2,  # x coordinate + 2 to the right to avoid overlap
#     #             embedding[i, 1] + 2,  # y coordinate + 2 to the right to avoid overlap
#     #             embedding[i, 2] + 2,  # z coordinate + 2 to the right to avoid overlap
#     #             f"{i} {smiles[i]}",  # text
#     #             size=5,
#     #             zorder=1,
#     #             color="k",
#     #         )
#     #     fig.canvas.draw_idle()  # redraw the canvas
#     #     return None

#     # plot umap in three dimensions with adjustable viewing angle ================================
#     fig = plt.figure(figsize=(3, 3))
#     ax = fig.add_subplot(111, projection="3d")
#     s = ax.scatter(
#         embedding[:, 0],
#         embedding[:, 1],
#         embedding[:, 2],
#         # c=labels_i,
#         # cmap="Spectral",
#         c=colors,
#         alpha=0.5,
#         s=3,
#         edgecolors="none",
#         picker=True,
#     )
#     ax.view_init(
#         30, 45
#     )  # ==================== change viewing angle here ====================

#     fig.canvas.mpl_connect("pick_event", onpick3)

#     # add hover information within the scatterplot, at (10,10,10) from the data point that is hovered over

#     # also show the label in the left corner
#     ax.text(
#         0.05,
#         0.95,
#         f"{labels[ind[0]]}",
#         transform=ax.transAxes,
#         size=10,
#         zorder=1,
#         color="k",
#     )
#     fig.canvas.draw_idle()  # redraw the canvas

#     fig.canvas.mpl_connect("pick_event", onpick)
#     # press 'e' to erase all annotations
#     fig.canvas.mpl_connect(
#         "key_press_event",
#         lambda event: [
#             (annotation.remove() for annotation in annotations),
#             fig.canvas.draw_idle() if event.key == "e" else None,
#         ],
#     )

#     # for i, txt in enumerate(labels):
#     #     ax.annotate(txt, (embedding[i, 0], embedding[i, 1], embedding[i, 2]))
#     #

#     # # # Add annotations and legends

#     # # plt.xlabel("UMAP 1")
#     # # plt.ylabel("UMAP 2")
#     # # ax.set_zlabel("UMAP 3")
#     # # plt.title(f"UMAP of {args.fingerprint.split('/')[-1]}")
#     # # plt.savefig("umap.png")
#     # # add annotations and legends, making the legend max 1/10th of the plot
#     # ax.legend(
#     #     loc="lower right",
#     #     title="Classes",
#     #     # draw far outside the plot
#     #     bbox_to_anchor=(1.1, 0.0),
#     #     borderaxespad=0,
#     # )
#     # ax.set_xlabel("UMAP 1")
#     # ax.set_ylabel("UMAP 2")
#     # ax.set_zlabel("UMAP 3")
#     return None


def main():
    set_style()

    args = cli()

    colour_dict = col_dicts["class"]

    setname = args.fingerprint.split("/")[-2]
    fp_name = "_".join(args.fingerprint.split("/")[-1].split(".")[0].split("_")[1:])

    iwd = os.getcwd()
    os.makedirs(f"{iwd}/umap/{setname}/{fp_name}", exist_ok=True)
    os.chdir(f"{iwd}/umap/{setname}/{fp_name}")

    fp = np.loadtxt(args.fingerprint, delimiter=",", dtype=int)
    # fp[fp > 1] = 1  # binarize
    labels = np.loadtxt(args.labels, delimiter="\t", dtype="str", usecols=0)
    # replace '' with 'none' and 'fatty_acid,isoprenoid' with 'isoprenoid'
    labels = np.where(labels == "", "None", labels)
    labels = np.where(labels == "fatty_acid,isoprenoid", "isoprenoid", labels)
    smiles = np.zeros(labels.shape[0]).astype(str)
    labels2 = np.zeros(labels.shape[0]).astype(str)

    if args.smiles:
        smiles = np.loadtxt(args.smiles, delimiter="\t", dtype="str", usecols=(0,))
    if args.additional:
        #  reaplace the values of labels2 with the values from the additional file
        labels2 = np.loadtxt(
            args.additional,
            delimiter="\t",
            dtype="str",
            usecols=1,
        )
        labels2 = np.where(labels2 == "", "No NP-Classifier prediction", labels2)
        labels2 = np.where(labels2 == "no_prediction", "None", labels2)
        # change any labels with a "," in them to "None"
        labels2 = np.where(
            np.core.defchararray.find(labels2, ",") != -1, "Multiple", labels2
        )

        print(labels2.shape, labels.shape, fp.shape, smiles.shape)

    if args.synthetic:
        synthetic_fp = np.loadtxt(args.synthetic, delimiter=",", dtype=int)
        ind = np.random.choice(synthetic_fp.shape[0], fp.shape[0], replace=False)
        synthetic_fp = synthetic_fp[ind]
        synthetic_labels = np.array(["synthetic" for _ in range(synthetic_fp.shape[0])])
        synthetic_smiles = np.array(["synthetic" for _ in range(synthetic_fp.shape[0])])
        fp = np.concatenate((fp, synthetic_fp))
        labels = np.concatenate((labels, synthetic_labels))
        smiles = np.concatenate((smiles, synthetic_smiles))
        labels2 = np.concatenate((labels2, synthetic_labels))

    # remove any where there is a , in the label
    idx = np.where(["," not in label for label in labels])
    fp = fp[idx]
    labels = labels[idx]
    smiles = smiles[idx]
    labels2 = labels2[idx]

    # remove any where there is a * in the smiles
    if args.smiles is not None:
        idx = np.where(["*" not in smile for smile in smiles])
        fp = fp[idx]
        labels = labels[idx]
        smiles = smiles[idx]
        labels2 = labels2[idx]

        # remove any where the smiles is the same as others
        idx = np.where([smiles[i] not in smiles[:i] for i in range(smiles.shape[0])])
        fp = fp[idx]
        labels = labels[idx]
        smiles = smiles[idx]
        labels2 = labels2[idx]

    logging.info(fp.shape, labels.shape, smiles.shape)
    # # random subsample
    # np.random.seed(42)
    # idx = np.random.choice(fp.shape[0], 100000, replace=False)
    # fp = fp[idx]
    # labels = labels[idx]
    # smiles = smiles[idx]

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
        tracemalloc.start()
        tic = perf_counter()
        # fit an embedding to data except class 6
        embedding = reducer.fit_transform(fp[labels_i != s_label])
        toc = perf_counter()
        tracemalloc.stop()
        # show class 6 on that embedding
        embedding = reducer.transform(fp)
    else:
        tracemalloc.start()
        tic = perf_counter()
        embedding = reducer.fit_transform(fp)
        toc = perf_counter()
        tracemalloc.stop()

    curr, peak = tracemalloc.get_traced_memory()

    logging.info(f"umap took {toc - tic:0.4f} seconds")
    logging.info(f"Current memory usage is {curr / 10**6}MB; Peak was {peak / 10**6}MB")
    np.savetxt(
        "time_peakmem.tsv", [(toc - tic), (peak / 10**6)], delimiter="\t", fmt="%.4f"
    )
    # save loadings of the embedding
    np.savetxt("loadings.txt", reducer.embedding_, fmt="%s")

    # make df out of embedding
    df = pd.DataFrame(
        embedding, columns=[f"UMAP {i}" for i in range(embedding.shape[1])]
    )
    df["class"] = labels
    df["smiles"] = smiles

    sb = scatter(
        df,
        col_x="UMAP 0",
        col_y="UMAP 1",
        color_by="class",
        figtitle=f"UMAP of {fp_name.replace('_', ' ')}",
        # s=3,
        # figsize=(6, 10),
    )
    savefig(sb, args.output)
    if args.additional:
        df["class"] = labels2
        sb = scatter(
            df,
            col_x="UMAP 0",
            col_y="UMAP 1",
            color_by="class",
            figtitle=f"UMAP of {fp_name.replace('_', ' ')} labeled by npclassifier predictions",
            # s=3,
            # figsize=(6, 10),
        )
        savefig(sb, f"umap_additional.svg")
        df["class"] = labels
        df["class_2"] = labels2
    df.to_csv("embedding.tsv", sep="\t", index=True)

    exit(0)
    # # # # show loadings of the embedding
    # # # print(reducer.embedding_)

    # # annotations = []

    # # def onpick(event):
    # #     print("onpick scatter")
    # #     ind = event.ind
    # #     print(
    # #         "onpick scatter:",
    # #         ind,
    # #         fp[ind],
    # #         labels[ind],
    # #         smiles[ind],
    # #         # np.take(embedding[:, 0], ind),
    # #         # np.take(embedding[:, 1], ind),
    # #         # np.take(embedding[:, 2], ind),
    # #     )
    # #     # annotations = []  # make list for removing annotations
    # #     for i in ind:
    # #         annotation = ax.text(
    # #             (
    # #                 embedding[i, 0] + 0.1
    # #             ),  # x coordinate + 2 to the right to avoid overlap
    # #             (
    # #                 embedding[i, 1] + 0.05
    # #             ),  # y coordinate + 2 to the right to avoid overlap
    # #             f"{i} {smiles[i]}",  # text
    # #             size=2,
    # #             zorder=1,
    # #             color="k",
    # #         )
    # #         annotations.append(annotation)
    # #     # force redraw
    # #     fig.canvas.draw_idle()
    # #     return annotations

    # # # plot umap with different colours for each label, and a legend on the right side
    # # fig = plt.figure(figsize=(3, 3))
    # # ax = fig.add_subplot(111)
    # # s = plt.scatter(
    # #     embedding[:, 0],
    # #     embedding[:, 1],
    # #     # c=labels_i,
    # #     # cmap="Spectral",
    # #     c=colors,
    # #     # alpha=0.5,
    # #     edgecolors="none",
    # #     picker=True,
    # # )
    # # s.set_alpha(0.5)  # set afterwards

    # # # set background color
    # # ax.set_facecolor("white")

    # mols_info = {
    #     "fp": fp,
    #     "labels": labels,
    #     "smiles": smiles,
    #     "colors": colors,
    #     "labels_cl": np.array(list(label_to_idx.keys())),
    # }
    # umap_2d(embedding, mols_info)

    # plt.title(f"UMAP of {args.fingerprint.split('/')[-1]}")
    # # plt.show()
    # if args.output:
    #     plt.savefig(args.output)
    # else:
    #     plt.savefig(f"umap_{args.fingerprint.split('/')[-1]}.png")

    # # save embedding and labels
    # np.savetxt("embedding.txt", embedding, fmt="%s")
    # np.savetxt("labels.txt", labels, fmt="%s")


if __name__ == "__main__":
    main()
