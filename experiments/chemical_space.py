# -*- coding: utf-8 -*-
"""
Created on Sun Oct  1 21:53:10 2023

@author: lucina-may
"""

import sys, os
from sys import argv

from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
import numpy as np
import pandas as pd

sys.path.append(os.path.abspath(os.path.join(sys.path[0], os.pardir, "src")))
from utils import figuremaking as fm
from utils.figuremaking import df_scatterplot

# for intra-biosynfoni-code running
sys.path.append(
    os.path.abspath(os.path.join(sys.path[0], os.pardir, "src", "biosynfoni"))
)
from biosynfoni.inoutput import outfile_namer


def biosyfonis_to_array(fingerprint_file: str) -> np.array:
    data = []
    with open(fingerprint_file, "r") as f:
        for line in f:
            data.append([int(x) for x in line.strip().split(",")])
    data = np.array(data)
    return data


def pcaer(data: np.array, n_components: int = 50, random_state=1):
    "random_state: Pass an int for reproducible results across multiple function calls."
    pca = PCA(n_components=n_components, random_state=random_state)
    pca.fit(data)
    pcaed = pca.transform(data)
    var = pca.explained_variance_
    return pcaed, var


def tsner(
    data: np.array,
    n_components: int = 2,
    perplexity: int = 40,
    n_iter: int = 300,
    *args,
    **kwargs,
) -> pd.DataFrame:
    # tsne = TSNE(n_components=n_components, *args, **kwargs)
    tsne = TSNE(
        n_components=n_components, perplexity=perplexity, n_iter=n_iter, *args, **kwargs
    )
    z = tsne.fit_transform(data)
    df = pd.DataFrame()
    df["comp-1"] = z[:, 0]
    df["comp-2"] = z[:, 1]
    return df


def get_first_ncps_est(annotfile: str):
    npcs = []
    # with open('../arch/0914_COCONUT_DB_rdk_npcs.tsv','r') as cl:
    with open(annotfile, "r") as cl:
        for line in cl:
            npcs.append(line.strip().split("\t")[0].split(",")[0])
    return npcs


def annotate_df(df: pd.DataFrame, col: str, annot: list) -> pd.DataFrame:
    df[col] = annot
    df[col] = df[col].replace("", "None")
    return df


def get_subset(df: pd.DataFrame, n: int = 10000) -> pd.DataFrame:  # remove
    subset_df = df.sample(n=n)
    # print(subset_df.index().tolist)
    return subset_df


def _get_combi_annot(arr_np, arr_syn, np_annotfile: str = "npcs.tsv"):
    annot = []
    npcs = get_first_ncps_est(np_annotfile)#[: len(arr_np)]
    syns = ["synthetic" for x in range(len(arr_syn))]
    #annot.append(x for x in npcs)
    #annot.append("synthetic" for x in range(len(arr_syn)))
    return npcs+syns


def pca_plot(
    arr, annotation: list[str], filename="syncom_pca_bsf", randomseed=333
) -> None:
    # pca
    print("starting pca...")
    n_components = len(arr[0])
    n_components = 2
    print("running pca function...")
    pcaed, var = pcaer(arr, n_components=n_components, random_state=randomseed)
    print("making dataframe...")
    pca_df = pd.DataFrame(pcaed)
    pca_df.columns = [f"{_a} ({_b})" for _a, _b in zip(pca_df.columns, var)]
    print("annotating dataframe...")
    df = annotate_df(
        pca_df, "class", annotation
    )  # get_first_ncps_est(annotfile)[: len(arr)])
    print("plotting pca...")
    df_scatterplot(
        df,
        df.columns[0],
        df.columns[1],
        figtitle=f"pca of biosynfoni, coloured on npclassifier prediction {var[0]}{var[1]}",
        marginal_x="box",
        marginal_y="box",
        color="class",
        filename=filename,
        auto_open=False,
        color_discrete_map=fm.COLOUR_DICT["class"],
    )
    return None


def pcaed_tsne(
    arr: np.array,
    annot: list[str] = [],
    initial_pca_components: int = 10,
    verbose: int = 1,
    filename: str = "tsne_pcaed_bsf",
    perplexity: int = 50,
    n_iter: int = 500,
    randomseed=333,
    *args,
    **kwargs,
) -> None:
    """
    produces tsne plot of fingerprints, coloured on classification
    initial_pca_components: number of components to reduce to before running tsne, default 10, if higher than len(arr[0]), defaults to len(arr[0])
    """
    # run initial pca to reduce computational time
    if initial_pca_components > len(arr[0]):
        initial_pca_components = len(arr[0])
    pcaed, _ = pcaer(arr, n_components=initial_pca_components, random_state=randomseed)

    print("running tsne...")
    tsne_comp = tsner(
        pcaed,
        n_components=2,
        perplexity=50,
        n_iter=500,
        verbose=verbose,
        *args,
        **kwargs,
    )
    if annot:
        df = annotate_df(tsne_comp, "class", annot)
    df_scatterplot(
        df,
        df.columns[0],
        df.columns[1],
        figtitle="tsne of pca'ed biosynfoni of coconut & zinc (synthetic only), coloured on npclassifier prediction / synthetic status",
        marginal_x="box",
        marginal_y="box",
        color="class",
        filename=filename,
        auto_open=False,
        hover_data={"index": ("|%B %d, %Y", df.index)},
        color_discrete_map=fm.COLOUR_DICT["class"],
    )
    return None


def save_tsne_settings(
    dic: dict, filename: str = "tsne_settings", extra: str = ""
) -> None:
    with open(f"{filename}.txt", "w") as s:
        for key, val in dic.items():
            s.write(f"{key}:\t{val}\n")
        if extra:
            s.write(f"extra:\t{extra}\n")
    return None


def main():
    print("hello")
    fingerprintfile_coco = argv[1]  # natural products
    fingerprintfile_zinc = argv[2]  # synthetic compounds
    coco = biosyfonis_to_array(fingerprintfile_coco)
    zinc_toolarge = biosyfonis_to_array(fingerprintfile_zinc)
    # select random subset of zinc, with seed for reproducibility
    np.random.seed(333)
    zinc = zinc_toolarge[np.random.choice(zinc_toolarge.shape[0], len(coco), replace=False)]
    #zinc = np.random.choice(zinc_toolarge, size=len(coco), replace=False)
    assert len(coco[0]) == len(zinc[0]), "fingerprint lengths not equal, check input"

    # now, we concatenate the list
    arr = np.concatenate((coco, zinc))

    # annotfile: the npcs.tsv or other classification infromation for colour
    annotfile = argv[3]  # for natural products classification
    annotation = _get_combi_annot(coco, zinc, annotfile)
    assert len(annotation) == len(arr), "annotation number not equal to molecule number"

    # configuration for tsne:--------------------------
    tsne_settings = {
        "perplexity": 50,
        "n_iter": 2000,
        "initial_pca_components": 10,
        "max_initial_pca_components": len(arr[0]),
        "random_seed": 333,
    }

    # -------------------------------------------------

    save_tsne_settings(
        tsne_settings,
        filename=outfile_namer("tsne_settings"),
        extra=f"{fingerprintfile_coco}_{fingerprintfile_zinc}",
    )

    pca_plot(
        arr,
        annotation,
        filename=outfile_namer("syncom_pca_bsf"),
        randomseed=tsne_settings["random_seed"],
    )
    pcaed_tsne(
        arr,
        annotation,
        initial_pca_components=tsne_settings["initial_pca_components"],
        verbose=1,
        perplexity=tsne_settings["perplexity"],
        n_iter=tsne_settings["n_iter"],
        filename=outfile_namer("syncom_tsne_pcaed_bsf"),
        randomseed=tsne_settings["random_seed"],
    )

    print("done")
    sys.exit()
    return None


if __name__ == "__main__":
    main()
