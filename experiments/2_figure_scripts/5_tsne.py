# -*- coding: utf-8 -*-
"""
Created on Sun Oct  1 21:53:10 2023

@author: lucina-may
"""

import sys, os, logging, tracemalloc
from sys import argv
from datetime import datetime
from time import perf_counter

from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
import numpy as np
import pandas as pd

sys.path.append(os.path.abspath(os.path.join(sys.path[0], os.pardir, "src")))
from experiments.figure_producing.utils.figures import scatter, savefig
from utils import set_style


def __counted_tanimoto_sim(fp1: np.array, fp2: np.array) -> float:
    """
    Counted Tanimoto similarity taken from biosynfoni.fingerprints
    """
    nom = sum(np.minimum(fp1, fp2))  # overlap
    denom = float(sum(np.maximum(fp1, fp2)))  # all bits that are 'on'
    if denom == 0:
        return -1  # fingerprint does not describe either molecule
    else:
        return nom / denom


def _countanimoto(fps: list) -> float:
    """Counted Tanimoto similarity of list of two fingerprints"""
    return __counted_tanimoto_sim(np.array(fps[0]), np.array(fps[1]))


def outfile_namer(filename: str) -> str:
    "returns a filename with MMDD_ prefix"
    # get MMDD
    now = datetime.now()
    mmdd = now.strftime("%m%d")
    return f"{mmdd}_{filename}"


def fps_to_array(fingerprint_file: str) -> np.array:
    """returns numpy array of fingerprint file"""
    return np.loadtxt(fingerprint_file, delimiter=",", dtype=int)


def labels_to_array(label_file: str) -> np.array:
    """returns numpy array of label file"""
    return np.loadtxt(label_file, delimiter="\t", dtype=str)


def filter_fps(fps: np.array, labels: np.array) -> tuple[np.array, np.array]:
    """removes empty labels and fingerprints"""
    labels = np.where(labels == "", "None", labels)
    labels = np.where(labels == "fatty_acid,isoprenoid", "isoprenoid", labels)
    idx = np.where(["," not in label for label in labels])
    fps = fps[idx]
    labels = labels[idx]
    return fps, labels


def pcaer(data: np.array, n_components: int = 50, random_state=1):
    """
    PCA of data

        Args:
            data: numpy array of data
            n_components: int, default=50
            random_state: int, default=1
        Returns:
            pcaed: numpy array of PCA transformed data
            var: numpy array of explained variance ratios

    Remarks:
        - loadings are written to file
        -  random_state: Pass an int for reproducible results across multiple function calls.
    """
    pca = PCA(n_components=n_components, random_state=random_state)
    pca.fit(data)
    pcaed = pca.transform(data)
    var = pca.explained_variance_ratio_
    loadings = pca.components_
    with open(outfile_namer("pca_loadings.tsv"), "w") as f:
        for i, row in enumerate(loadings):
            f.write(f"{i}\t{row}\n")
    return pcaed, var


def distance_matrix(arr: np.array, metric: str = "euclidean") -> np.array:
    """returns distance matrix of array"""
    if metric == "euclidean":
        # return np.array([[np.linalg.norm(i - j) for j in arr] for i in arr])
        return arr  # array will automatically be converted to distance matrix euclidean in tsne
    if metric == "tanimoto":
        return np.array([[1.0 - _countanimoto([i, j]) for j in arr] for i in arr])
    if metric == "manhattan":
        # not checked!!!!!!!!!!!!
        print("warning: not checked!! check formula for manhattan distance")
        return np.array([[np.linalg.norm(i - j, ord=1) for j in arr] for i in arr])


def tsner(
    data: np.array,
    n_components: int = 2,
    perplexity: int = 40,
    n_iter: int = 300,
    *args,
    **kwargs,
) -> pd.DataFrame:
    """
    tSNE of data

    Args:
        data: numpy array of data
        n_components: int, default=2
        perplexity: int, default=40
        n_iter: int, default=300
        *args: additional arguments to pass to TSNE
        **kwargs: additional keyword arguments to pass to TSNE
    Returns:
        df: pandas DataFrame of tSNE transformed data's first two components
    """
    # tsne = TSNE(n_components=n_components, *args, **kwargs)
    tsne = TSNE(
        n_components=n_components, perplexity=perplexity, n_iter=n_iter, *args, **kwargs
    )
    z = tsne.fit_transform(data)
    df = pd.DataFrame()
    df["comp-1"] = z[:, 0]
    df["comp-2"] = z[:, 1]
    df.to_csv("tsne.tsv", index=True, header=True, sep="\t")
    return df


def get_first_ncps_est(annotfile: str):
    """Get the first class from first column of the annotation file"""
    npcs = []
    # with open('../arch/0914_COCONUT_DB_rdk_npcs.tsv','r') as cl:
    with open(annotfile, "r") as cl:
        for line in cl:
            npcs.append(line.strip().split("\t")[0].split(",")[0])
    return npcs


def annotate_df(df: pd.DataFrame, col: str, annot: list) -> pd.DataFrame:
    """Annotate a dataframe with a list of annotations"""
    df[col] = annot
    df[col] = df[col].replace("", "None")
    return df


def pca_plot(arr, labels: list[str], fp_name="biosynfoni") -> None:
    """
    Plot PCA of array

        Args:
            arr: numpy array of data
            labels: list of labels
            fp_name: str, default="biosynfoni"
        Returns:
            None

    Remarks:
        - saves pca plot to file with util.figuremaking.savefig
    """
    # pca
    logging.info("starting pca...")
    n_components = len(arr[0])
    n_components = 2
    logging.info("running pca function...")
    pcaed, var = pcaer(arr, n_components=n_components, random_state=333)
    logging.info("making dataframe...")
    pca_df = pd.DataFrame(pcaed)
    pca_df.columns = [f"{_a} ({_b})" for _a, _b in zip(pca_df.columns, var)]
    logging.info("annotating dataframe...")
    # df = annotate_df(pca_df, "npclassifier", get_first_ncps_est(annotfile)[: len(arr)])
    df = annotate_df(pca_df, "class", labels)
    logging.info("plotting pca...")
    sc = scatter(
        df,
        df.columns[0],
        df.columns[1],
        figtitle=f"pca of {fp_name}",
        color_by="class",
    )
    savefig(sc, "pca")
    return None


def pcaed_tsne(
    arr: np.array,
    labels: list[str],
    metric: str = "euclidean",
    verbose: int = 0,
    fp_name: str = "biosynfoni",
    perplexity: int = 50,
    n_iter: int = 500,
    *args,
    **kwargs,
) -> None:
    """
    Plot tSNE of array (2 components, unchangeable)

        Args:
            arr: numpy array of data
            labels: list of labels
            metric: str, default="euclidean"
            initial_pca_components: int, default=10
            verbose: int, default=0
            fp_name: str, default="biosynfoni"
            perplexity: int, default=50
            n_iter: int, default=500
            *args: additional arguments to pass to TSNE
            **kwargs: additional keyword arguments to pass to TSNE
        Returns:
            None

    Remarks:
        - saves tSNE plot to file with util.figuremaking.savefig
        - will run pca if arr.shape[1] > 39 to have same number of components as biosynfoni
    """

    n_components = 39  # to mimic biosynfoni
    arr = distance_matrix(arr, metric="euclidean")
    time_add = 0
    if arr.shape[1] > n_components:
        logging.info(f"running pca {arr.shape[1]}>{n_components} for initialisation...")
        tic0 = perf_counter()
        arr, _ = pcaer(arr, n_components=n_components, random_state=333)
        toc0 = perf_counter()
        time_add = toc0 - tic0
        logging.info(f"pca took {time_add:0.4f} seconds")

    logging.info("running tsne on precomputed distance matrix...")

    tic = perf_counter()
    tsne_comp = tsner(
        arr,
        n_components=2,
        perplexity=perplexity,
        n_iter=n_iter,
        verbose=verbose,
        metric="euclidean" if metric == "euclidean" else "precomputed",
        init="random",  # pca init cannot be used with precomputed distance matrices
        *args,
        **kwargs,
    )
    toc = perf_counter()
    logging.info(f"tsne took {toc - tic + time_add:0.4f} seconds")

    df = annotate_df(tsne_comp, "class", labels)
    sc = scatter(
        df,
        df.columns[0],
        df.columns[1],
        figtitle=f"tSNE of {fp_name}",
        color_by="class",
    )
    savefig(sc, f"tsne_{metric}")
    return None


def save_tsne_settings(
    dic: dict, filename: str = "tsne_settings", extra: str = ""
) -> None:
    """Save settings to file"""
    with open(f"{filename}.txt", "w") as s:
        for key, val in dic.items():
            s.write(f"{key}:\t{val}\n")
        if extra:
            s.write(f"extra:\t{extra}\n")
    return None


def main():
    set_style()
    tracemalloc.start()
    logging.info("hello")
    logging.getLogger().setLevel(logging.INFO)
    fingerprintfile = argv[1]
    annotfile = argv[2]
    fingerprintfile = os.path.abspath(fingerprintfile)
    annotfile = os.path.abspath(annotfile)

    setname = fingerprintfile.split("/")[-2]
    fp_name = "_".join(fingerprintfile.split("/")[-1].split(".")[0].split("_")[1:])

    iwd = os.getcwd()
    os.makedirs(f"{iwd}/tsnes/{setname}/{fp_name}", exist_ok=True)
    os.chdir(f"{iwd}/tsnes/{setname}/{fp_name}")

    # annotfile: the npcs.tsv or other classification infromation for colour
    arr = fps_to_array(fingerprintfile)
    labels = labels_to_array(annotfile)
    arr, labels = filter_fps(arr, labels)
    # configuration for tsne:--------------------------
    tsne_settings = {
        "perplexity": 50,
        "n_iter": 2000,
        # "initial_pca_components": 10,
        "metric": "euclidean",
    }
    # -------------------------------------------------
    logging.info(tsne_settings)

    save_tsne_settings(
        tsne_settings,
        filename="tsne_settings",
    )

    pca_plot(arr, labels=labels, fp_name=fp_name)
    pcaed_tsne(
        arr,
        metric=tsne_settings["metric"],
        labels=labels,
        # initial_pca_components=tsne_settings["initial_pca_components"],
        verbose=1,
        perplexity=tsne_settings["perplexity"],
        n_iter=tsne_settings["n_iter"],
        fp_name=fp_name,
    )
    tracemalloc.stop()
    logging.info("done")
    _, peak = tracemalloc.get_traced_memory()
    np.savetxt("peak_memory.txt", [peak / 10**6], fmt="%.4f")
    os.chdir(iwd)
    exit(0)
    return None


if __name__ == "__main__":
    main()
