#!/usr/bin/env python3
import logging, time, tracemalloc, sys, re
from typing import Generator
from pathlib import Path


import numpy as np
import networkx as nx
import pandas as pd
import tmap
from sklearn.manifold import TSNE
import tmap.api
import tmap.api.general
import tmap.tda
import tmap.tda.mapper
from umap import UMAP
import joblib
from tqdm import tqdm
from sklearn.metrics import (
    confusion_matrix,
    classification_report,
    roc_curve,
    multilabel_confusion_matrix,
)
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import KFold, GroupKFold
from sklearn.tree import export_graphviz
from scipy.spatial.distance import pdist, squareform
from tmap.tda.metric import Metric

from helper import ChangeDirectory
from biosynfoni.subkeys import get_names


def _read_fps(folder_path) -> Generator:
    for file in folder_path.glob("*.csv"):
        if "chebi_" in file.stem and not len(file.stem.split("_")) > 2:
            yield file.stem.replace("chebi_", ""), np.loadtxt(
                file, delimiter=",", dtype=int
            )


def _read_sim(folder_path) -> Generator:
    for file in folder_path.glob("*_sim.csv"):
        yield file.stem.replace("chebi_", "").replace("_sim", ""), np.loadtxt(
            file, delimiter=",", dtype=float
        )


# def unilabels_and_dict(classifications: np.array) -> tuple[np.array, dict]:
#     classifications = np.array(
#         [";".join(sorted(i.split(";"))) for i in classifications]
#     )
#     id_to_class = {i: class_ for i, class_ in enumerate(classifications.unique())}
#     id_to_class = {
#         i: "none" if class_ == "" else class_ for i, class_ in id_to_class.items()
#     }
#     return classifications.map(id_to_class), id_to_class


def multilabel_and_dict(classifications: np.array) -> tuple[np.array, dict]:
    classes = set(";".join(map(str, classifications)).split(";"))
    class_to_id = {class_: i for i, class_ in enumerate(sorted(classes))}
    id_to_class = {i: class_ for class_, i in class_to_id.items()}
    classification_array = np.zeros((len(classifications), len(classes)), dtype=int)
    for i, classification in enumerate(classifications):
        for class_ in re.split(r"[;,]", classification):
            classification_array[i, class_to_id[class_]] = 1
    return classification_array, id_to_class


def separate_nones(class_index: dict, x_y_andco: list):
    """
    Separate "None" classifications from the data. x_y_andco is a list of arrays to separate, with at least X and y.
    """
    # get the key for which the value is "None"
    ind_none = class_index["None"]

    # get instance indices of "None" in y

    x_y_andco = list(x_y_andco)
    X = x_y_andco[0]
    y = x_y_andco[1]
    if y.shape[1] > 1:
        none_entries = np.where(y[:, ind_none] == 1)[0]
        # remove the "None" column from y
        y = np.delete(y, ind_none, 1)
        for key in class_index.keys():
            if class_index[key] > ind_none:
                class_index[key] -= 1
    else:  # unilabel
        none_entries = np.where(y == ind_none)[0]

    nones = x_y_andco.copy()
    for i, array in enumerate(x_y_andco):
        x_y_andco[i] = array[~none_entries]
        nones[i] = array[none_entries]
    return class_index, x_y_andco, nones


def predict_one(
    X_test: np.array,
    classifier: RandomForestClassifier,
) -> np.array:
    return np.array(classifier.predict(X_test))


def kfold_performance(X, y, *args, **kwargs) -> tuple:
    """
    *args and **kwargs are passed to RandomForestClassifier
    """
    importances = []

    probabilities, ks = (np.empty(y.shape), np.empty(y.shape[0]))

    # random state has to be same to have same splits for each fingerprint
    for k, (train, test) in enumerate(KFold(shuffle=True, random_state=42).split(X)):
        X_train, X_test = X[train], X[test]
        y_train, y_test = y[train], y[test]
        classifier = RandomForestClassifier(*args, **kwargs).fit(X_train, y_train)
        importances.append(classifier.feature_importances_)

        # only probabilities of being in class
        y_proba = np.transpose(np.array(classifier.predict_proba(X_test))[:, :, 1])
        y_pred = (y_proba >= 0.5).astype(int)
        probabilities[test] = y_proba
        ks[test] = k

    return importances, probabilities, ks


def write_comb_confusion_matrix(
    cm: np.array,
    classification_types: list,
    title: str = "confusion_matrix.txt",
    unilabel: bool = False,
) -> None:
    """
    Write confusion matrix to file

        Args:
            cm (np.array): confusion matrix
            classification_types (list): list of classification types
            title (str): title of file. Default: "confusion_matrix.txt"
            unilabel (bool): if True, will write unilabel confusion matrix, if False, will write multilabel confusion matrix. Default: False
        Returns:
            None
    """
    with open(title, "w") as fo:
        if unilabel:
            fo.write("\t".join(classification_types))
        else:
            fo.write("classification_types\tTP\tFP\tFN\tTN\n")
        for i, row in enumerate(cm):
            fo.write(f"{classification_types[i]}\t")
            fo.write("\t".join(["\t".join(str(j).strip("[]").split()) for j in row]))
            fo.write("\n")
    return None


def decision_tree_visualisation(classifier: RandomForestClassifier) -> None:
    feature_names = get_names()
    export_graphviz(
        classifier.estimators_[0],
        out_file="bsf_decision.dot",
        feature_names=feature_names,
    )


def rf_classify(ids_classifications, fp_folder) -> None:
    n_estimators, max_depth = 1000, 100
    # separate_none = False  # for NPClassifier based ones that have empty predictions

    ids, classifications = ids_classifications[:, 0], ids_classifications[:, 1]

    y, cl_idx = multilabel_and_dict(classifications)
    print(y.shape, len(cl_idx))
    with ChangeDirectory(fp_folder.parent / "output"):
        np.savetxt("ids.tsv", ids, delimiter="\t", fmt="%s")
        np.savetxt("classifications.tsv", classifications, delimiter="\t", fmt="%s")
        np.savetxt(
            "class_labels.tsv", np.array(list(cl_idx.keys())), delimiter="\t", fmt="%s"
        )

    for fp_name, X in _read_fps(fp_folder):
        logging.info(f"{fp_name} X, y:{X.shape}, {y.shape}")
        # k-fold cross validation. ------------------------------------------------------
        classes = list(cl_idx.keys())
        importances, probas, ks = kfold_performance(
            X, y, n_estimators=n_estimators, max_depth=max_depth
        )

        with ChangeDirectory(fp_folder.parent / "output"):
            np.savetxt(f"ks.tsv", ks, delimiter="\t", fmt="%s")
            np.savetxt(f"{fp_name}_proba.tsv", probas, delimiter="\t", fmt="%.3f")
            np.savetxt(
                f"{fp_name}_importances.tsv",
                importances,
                delimiter="\t",
                fmt="%.3f",
            )
            tracemalloc.start()
            tic = time.perf_counter()
            full_classifier = RandomForestClassifier(
                n_estimators=n_estimators, max_depth=max_depth
            ).fit(X, y)
            toc = time.perf_counter()
            current, peak = tracemalloc.get_traced_memory()
            tracemalloc.stop()
            logging.info(f"Current:{current / 10**6}MB; Peak:{peak / 10**6}MB")
            logging.info(f"training took {toc-tic} seconds")

            joblib.dump(full_classifier, f"{fp_name}_model.joblib", compress=3)
            np.savetxt(
                f"{fp_name}_time.tsv",
                np.array([toc - tic]),
                delimiter="\t",
                fmt="%s",
            )
            np.savetxt(
                f"{fp_name}_mem.tsv",
                np.array([peak / 10**6]),
                delimiter="\t",
                fmt="%s",
            )
            np.savetxt("model_labels.tsv", classes, delimiter="\t", fmt="%s")
            if fp_name == "bsf":
                decision_tree_visualisation(full_classifier)
    return None


def write_similarities(fp_folder: np.array):
    for fp_name, fps in _read_fps(fp_folder):
        sims = np.empty((fps.shape[0], fps.shape[0]), dtype="f")  # float32
        with ChangeDirectory(fp_folder.parent / "output"):
            filename = f"{fp_name}_sim.csv"
            if Path(filename).exists():
                continue

            if fp_name == "bsf":
                with open(filename, "w") as f:
                    pass
                for i, fp1 in tqdm(
                    enumerate(fps),
                    total=fps.shape[0],
                    desc="similarity matrix",
                    unit="fp",
                ):
                    i_sims = np.zeros(len(fps), dtype="f")
                    if not all(fp1 == 0):
                        for j, fp2 in enumerate(fps[i:], start=i):
                            if np.sum(fp1 | fp2) == 0:
                                i_sims[j] = -1.0
                            else:
                                i_sims[j] = np.sum(np.minimum(fp1, fp2)) / np.sum(
                                    np.maximum(fp1, fp2)
                                )
                        # append to file
                        # write as %.3f
                    with open(filename, "a") as f:
                        # round to 3 decimal places
                        i_sims = np.round(i_sims, 3)
                        f.write(f"{','.join(map(str, i_sims))}\n")

            else:
                continue
                sims = 1 - squareform(pdist(fps, metric="jaccard"))
                print(sims.shape)
                np.savetxt(filename, sims, delimiter=",", fmt="%.3f")


def dimensionality_reduction(output_folder) -> None:
    for fp_name, sim in _read_sim(output_folder):
        print(fp_name, sim.shape)
        if fp_name == "maccs":
            continue
        # duplicte the upper triangle to the lower triangle
        # make lower triangle 0
        if fp_name == "bsf":
            sim = np.triu(sim)
            sim = np.triu(sim) + np.triu(sim, 1).T
            # fill the diagonal with 1s
            np.fill_diagonal(sim, 1)
            sim = np.where(np.isnan(sim), 0, sim)

        dist = 1 - sim
        dist = np.where(np.isnan(dist), 1, dist)
        dist = np.where(dist < 0, 1, dist)
        assert np.all(dist >= 0), f"negative values: [{dist.min()}, {dist.max()}]"

        with ChangeDirectory(output_folder):
            # make umap
            umap = UMAP(n_components=2, metric="precomputed")
            embeddings = umap.fit_transform(dist, force_all_finite=False)
            np.savetxt(f"{fp_name}_umap.csv", embeddings, delimiter=",", fmt="%.3f")
            print(Path.cwd() / f"{fp_name}_umap.csv")

            # make tsne
            tsne = TSNE(n_components=2, metric="precomputed", init="random")
            embeddings = tsne.fit_transform(dist)
            np.savetxt(f"{fp_name}_tsne.csv", embeddings, delimiter=",", fmt="%.3f")

    return None


def main():
    input_folder = Path(sys.argv[1]).resolve(strict=True)  # root/data/input
    fp_folder = input_folder.parent.parent / "fps"  # root/fps

    ids_classifications = np.loadtxt(
        input_folder / "chebi_classes.csv", delimiter=",", dtype=str
    )
    # rf_classify(ids_classifications, fp_folder)
    # write_similarities(fp_folder)
    dimensionality_reduction(fp_folder.parent / "output")

    exit(0)
    return None


if __name__ == "__main__":
    main()
