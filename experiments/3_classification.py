#!/usr/bin/env python3
import logging, sys, re, time, tracemalloc
from typing import Generator
from pathlib import Path


import joblib
import numpy as np
from scipy.spatial.distance import pdist, squareform
from sklearn.ensemble import RandomForestClassifier
from sklearn.manifold import TSNE
from sklearn.model_selection import KFold
from sklearn.tree import export_graphviz
from tqdm import tqdm
from umap import UMAP


from helper import ChangeDirectory
from biosynfoni.subkeys import get_names


class TrackMemoryAndTime:
    def __init__(self, name):
        self.fp_name = name

    def __enter__(self):
        tracemalloc.start()
        self.tic = time.perf_counter()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        toc = time.perf_counter()
        current, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()

        # Logging memory and time statistics
        logging.info(f"current:{current / (1024*1024)}MB; peak:{peak / (1024*1024)}MB")
        logging.info(f"training took {toc-self.tic} seconds")

        np.savetxt(f"{self.fp_name}_time.csv", np.array([toc - self.tic]))
        np.savetxt(f"{self.fp_name}_mem.csv", np.array([peak / (1024 * 1024)]))


def _read_fps(folder_path, data="chebi") -> Generator:
    for file in folder_path.glob("*.csv"):
        if (
            f"{data}_" in file.stem
            and not len(file.stem.split("_")) > 2
            and not "coverages" in file.stem
        ):
            yield file.stem.replace(f"{data}_", ""), np.loadtxt(
                file, delimiter=",", dtype=int
            )


def _read_sim(folder_path) -> Generator:
    for file in folder_path.glob("*_sim.csv"):
        if "mini" in file.stem:
            continue
        yield file.stem.replace("chebi_", "").replace("_sim", ""), np.loadtxt(
            file, delimiter=",", dtype=float
        )


def get_ids_and_labels(labels_path) -> np.array:
    ids_classifications = np.loadtxt(labels_path, delimiter=",", dtype=str)
    # ids_classifications[:, 1] = np.array(
    #     [
    #         x.replace("isoprenoid;fatty_acid", "isoprenoid")
    #         for x in ids_classifications[:, 1]
    #     ]
    # )
    return ids_classifications


def multilabel_and_dict(classifications: np.array) -> tuple[np.array, dict]:
    """
    Convert classifications to a binary array and a dictionary of classes to indices
    """
    classes = set(";".join(map(str, classifications)).split(";"))
    class_to_id = {class_: i for i, class_ in enumerate(sorted(classes))}
    id_to_class = {i: class_ for class_, i in class_to_id.items()}
    classification_array = np.zeros((len(classifications), len(classes)), dtype=int)
    for i, classification in enumerate(classifications):
        for class_ in re.split(r"[;,]", classification):
            classification_array[i, class_to_id[class_]] = 1
    return classification_array, id_to_class


def predict_one(X_test: np.array, classifier) -> np.array:
    return np.array(classifier.predict(X_test))


def kfold_performance(X, y, *args, **kwargs) -> dict:
    """
    *args and **kwargs are passed to RandomForestClassifier
    """
    importances = []
    k = 5
    # train_probabilities, train_ks = (np.empty(y.shape, k),np.empty(y.shape[0], k))
    train_probabilities = np.empty((y.shape[0], y.shape[1], k))
    probabilities, ks = (np.empty(y.shape), np.empty(y.shape[0]))

    for i, (train, test) in enumerate(KFold(shuffle=True, random_state=42).split(X)):
        X_train, X_test = X[train], X[test]
        y_train, y_test = y[train], y[test]
        classifier = RandomForestClassifier(*args, **kwargs).fit(X_train, y_train)
        importances.append(classifier.feature_importances_)

        # only probabilities of being in class
        y_train_proba = np.transpose(
            np.array(classifier.predict_proba(X_train))[:, :, 1]
        )
        y_proba = np.transpose(np.array(classifier.predict_proba(X_test))[:, :, 1])
        train_probabilities[train, :, i] = y_train_proba
        # train_ks[train,i] = i
        probabilities[test] = y_proba
        ks[test] = i

    return {
        "importances": importances,
        "probabilities": probabilities,
        "ks": ks,
        "train_probabilities": train_probabilities,
        # "train_ks": train_ks,
    }
    return importances, probabilities, ks


def cluster_performace(X, y, cluster_labels, *args, **kwargs) -> tuple:
    # assess performance based on cluster stratification
    importances = []
    probabilities = []
    ks = cluster_labels
    # do a scaffold cluster based split
    for cluster in np.unique(cluster_labels):
        train = np.where(cluster_labels != cluster)[0]
        test = np.where(cluster_labels == cluster)[0]
        X_train, X_test = X[train], X[test]
        y_train, y_test = y[train], y[test]
        classifier = RandomForestClassifier(*args, **kwargs).fit(X_train, y_train)
        importances.append(classifier.feature_importances_)
        y_proba = np.transpose(np.array(classifier.predict_proba(X_test))[:, :, 1])
        probabilities[test] = y_proba
        ks.append(cluster)
    return importances, probabilities, ks


def export_tree(classifier: RandomForestClassifier) -> None:
    feature_names = get_names()
    export_graphviz(
        classifier.estimators_[0],
        out_file="bsf_decision.dot",
        feature_names=feature_names,
    )


def undersample(multilabels) -> np.array:
    """
    Undersamples labels to the smallest label class size, ignoring and keeping multilabeled instances.
    """
    np.random.seed(42)
    single_labels = set(";".join(map(str, multilabels)).split(";"))
    multilabel_indices = np.nonzero(";" in multilabels)[0]

    # get the indices of each instance of a single label
    label_idxs = {lbl: np.nonzero(multilabels == lbl)[0] for lbl in single_labels}

    # make sure that none of the indices are in the other classes
    for label in single_labels:
        assert all(
            [
                label not in multilabels[label_idxs[label_]]
                for label_ in single_labels
                if label_ != label
            ]
        )

    minsize = min([len(label_idxs[lbl]) for lbl in single_labels])
    keep = np.zeros(multilabels.shape[0], dtype=bool)  # initiate
    keep[multilabel_indices] = True
    for label, idxs in label_idxs.items():
        if len(idxs) == minsize:
            keep[idxs] = True
        else:
            # randomly choose n=minsize instances indices
            keeping_inds = np.random.choice(idxs, minsize, replace=False)
            assert len(keeping_inds) == minsize
            keep[keeping_inds] = True
    return keep


def random_forest(ids_classifications, fp_folder) -> None:
    parameters = {
        "n_estimators": 1000,
        "max_depth": 100,
    }
    ids, classifications = ids_classifications.T
    y, cl_idx = multilabel_and_dict(classifications)

    logging.info(y.shape, len(cl_idx))

    with ChangeDirectory(fp_folder.parent / "output"):
        np.savetxt("ids.tsv", ids, delimiter="\t", fmt="%s")
        np.savetxt("classifications.tsv", classifications, delimiter="\t", fmt="%s")
        np.savetxt("class_labels.tsv", np.array(list(cl_idx)), delimiter="\t", fmt="%s")

    for fp_name, X in _read_fps(fp_folder):
        logging.info(f"{fp_name} X, y:{X.shape}, {y.shape}")

        # k-fold cross validation. ------------------------------------------------------
        # importances, probas, ks = kfold_performance(X, y, **parameters)
        res = kfold_performance(X, y, **parameters)
        importances, probas, ks = res["importances"], res["probabilities"], res["ks"]
        train_probas = res["train_probabilities"]

        with ChangeDirectory(fp_folder.parent / "output"):
            np.savetxt("ks.csv", ks, fmt="%d")
            np.savetxt(f"{fp_name}_proba.tsv", probas, delimiter="\t", fmt="%.3f")
            np.savetxt(
                f"{fp_name}_importances.tsv",
                importances,
                delimiter="\t",
                fmt="%.3f",
            )
            np.save(f"{fp_name}_trainproba.npy", train_probas)

            with TrackMemoryAndTime(fp_name):
                full_classifier = RandomForestClassifier(**parameters).fit(X, y)

            joblib.dump(full_classifier, f"{fp_name}_model.joblib", compress=3)

            np.savetxt(
                "model_labels.tsv", list(cl_idx.keys()), delimiter="\t", fmt="%s"
            )
            if fp_name == "bsf":
                export_tree(full_classifier)
    return


def write_similarities(fp_folder: np.array):
    for fp_name, fps in _read_fps(fp_folder):
        sims = np.empty((fps.shape[0], fps.shape[0]), dtype="f")
        with ChangeDirectory(fp_folder.parent / "output"):
            filename = f"{fp_name}_sim.csv"

            if fp_name == "bsf":
                with open(filename, "w") as f:
                    pass
                for i, fp1 in tqdm(
                    enumerate(fps),
                    # total=fps.shape[0],
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
                    with open(filename, "a") as f:
                        i_sims = np.round(i_sims, 3)
                        f.write(f"{','.join(map(str, i_sims))}\n")

            else:
                sims = 1 - squareform(pdist(fps, metric="jaccard"))
                np.savetxt(filename, sims, delimiter=",", fmt="%.3f")


def dimensionality_reduction(output_folder) -> None:
    for fp_name, sim in _read_sim(output_folder):
        # duplicte the upper triangle to the lower triangle

        sim = np.triu(sim)
        sim = np.triu(sim) + np.triu(sim, 1).T
        # fill the diagonal with 1s
        np.fill_diagonal(sim, 1)
        sim = np.where(np.isnan(sim), 0, sim)  # any nan --> 'no similarity'

        dist = 1 - sim
        dist = np.where(np.isnan(dist), 1, dist)
        dist = np.where(dist < 0, 1, dist)
        assert np.all(dist >= 0), f"negative values: [{dist.min()}, {dist.max()}]"

        with ChangeDirectory(output_folder):
            # make umap
            umap = UMAP(n_components=2, metric="precomputed")
            embeddings = umap.fit_transform(dist)
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

    ids_labels = get_ids_and_labels(input_folder / "chebi_classes.csv")

    random_forest(ids_labels, fp_folder)
    write_similarities(fp_folder)
    dimensionality_reduction(fp_folder.parent / "output")

    exit(0)
    return None


if __name__ == "__main__":
    main()
