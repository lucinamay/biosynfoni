import argparse, os
from collections import Counter
import pickle


import numpy as np
from tqdm import tqdm
from sklearn.metrics import (
    confusion_matrix,
    classification_report,
    accuracy_score,
    multilabel_confusion_matrix,
)
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import KFold, GroupKFold
from sklearn.tree import export_graphviz


def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Multi-label Random forest classifier. "
    )
    parser.add_argument(
        "-f",
        "--fingerprints",
        type=str,
        help="Path to CSV file containing fingerprints.",
    )
    parser.add_argument(
        "-c",
        "--classifications",
        type=str,
        # nargs="+", # one or more
        help="Path to CSV file containing NPClassifier classifications.",
    )
    parser.add_argument(
        "-m",
        "--model",
        type=str,
        help="Path to model.",
    )
    parser.add_argument(
        "-n",
        "--names",
        type=str,
        help="Path to file containing names of compounds.",
    )

    args = parser.parse_args()
    return args


def main():
    args = cli()

    # load data
    fingerprints = np.loadtxt(args.fingerprints, delimiter=",")
    # classifications = np.loadtxt(args.classifications, delimiter=",")
    names = np.loadtxt(args.names, delimiter=",", dtype=str, usecols=0)

    # load model from .npy file
    model = pickle.load(open(args.model, "rb"))

    # predict
    predictions = model.predict(fingerprints)
    probabilites = model.predict_proba(fingerprints)

    # print probabilites
    print(probabilites[:100])

    # print results
    print(predictions[:100])

    # # print results
    # print(classification_report(classifications, predictions, target_names=names))
    # print(confusion_matrix(classifications, predictions))
    # print(accuracy_score(classifications, predictions))

    # predict


if __name__ == "__main__":
    main()
