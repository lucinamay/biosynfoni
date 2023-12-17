#!/usr/bin/env python3
import argparse, os
from collections import Counter


import numpy as np

# import umap
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
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "fingerprints",
        type=str,
        help="Path to CSV file containing fingerprints.",
    )
    parser.add_argument(
        "classifications",
        type=str,
        help="Path to CSV file containing NPClassifier classifications.",
    )
    parser.add_argument(
        "-n",
        "--names",
        "--ids",
        type=str,
        default=None,
        help="Path to file containing names of molecules. Will read first column in as names.",
    )
    parser.add_argument(
        "-s",
        "--subsample",
        type=int,
        default=None,
        help="Subsample size. Default: None., does not subsample.",
    )
    parser.add_argument(
        "-r",
        "--random_seed",
        type=int,
        default=None,
        help="Random seed for reproducibility. Default: None, does not set random seed.",
    )
    parser.add_argument(
        "-u",
        "--unilabel",
        "--singlelabel",
        action="store_true",
        default=False,
        help="If set, will run single label classification instead of multilabel classification.",
    )
    parser.add_argument(
        "-c",
        "--cutoff",
        type=float,
        default=0.5,
        help="Cutoff for classification. Default: 0.5",
    )
    parser.add_argument(
        "--include_none",
        action="store_true",
        default=False,
        help=(
            "If set, will include None class in classification training."
            "Default: False -- will exclude compounds with no classification from the training set, "
            "and will instead use them only as test set for wrong_prediction output"
            "(not for k-fold cross validation))"
        ),
    )
    return parser.parse_args()


def get_unilabels(classifications_path: str) -> tuple[list, dict]:
    """returns list of classification types"""
    # parse labels from input file
    classification_types = Counter()
    labels_names = []
    with open(classifications_path, "r") as fo:
        for line in tqdm(fo):
            line = line.strip().split("\t")  # if more than one type of class
            classification = line[0]
            if classification == "":
                classification = "None"
            if isinstance(classification, list):
                classification = classification.split(",").sort()
                classification = ",".join(classification)
            labels_names.append(classification)
            classification_types[classification] += 1

    class_index = {cl: i for i, cl in enumerate(sorted(classification_types.keys()))}
    if "" in class_index.keys():
        class_index["None"] = class_index.pop("")
    for j in range(len(labels_names)):
        if labels_names[j] == "":
            labels_names[j] = "None"
    return labels_names, class_index


def get_multilabels(classifications_path: str) -> tuple[list, dict]:
    """from the file, parse the classifications for the first tab-separated field,
    and extract the individual labels separated by commas, putting them in a multilabel
    classification table
    """
    classification_types = Counter()
    labels_names = []
    with open(classifications_path, "r") as fo:
        for line in tqdm(fo):
            line = line.strip().split("\t")  # if more than one type of class
            classifications = line[0]
            if classifications == "":
                classifications = "None"
            classifications = classifications.split(",")
            for classification in classifications:
                classification_types[classification] += 1
            labels_names.append(classifications)

    class_index = {cl: i for i, cl in enumerate(sorted(classification_types.keys()))}
    return labels_names, class_index


def unilabel_to_numeric(labels_names: list, class_index: dict) -> tuple:
    """for a list labels, return integer values assigned to each label"""
    # Assign every key from classification type to an integer.
    y = np.array([class_index[label] for label in labels_names])
    return y


def multilabels_to_binary(labels_names: list, class_index: dict) -> tuple:
    """from the list of labels, make a binary table of labels"""
    # Labels.
    y = np.zeros((len(labels_names), len(class_index.keys())), dtype=int)
    for i, labels in enumerate(labels_names):
        for label in labels:
            y[i, class_index[label]] = 1
    return y


def subsampler(
    X: np.array,
    y: np.array,
    size: int = 10000,
    random_seed: int = None,
    verbose: bool = True,
) -> tuple:
    if verbose:
        print(f"Subsampling {args.subsample}. Random seed: {args.random_seed}")
    if isinstance(random_seed, int):
        np.random.seed(random_seed)
    indices = np.random.choice(X.shape[0], size=size, replace=False)
    X = X[indices]
    y = y[indices]
    if verbose:
        print(X.shape, y.shape)
    return X, y


def subsampler_ids(
    X: np.array,
    y: np.array,
    ids: np.array,
    size: int = 10000,
    random_seed: int = None,
    verbose: bool = True,
) -> tuple:
    """subsample including ids"""
    if verbose:
        print(f"Subsampling {args.subsample}. Random seed: {args.random_seed}")

    if len(ids) == 0:
        ids = np.arange(X.shape[0])
    indices = np.random.choice(X.shape[0], size=size, replace=False)
    X = X[indices]
    y = y[indices]
    ids = ids[indices]

    if verbose:
        print(X.shape, y.shape)
    return X, y, ids


def train_test_split(X: np.array, y: np.array, fract: float = 0.8) -> tuple:
    """regular train test split"""
    # Split into train and test.
    indices = np.random.permutation(X.shape[0])
    X = X[indices]
    y = y[indices]
    train_size = int(fract * X.shape[0])
    X_train, X_test = X[:train_size], X[train_size:]
    y_train, y_test = y[:train_size], y[train_size:]
    print(X_train.shape, X_test.shape, y_train.shape, y_test.shape)
    return X_train, X_test, y_train, y_test


def tr_te_split_ids(
    X: np.array, y: np.array, ids: np.array = np.array([]), fract: float = 0.8
) -> tuple:
    """train test split including ids"""
    # Split into train and test.
    if len(ids) == 0:
        ids = np.arange(X.shape[0])

    indices = np.random.permutation(X.shape[0])
    X = X[indices]
    y = y[indices]
    ids = ids[indices]

    train_size = int(fract * X.shape[0])
    X_train, X_test = X[:train_size], X[train_size:]
    y_train, y_test = y[:train_size], y[train_size:]
    ids_train, ids_test = ids[:train_size], ids[train_size:]

    print(
        X_train.shape,
        X_test.shape,
        y_train.shape,
        y_test.shape,
        ids_train.shape,
        ids_test.shape,
    )
    return X_train, X_test, y_train, y_test, ids_train, ids_test


def random_forester(
    X_train: np.array,
    X_test: np.array,
    y_train: np.array,
    y_test: np.array,
    n_estimators: int = 1000,  # number of trees
    max_depth: int = 100,
) -> tuple[list, list]:
    # Train random forest classifier.
    clf = RandomForestClassifier(
        n_estimators=n_estimators, max_depth=max_depth, n_jobs=-1
    )
    clf.fit(X_train, y_train)
    # y_pred = clf.predict_proba(X_test)
    y_pred = clf.predict(X_test)
    return np.array(y_pred), clf.feature_importances_


def rf_proba(
    X_train: np.array,
    X_test: np.array,
    y_train: np.array,
    y_test: np.array,
    n_estimators: int = 1000,  # number of trees
    max_depth: int = 100,
) -> tuple[list, list]:
    # Train random forest classifier.
    clf = RandomForestClassifier(
        n_estimators=n_estimators, max_depth=max_depth, n_jobs=-1
    )
    clf.fit(X_train, y_train)
    probabilities = np.array(clf.predict_proba(X_test))
    # only probability of "yes class"
    only_class_membership = probabilities[:, :, 1]
    return np.transpose(only_class_membership)


def cutoff(y_proba: np.array, cutoff: float = 0.5) -> np.array:
    # change all values below cutoff to 0
    y_proba[y_proba < cutoff] = 0
    # change all values above cutoff to 1
    y_proba[y_proba >= cutoff] = 1
    return y_proba.astype(int)


def kfold_yielder(X: np.array, y: np.array, k: int = 5) -> tuple:
    k = 5
    kf = KFold(n_splits=k, shuffle=True)
    for train_index, test_index in kf.split(X):
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]
        # print(X_train.shape, X_test.shape, y_train.shape, y_test.shape)
        yield X_train, X_test, y_train, y_test


def get_wrong_predictions(
    y_test: np.array, y_pred: np.array, ids_test: np.array
) -> tuple:
    y_test, y_pred, ids_test = np.array(y_test), np.array(y_pred), np.array(ids_test)
    wrong = np.where(y_test != y_pred)
    if isinstance(wrong, tuple):
        wrong = wrong[0]
    wrong_ids = ids_test[wrong]
    return wrong_ids, y_test[wrong], y_pred[wrong]


def write_comb_confusion_matrix(
    cm: np.array,
    classification_types: list,
    title: str = "confusion_matrix.txt",
    unilabel: bool = False,
) -> None:
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


def write_unilab_confusion_matrices(
    cms: list, classification_types: list, title: str = "confusion_matrices.txt"
) -> None:
    with open(title, "w") as fo:
        fo.write(f"\n{classification_types}\n")
        fo.write("\n")
        for cm in cms:
            for i, row in enumerate(cm):
                fo.write("\t".join(str(j) for j in row))
                fo.write("\n")
            fo.write("\n")
    return None


def write_multilab_confusion_matrices(
    cms: list, classification_types: list, title: str = "confusion_matrices.txt"
) -> None:
    with open(title, "w") as fo:
        fo.write(f"\n{classification_types}\n")
        fo.write("\n")
        for cm in cms:
            for i, row in enumerate(cm):
                fo.write(f"\n{str(row)}")
                fo.write("\n")
            fo.write("\n")
    return None


def combine_cms(cms: list[np.ndarray]) -> np.ndarray:
    combined = np.zeros(cms[0].shape)
    for cm in cms:
        combined += cm
    return combined.astype(int)


def write_classification_report(
    classification_results: list, classes: list, title="classification_report.txt"
) -> None:
    with open(title, "w") as fo:
        fo.write("\n")
        fo.write("\t".join(classes))
        fo.write("\n")
        for report in classification_results:
            fo.write("\n")
            for line in report.split("\n"):
                fo.write(f"{line}\n")
    return None


def handle_outdirs(db_name: str, fp_name: str, unilabel: bool) -> str:
    iwd = os.getcwd()
    if not os.path.exists("./rf"):
        os.mkdir("./rf")
    os.chdir("./rf")
    if not os.path.exists(f"./{db_name}"):
        os.mkdir(f"./{db_name}")
    os.chdir(f"./{db_name}")
    if not os.path.exists(f"./{fp_name}"):
        os.mkdir(f"./{fp_name}")
    os.chdir(f"./{fp_name}")
    if unilabel:
        if not os.path.exists("./unilabel"):
            os.mkdir("./unilabel")
        os.chdir("./unilabel")
    return iwd


def main() -> None:
    n_estimators = 1000  # number of trees
    max_depth = 100
    args = cli()

    # get absolute path from relative path
    fp_path = os.path.abspath(args.fingerprints)
    args.classifications = os.path.abspath(args.classifications)
    if args.names:
        args.names = os.path.abspath(args.names)

    folder = os.path.basename(os.path.dirname(fp_path))
    # get fingerprint name
    fp_name = fp_path.split("/")[-1].split(".")[0]
    # get database name
    db_name = folder.split("/")[-1]

    # make directory for output
    iwd = handle_outdirs(db_name, fp_name, args.unilabel)

    if args.unilabel and args.cutoff != 0.5:
        print(f"WARNING: cutoff {args.cutoff} is ignored for unilabel classification.")

    # Parse fingerprints from input file.
    delimiter = "\t" if args.fingerprints.endswith(".tsv") else ","
    X = np.loadtxt(fp_path, delimiter=delimiter, dtype=int)

    if args.unilabel:
        # Parse labels from input file.
        label_names, class_index = get_unilabels(args.classifications)
        # change empty key to "None"
        # Labels.
        # Assign every key from classification type to an integer.
        y = unilabel_to_numeric(label_names, class_index)
    else:
        label_names, class_index = get_multilabels(args.classifications)
        y = multilabels_to_binary(label_names, class_index)
    print("X, y:", X.shape, y.shape)

    if args.names:
        ids = np.loadtxt(args.names, dtype=str, delimiter="\t", usecols=0)

    else:
        ids = np.array(range(X.shape[0]))

    # print(label_names[:100])
    # exit()

    none_indices = []
    if not args.include_none:
        for i, name in enumerate(label_names):
            if name == ["None"]:
                none_indices.append(int(i))
        if "None" in class_index.keys():
            none_index = class_index.pop("None")
            if not args.unilabel:
                y = np.delete(y, none_index, 1)  # drop column
        print(f'separating {len(none_indices)} "nones" from training/testing set')

    none_indices = np.array(none_indices, dtype=int)
    nones_X, nones_y, nones_ids = X[none_indices], y[none_indices], ids[none_indices]

    X = np.delete(X, none_indices, 0)
    y = np.delete(y, none_indices, 0)
    ids = np.delete(ids, none_indices, 0)

    classes = list(class_index.keys())

    # Subsample 10.000 randomly from X and y.
    if isinstance(args.subsample, int):
        X, y, ids = subsampler_ids(
            X, y, ids, size=args.subsample, random_seed=args.random_seed
        )

    # Split into train and test.
    # X_train, X_test, y_train, y_test = train_test_split(X, y, fract=0.8)
    X_train, X_test, y_train, y_test, ids_train, ids_test = tr_te_split_ids(
        X, y, ids, fract=0.8
    )
    # print(y_train, y_test)

    print(
        "starting one RF to get importances, wrong predictions and/or predictions for Nones"
    )
    if args.unilabel:
        y_pred, _ = random_forester(
            X_train,
            X_test,
            y_train,
            y_test,
            n_estimators=n_estimators,
            max_depth=max_depth,
        )
        cm = confusion_matrix(y_test, y_pred)

        if len(nones_X) > 0:
            # get prediction for the "Nones"
            ny_pred, _ = random_forester(
                X,
                nones_X,
                y,
                nones_y,
                n_estimators=n_estimators,
                max_depth=max_depth,
            )

    else:
        # multilabel, consider cutoff, no feature importances
        y_proba = rf_proba(
            X_train,
            X_test,
            y_train,
            y_test,
            n_estimators=n_estimators,
            max_depth=max_depth,
        )
        # test = np.transpose(y_proba)
        # print(y_proba, args.cutoff)
        print(y_proba.shape)
        y_pred = cutoff(
            y_proba,
            cutoff=args.cutoff,
        )
        if len(nones_X) > 0:
            print('get prediction for the "Nones"')
            ny_proba = rf_proba(
                X,
                nones_X,
                y,
                nones_y,
                n_estimators=n_estimators,
                max_depth=max_depth,
            )
            ny_pred = cutoff(
                ny_proba,
                cutoff=args.cutoff,
            )

        print(y_test.shape, y_pred.shape)
        cm = multilabel_confusion_matrix(y_test, y_pred)

    if len(nones_X) > 0:
        ny_pred = [class_array[np.where(i == 1)] for i in ny_pred]
        ny_pred = [",".join(i) for i in ny_pred]
        nones_preds = np.concatenate([nones_ids, ny_pred], axis=1)
        np.savetxt("nones_predictions.tsv", nones_preds, delimiter="\t")

        nones_probas = np.concatenate([nones_ids, ny_proba], axis=1)
        np.savetxt("nones_probas.tsv", nones_probas, delimiter="\t")

    try:
        wrong_ids, wrongs_test, wrongs_pred = get_wrong_predictions(
            y_test, y_pred, ids_test
        )
        # for each row, get indexes of column with value 1:
        class_array = np.array(classes)
        wrongs_test_names = [class_array[np.where(i == 1)] for i in wrongs_test]
        wrongs_pred_names = [class_array[np.where(i == 1)] for i in wrongs_pred]

        wrongs_test_names = [",".join(i) for i in wrongs_test_names]
        wrongs_pred_names = [",".join(i) for i in wrongs_pred_names]

        # write to tsv
        np.savetxt(
            "wrong_predictions.tsv",
            np.array([wrong_ids, wrongs_test_names, wrongs_pred_names]).T,
            delimiter="\t",
            fmt="%s",
        )
    except:
        print("wrong predictions not (/not properly) detected")

    # Perform k-fold cross validation.
    confusion_matrices = []
    classification_results = []
    importances = []
    k = 5
    for X_train, X_test, y_train, y_test in tqdm(kfold_yielder(X, y, k=k), total=k):
        if args.unilabel:
            y_pred, importance = random_forester(
                X_train,
                X_test,
                y_train,
                y_test,
                n_estimators=n_estimators,
                max_depth=max_depth,
            )
            # print(f"importances:\n{importance}")
            classification_results.append(
                classification_report(y_test, y_pred, zero_division=np.nan)
            )
            cm = confusion_matrix(y_test, y_pred)

        else:
            # multilabel, consider cutoff
            y_proba = rf_proba(
                X_train,
                X_test,
                y_train,
                y_test,
                n_estimators=n_estimators,
                max_depth=max_depth,
            )
            y_pred = cutoff(y_proba, cutoff=args.cutoff)
            classification_results.append(
                classification_report(
                    y_test, y_pred, target_names=classes, zero_division=np.nan
                )
            )
            cm = multilabel_confusion_matrix(y_test, y_pred)
            importance = []

        # Get individual confusion matrix.
        # print(f"confusion_matrix:\n{cm}")
        # append
        confusion_matrices.append(cm)
        importances.append(importance)

    write_classification_report(classification_results, classes)
    write_comb_confusion_matrix(combine_cms(confusion_matrices), classes)

    # Write importances.
    np.savetxt(
        "importances.tsv", np.array(importances), delimiter="\t", fmt="%s"
    )  # graph @ fp_avger

    # get argument dictionary
    args_dict = vars(args)
    # Write arguments as yaml.
    with open("arguments.yaml", "w") as fo:
        for key, value in args_dict.items():
            fo.write(f"{key}: {value}\n")

    # Write decision tree.
    # export_graphviz(clf.estimators_[0], out_file="tree.dot", feature_names=classes)

    # Write wrong predictions.

    # for each category in combined_cms, write to file
    # for cat in range(combined_cms.shape[0]):
    #     np.savetxt(
    #         f"combined_confusion_matrix_{classes[cat]}.tsv",
    #         combined_cms[cat],
    #         delimiter="\t",
    #         fmt="%s",
    #     )

    os.chdir(iwd)
    exit(0)

    # Get confusion matrices.

    # # for each category in combined_cms, write to file
    # for cat in range(combined_cms.shape[0]):
    #     np.savetxt(
    #         f"combined_confusion_matrix_{classes[cat]}.tsv",
    #         combined_cms[cat],
    #         delimiter="\t",
    #         fmt="%s",
    #     )

    os.chdir(cwd)
    exit(0)


if __name__ == "__main__":
    main()
