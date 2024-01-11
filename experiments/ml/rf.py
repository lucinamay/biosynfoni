#!/usr/bin/env python3
import argparse, os
from collections import Counter
import pickle

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
    parser = argparse.ArgumentParser(
        description="Multi-label Random forest classifier. "
    )
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
    parser.add_argument(
        "-e",
        "--export",
        action="store_true",
        help="If set, will export model as pickle file.",
    )
    args = parser.parse_args()

    if args.unilabel and args.cutoff != 0.5:
        print(f"WARNING: cutoff {args.cutoff} is ignored for unilabel classification.")

    return args


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


def read_in_classifications(cl_path: str, unilabel: bool) -> tuple:
    if unilabel:
        label_names, class_index = get_unilabels(cl_path)
        y = unilabel_to_numeric(label_names, class_index)
    else:
        label_names, class_index = get_multilabels(cl_path)
        y = multilabels_to_binary(label_names, class_index)
    return label_names, class_index, y


def i_to_cl(y: np.array, class_index: dict, unilabel: bool = False) -> tuple:
    """for each row, get indexes of column with value 1"""
    if unilabel:
        # single prediction, with integer corresponding to class
        ind_cl = {i: cl for cl, i in class_index.items()}
        y_classes = np.array([ind_cl[i] for i in y])
    else:
        # multilabel prediction, with binary values in an array
        # re-sort classes on index values in case something changed
        classes = [k for k, v in sorted(class_index.items(), key=lambda item: item[1])]
        class_array = np.array(classes)

        y_classes = [class_array[np.where(i == 1)] for i in y]
        y_classes = np.array([",".join(i) for i in y_classes])
    return y_classes


def separate_nones(label_names, class_index, X, y, ids, unilabel=False):
    none_indices = []
    # get indexes of "None" annotations
    for i, name in enumerate(label_names):
        if name == ["None"]:
            none_indices.append(int(i))

    # remove "None" class index dictionary (integer values)
    if "None" in class_index.keys():
        none_index = class_index.pop("None")
        if not unilabel:
            y = np.delete(y, none_index, 1)  # drop column
    print(f'separating {len(none_indices)} "nones" from training/testing set')

    none_indices = np.array(none_indices, dtype=int)
    nones_X, nones_y, nones_ids = X[none_indices], y[none_indices], ids[none_indices]

    X = np.delete(X, none_indices, 0)
    y = np.delete(y, none_indices, 0)
    ids = np.delete(ids, none_indices, 0)
    return label_names, class_index, X, y, ids, nones_X, nones_y, nones_ids


def subsampler(
    X: np.array,
    y: np.array,
    size: int = 10000,
    random_seed: int = None,
    verbose: bool = True,
) -> tuple:
    if verbose:
        print(f"Subsampling {size}. Random seed: {random_seed}")
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
        print(f"Subsampling {size}. Random seed: {random_seed}")

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
    """predict probabilities with RF"""
    # Train random forest classifier.
    clf = RandomForestClassifier(
        n_estimators=n_estimators, max_depth=max_depth, n_jobs=-1
    )
    clf.fit(X_train, y_train)
    probabilities = np.array(clf.predict_proba(X_test))
    # only probability of "yes class"
    only_class_membership = probabilities[:, :, 1]
    return np.transpose(only_class_membership)


def cutoffr(y_proba: np.array, cutoff: float = 0.5) -> np.array:
    """bitifies the probabilities with given cutoff (>= cutoff = 1, < cutoff = 0)"""
    # make copy of y_proba to not affect the probability array
    y_rounded = np.copy(y_proba)
    # change all values below cutoff to 0 and above/equal cutoff to 1
    y_rounded[y_rounded < cutoff] = 0
    y_rounded[y_rounded >= cutoff] = 1
    return y_rounded.astype(int)


def kfold_yielder(X: np.array, y: np.array, k: int = 5) -> tuple:
    """splits X and y into k folds and yields train and test sets for each fold"""
    k = 5
    kf = KFold(n_splits=k, shuffle=True)
    for train_index, test_index in kf.split(X):
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]
        # print(X_train.shape, X_test.shape, y_train.shape, y_test.shape)
        yield X_train, X_test, y_train, y_test


def kfold_yielder_ids(X: np.array, y: np.array, ids: np.array, k: int = 5) -> tuple:
    """splits X, y and ids into k folds and yields train and test sets for each fold"""
    k = 5
    kf = KFold(n_splits=k, shuffle=True)
    for train_index, test_index in kf.split(X):
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]
        ids_train, ids_test = ids[train_index], ids[test_index]
        # print(X_train.shape, X_test.shape, y_train.shape, y_test.shape)
        yield X_train, X_test, y_train, y_test, ids_train, ids_test


def kfold_preds(
    X,
    y,
    k: int = 5,
    unilabel: bool = False,
    n_estimators: int = 1000,
    max_depth: int = 100,
    cutoff: float = 0.5,
    target_names=None,
):
    # Perform k-fold cross validation.
    cms = []  # confusion matrices
    cl_reps = []  # classification reports
    importances = []
    for X_train, X_test, y_train, y_test in tqdm(kfold_yielder(X, y, k=k), total=k):
        if unilabel:
            y_pred, importance = random_forester(
                X_train,
                X_test,
                y_train,
                y_test,
                n_estimators=n_estimators,
                max_depth=max_depth,
            )
            cl_rep = classification_report(y_test, y_pred, zero_division=np.nan)
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
            y_pred = cutoffr(y_proba, cutoff=cutoff)
            # classification report with classes as target names for multilabel
            cl_rep = classification_report(
                y_test, y_pred, target_names=target_names, zero_division=np.nan
            )
            cm = multilabel_confusion_matrix(y_test, y_pred)
            # feature importance not possible for predict_proba,
            # therefore not included if multilabel (i.e. proba is possible)
            importance = []
        # append results
        cms.append(cm)
        cl_reps.append(cl_rep)
        importances.append(importance)
    return cms, cl_reps, importances


def get_wrong_predictions(
    y_test: np.array, y_pred: np.array, ids_test: np.array
) -> tuple:
    y_test, y_pred, ids_test = np.array(y_test), np.array(y_pred), np.array(ids_test)
    wrong = np.where(y_test != y_pred)

    if isinstance(wrong, tuple):
        # as the np.where will return two arrays (if multilabel?)
        wrong = wrong[0]
    # get array with only unique wrong
    wrong = np.unique(wrong)
    return ids_test[wrong], y_test[wrong], y_pred[wrong]


def wrongs_array(
    ids_test: np.array, y_test: np.array, y_pred: np.array, class_index: dict
) -> np.array:
    # get wrong predictions.
    w_id, w_test, w_pred = get_wrong_predictions(y_test, y_pred, ids_test)
    # get the classification predictions in strings:
    w_test_str = i_to_cl(w_test, class_index)
    w_pred_str = i_to_cl(w_pred, class_index)
    wrong = np.array([w_id, w_test_str, w_pred_str]).T
    print(wrong.shape[0], " wrong predictions")
    return wrong


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


def args_yaml(args: argparse.Namespace, title: str = "arguments.yaml") -> None:
    args_dict = vars(args)
    # Write arguments as yaml.
    with open(title, "w") as fo:
        for key, value in args_dict.items():
            fo.write(f"{key}: {value}\n")
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
    print(args)
    include_none = args.include_none  # preserve value of args.include_none

    # get absolute path from relative path
    fp_path = os.path.abspath(args.fingerprints)
    classif = os.path.abspath(args.classifications)
    if args.names:
        args.names = os.path.abspath(args.names)

    folder = os.path.basename(os.path.dirname(fp_path))
    # get fingerprint name
    fp_name = fp_path.split("/")[-1].split(".")[0]
    # get database name
    db_name = folder.split("/")[-1]

    # make directory for output
    iwd = handle_outdirs(db_name, fp_name, args.unilabel)

    # Parse fingerprints from input file.
    delimiter = "\t" if args.fingerprints.endswith(".tsv") else ","
    X = np.loadtxt(fp_path, delimiter=delimiter, dtype=int)

    label_names, class_index, y = read_in_classifications(classif, args.unilabel)
    print("X, y:", X.shape, y.shape)

    if args.names:
        ids = np.loadtxt(args.names, dtype=str, delimiter="\t", usecols=0)
    else:
        ids = np.array(range(X.shape[0]))

    if not include_none:
        (
            label_names,
            class_index,
            X,
            y,
            ids,
            nones_X,
            nones_y,
            nones_ids,
        ) = separate_nones(label_names, class_index, X, y, ids, unilabel=False)

    if nones_X.shape[0] == 0:
        print("no nones")
        include_none = True
    classes = list(class_index.keys())

    # Subsample 10.000 randomly from X and y.
    if isinstance(args.subsample, int):
        X, y, ids = subsampler_ids(
            X, y, ids, size=args.subsample, random_seed=args.random_seed
        )

    # Split into train and test.
    X_train, X_test, y_train, y_test, ids_train, ids_test = tr_te_split_ids(
        X, y, ids, fract=0.8
    )

    # get initial predictions. ------------------------------------------------------
    print("RF for getting wrong predictions")
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
    else:
        # multilabel, considers cutoff, no feature importances
        y_proba = rf_proba(
            X_train,
            X_test,
            y_train,
            y_test,
            n_estimators=n_estimators,
            max_depth=max_depth,
        )
        y_pred = cutoffr(y_proba, cutoff=args.cutoff)
        cm = multilabel_confusion_matrix(y_test, y_pred)
        # print("proba, test, pred shapes:", y_proba.shape, y_test.shape, y_pred.shape)

    # try and get wrong predictions. ------------------------------------------------
    wrongs = wrongs_array(ids_test, y_test, y_pred, class_index)
    np.savetxt("wrong_predictions.tsv", wrongs, delimiter="\t", fmt="%s")
    # print("wrong predictions not (/not properly) detected")

    # none predictions. -------------------------------------------------------------
    if not include_none:
        print('get prediction for the "Nones"')
        if args.unilabel:
            # get prediction for the "Nones" by training with full not-none dataset
            ny_pred, _ = random_forester(
                X,
                nones_X,
                y,
                nones_y,
                n_estimators=n_estimators,
                max_depth=max_depth,
            )
        else:
            ny_proba = rf_proba(
                X,
                nones_X,
                y,
                nones_y,
                n_estimators=n_estimators,
                max_depth=max_depth,
            )
            ny_pred = cutoffr(ny_proba, cutoff=args.cutoff)

        ny_pred = i_to_cl(ny_pred, class_index, unilabel=args.unilabel)
        nones_preds = np.concatenate([nones_ids[:, None], ny_pred[:, None]], axis=1)
        np.savetxt("nones_predictions.tsv", nones_preds, delimiter="\t", fmt="%s")
        if not args.unilabel:
            nones_probas = np.concatenate([nones_ids[:, np.newaxis], ny_proba], axis=1)
            np.savetxt("nones_probas.tsv", nones_probas, delimiter="\t", fmt="%s")

    # k-fold cross validation. ------------------------------------------------------
    cms, cl_reps, importances = kfold_preds(
        X,
        y,
        unilabel=args.unilabel,
        n_estimators=n_estimators,
        max_depth=max_depth,
        cutoff=args.cutoff,
        target_names=classes,
    )

    write_classification_report(cl_reps, classes)
    write_comb_confusion_matrix(combine_cms(cms), classes)
    np.savetxt("importances.tsv", importances, delimiter="\t", fmt="%s")

    # write arguments as yaml. ------------------------------------------------------
    args_yaml(args, title="arguments.yaml")

    # save model
    if args.export:
        # Train random forest classifier.
        clf = RandomForestClassifier(
            n_estimators=n_estimators, max_depth=max_depth, n_jobs=-1
        )
        # clf.fit(X_train, y_train)
        clf.fit(X, y)
        # save model
        pickle.dump(clf, open("model.pkl", "wb"))
        # save labels for indexes
        np.savetxt("model_labels.tsv", classes, delimiter="\t", fmt="%s")

    # return to initial working directory
    os.chdir(iwd)
    exit(0)
    return None


if __name__ == "__main__":
    main()
