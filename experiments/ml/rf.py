#!/usr/bin/env python3
import argparse, os, logging, pickle, time, tracemalloc
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
        "--no_proba",
        action="store_true",
        default=False,
        help="If set, will not calculate probabilities for multilabel classification. Probabilties are not calculated for unilabel classification.",
    )
    parser.add_argument(
        "-e",
        "--export",
        action="store_true",
        default=False,
        help="If set, will export the model ",
    )

    args = parser.parse_args()
    if args.unilabel and args.cutoff != 0.5:
        logging.warning(f"cutoff {args.cutoff} is ignored for unilabel classification.")

    if args.unilabel:
        args.no_proba = True

    # get absolute path from relative path
    args.fingerprints = os.path.abspath(args.fingerprints)
    args.classifications = os.path.abspath(args.classifications)
    if args.names:
        args.names = os.path.abspath(args.names)
    return args


def _read_unilabels(classifications_path: str) -> tuple[list, dict]:
    """
    From the file, parse the classifications for the first tab-separated field and return list of classification types

        Args:
            classifications_path (str): path to file with classifications

        Returns:
            tuple[list, dict]: list of classification types and dictionary of classification types
    """
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


def _read_multilabels(classifications_path: str) -> tuple[list, dict]:
    """
    From the file, parse the classifications for the first tab-separated field and return list of classification types and dictionary of classification types

        Args:
            classifications_path (str): path to file with classifications

        Returns:
            tuple[list, dict]: list of classification types and dictionary of classification types


    extracts the individual labels separated by commas, putting them in a multilabel classification table
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


def _unilabel_to_numeric(labels_names: list, class_index: dict) -> np.array:
    """
    For a list of labels, return integer values assigned to each label

        Args:
            labels_names (list): list of labels
            class_index (dict): dictionary of classification types
        Returns:
            np.array: array of integer values assigned to each label
    """
    return np.array([class_index[label] for label in labels_names])


def _multilabels_to_binary(labels_names: list, class_index: dict) -> np.array:
    """
    From the list of labels, make a binary table of labels

        Args:
            labels_names (list): list of labels
            class_index (dict): dictionary of classification types
        Returns:
            np.array: binary table of labels
    """
    # Labels.
    y = np.zeros((len(labels_names), len(class_index.keys())), dtype=int)
    for i, labels in enumerate(labels_names):
        for label in labels:
            y[i, class_index[label]] = 1
    return y


def read_classifications(cl_path: str, unilabel: bool) -> tuple:
    """
    Read classifications from file and return label names, class index and y

        Args:
            cl_path (str): path to file with classifications
            unilabel (bool): if True, will read unilabels, if False, will read multilabels

        Returns:
            tuple: label names, class index, y
    """
    if unilabel:
        label_names, class_index = _read_unilabels(cl_path)
        y = _unilabel_to_numeric(label_names, class_index)
    else:
        label_names, class_index = _read_multilabels(cl_path)
        y = _multilabels_to_binary(label_names, class_index)
    return label_names, class_index, y


def i_to_cl(y: np.array, class_index: dict, unilabel: bool = False) -> tuple:
    """
    For an array of integer values, return the corresponding classification names

        Args:
            y (np.array): array of integer values (or binary values if multilabel)
            class_index (dict): dictionary of classification types
            unilabel (bool): if True, will return single classification, if False, will return multiple classifications. Default: False
        Returns:
            tuple: array of classification names
    """
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


def get_ids(names_path: str, size: int) -> np.array:
    """
    Get ids from file or create ids from size

        Args:
            names_path (str): path to file with ids
            size (int): size of fingerprints
        Returns:
            np.array: array of ids
    """
    if names_path:
        ids = np.loadtxt(names_path, dtype=str, delimiter="\t", usecols=0)
    else:
        ids = np.array(range(size))

    if ids.shape[0] != size:
        logging.warning(
            f"number of ids ({ids.shape[0]}) does not match number of fingerprints ({size})"
        )

    return ids


def separate_nones(label_names: list[list], class_index: dict, x_y_andco: list):
    """
    Separate "None" classifications from the data

        Args:
            label_names (list[list]): list of classification names
            class_index (dict): dictionary of classification types
            x_y_andco (list): list of arrays to separate
        Returns:
            tuple: label names, class index,

    * x_y_andco: (X, y, ids, ...) (...-> anything else that needs to be separated into nones set)
    """
    none_indices = []
    for i, name in enumerate(label_names):
        if name == ["None"]:
            none_indices.append(int(i))

    x_y_andco = list(x_y_andco)
    y = x_y_andco[1]  #!!! make sure it is a y
    if "None" in class_index.keys():
        none_column = class_index.pop("None")
        # for all numbers higher than none_column, subtract 1
        for key in class_index.keys():
            if class_index[key] > none_column:
                class_index[key] -= 1
        # if not unilabel:
        if y.shape[1] > 1:
            logging.debug(f"removing column {none_column} from y")
            x_y_andco[1] = np.delete(x_y_andco[1], none_column, 1)
    logging.info(f'separating {len(none_indices)} "nones" from training/testing set')

    none_indices = np.array(none_indices, dtype=int)

    nones = []
    for array in x_y_andco:
        nones.append(array[none_indices])  # add to nones
        np.delete(array, none_indices, 0)  # remove from not-nones

    logging.debug([array.shape for array in nones])
    return label_names, class_index, x_y_andco, nones


def subsampler(
    arrays_to_subsample: tuple[np.array],
    size: int = 10000,
    random_seed: int = None,
) -> tuple:
    """
    Subsample arrays in list, keeping same indices

        Args:
            arrays_to_subsample (tuple[np.array]): tuple of arrays to subsample
            size (int): size of subsample
            random_seed (int): random seed for reproducibility
        Returns:
            tuple: subsampled arrays

    """
    if isinstance(random_seed, int):
        np.random.seed(random_seed)
    indices = np.random.choice(
        arrays_to_subsample[0].shape[0], size=size, replace=False
    )
    subsampled_arrays = []
    for array in arrays_to_subsample:
        subsampled_arrays.append(array[indices])
    logging.info(f"Subsampled {size}. Random seed: {random_seed}")
    logging.info([array.shape for array in subsampled_arrays])
    return tuple(subsampled_arrays)


def train_test_split(arrays_to_split: tuple[np.array], fract: float = 0.8) -> tuple:
    """
    Split arrays into train and test sets

        Args:
            arrays_to_split (tuple[np.array]): tuple of arrays to split
            fract (float): fraction of data to use for training
        Returns:
            tuple: train and test sets
    """
    indices = np.random.permutation(arrays_to_split[0].shape[0])
    arrays = []
    for array in arrays_to_split:
        arrays.append(array[indices])
    train_size = int(fract * arrays[0].shape[0])
    train_arrays = []
    test_arrays = []
    for array in arrays:
        train_arrays.append(array[:train_size])
        test_arrays.append(array[train_size:])
    return tuple(train_arrays), tuple(test_arrays)


# def tr_te_split_ids(
#     X: np.array, y: np.array, ids: np.array = np.array([]), fract: float = 0.8
# ) -> tuple:
#     """train test split including ids"""
#     # Split into train and test.
#     if len(ids) == 0:
#         ids = np.arange(X.shape[0])

#     indices = np.random.permutation(X.shape[0])
#     X = X[indices]
#     y = y[indices]
#     ids = ids[indices]

#     train_size = int(fract * X.shape[0])
#     X_train, X_test = X[:train_size], X[train_size:]
#     y_train, y_test = y[:train_size], y[train_size:]
#     ids_train, ids_test = ids[:train_size], ids[train_size:]

#     print(
#         X_train.shape,
#         X_test.shape,
#         y_train.shape,
#         y_test.shape,
#         ids_train.shape,
#         ids_test.shape,
#     )
#     return X_train, X_test, y_train, y_test, ids_train, ids_test


# class rfData:
#     # X = data[0]
#     # y = data[1]
#     # ids = data[2]

#     def __init__(
#         self,
#         X: np.array,
#         cl_path: str,
#         ids: np.array,
#         unilabel: bool,
#         fract: float = 0.8,
#     ):
#         self.X = X
#         self.ids = ids
#         self.unilabel = unilabel
#         self.label_names, self.class_index, self.y = read_classifications(
#             cl_path, unilabel
#         )
#         self.data = (self.X, self.y, self.ids)
#         self.classes = None
#         self.nones_X = None
#         self.nones_y = None
#         self.nones_ids = None
#         self.nones_data = None

#         self.X_train = None
#         self.X_test = None
#         self.y_train = None
#         self.y_test = None
#         self.ids_train = None
#         self.ids_test = None
#         train_test_split(self.data, fract=0.8)

#     def subsample(self, size: int = 10000, random_seed: int = None):
#         if isinstance(random_seed, int):
#             np.random.seed(random_seed)
#         indices = np.random.choice(self.X.shape[0], size=size, replace=False)

#         self.X = self.X[indices]
#         self.y = self.y[indices]
#         self.ids = self.ids[indices]

#         logging.info(f"Subsampled {size}. Random seed: {random_seed}")
#         logging.debug(f"{self.X.shape}, {self.ids.shape}")
#         return None

#     def train_test_split(self, fract: float = 0.8):
#         """regular train test split"""
#         # Split into train and test.
#         indices = np.random.permutation(self.X.shape[0])
#         self.X = self.X[indices]
#         self.y = self.y[indices]
#         self.ids = self.ids[indices]

#         train_size = int(fract * self.X.shape[0])
#         self.X_train, self.X_test = self.X[:train_size], self.X[train_size:]
#         self.y_train, self.y_test = self.y[:train_size], self.y[train_size:]
#         self.ids_train, self.ids_test = self.ids[:train_size], self.ids[train_size:]

#         logging.debug(f"{self.X_train.shape}, {self.X_test.shape}"")
#         return None

#     def separate_nones(self):
#         """separate nones from data"""
#         none_indices = []
#         # get indexes of "None" annotations
#         for i, name in enumerate(self.label_names):
#             if name == ["None"]:
#                 none_indices.append(int(i))

#         # remove "None" class index dictionary (integer values)
#         if "None" in self.class_index.keys():
#             none_ind = self.class_index.pop("None")
#             # if not unilabel:
#             if self.y.shape[1] > 1:
#                 self.y = np.delete(self.y, none_ind, 1)

#         # separate


def train_classifier(
    X_train: np.array,
    y_train: np.array,
    *args,
    **kwargs,
) -> RandomForestClassifier:
    """
    Train random forest classifier

        Args:
            X_train (np.array): training data
            y_train (np.array): training labels
        Returns:
            RandomForestClassifier: trained classifier
    """
    # Train random forest classifier.
    n_jobs = -1  # use all cores
    clf = RandomForestClassifier(n_jobs=n_jobs, *args, **kwargs)
    clf.fit(X_train, y_train)
    return clf


def predict_one(
    X_test: np.array,
    classifier: RandomForestClassifier,
) -> np.array:
    """
    Predict label with RF

        Args:
            X_test (np.array): test data
            classifier (RandomForestClassifier): trained classifier
        Returns:
            np.array: array of predicted labels
    """
    y_pred = classifier.predict(X_test)
    return np.array(y_pred)


def proba(
    X_test: np.array,
    classifier: RandomForestClassifier,
) -> np.array:
    """
    Predict probabilities with RF.

        Args:
            X_test (np.array): test data
            classifier (RandomForestClassifier): trained classifier
        Returns:
            np.array: array of predicted probabilities

    Can cause errors if the sample is too small (<1000) and causes 1.0 probabilities for a class)
    """
    probabilities = np.array(classifier.predict_proba(X_test))
    # only probability of "yes class"
    only_class_membership = probabilities[:, :, 1]
    class_probabilities = np.transpose(only_class_membership)
    return class_probabilities


def cutoffr(y_proba: np.array, cutoff: float = 0.5) -> np.array:
    """
    Round probabilities to 0 or 1 based on cutoff

        Args:
            y_proba (np.array): array of probabilities
            cutoff (float): cutoff for classification. should be between 0.0 and 1.0. (default: 0.5)
        Returns:
            np.array: array of predicted labels
    """
    # make copy of y_proba to not affect the probability array
    y_rounded = np.copy(y_proba)
    # change all values below cutoff to 0 and above/equal cutoff to 1
    y_rounded[y_rounded < cutoff] = 0
    y_rounded[y_rounded >= cutoff] = 1
    return y_rounded.astype(int)


def kfold_yielder(X: np.array, y: np.array, k: int = 5) -> tuple:
    """
    Yield train and test sets for each fold

        Args:
            X (np.array): data
            y (np.array): labels
            k (int): number of folds
        Returns (generator):
            tuple: train and test sets for each fold (X_train, X_test, y_train, y_test)
    """
    k = 5
    kf = KFold(n_splits=k, shuffle=True)
    for train_index, test_index in kf.split(X):
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]
        # logging.info(X_train.shape, X_test.shape, y_train.shape, y_test.shape)
        yield X_train, X_test, y_train, y_test


def kfold_yielder_ids(X: np.array, y: np.array, ids: np.array, k: int = 5) -> tuple:
    """
    Yield train and test sets for each fold, including ids

        Args:
            X (np.array): data
            y (np.array): labels
            ids (np.array): ids
            k (int): number of folds
        Returns (generator):
            tuple: train and test sets for each fold (X_train, X_test, y_train, y_test, ids_train, ids_test)

    """
    k = 5
    kf = KFold(n_splits=k, shuffle=True)
    for train_index, test_index in kf.split(X):
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]
        ids_train, ids_test = ids[train_index], ids[test_index]
        # logging.info(X_train.shape, X_test.shape, y_train.shape, y_test.shape)
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
) -> tuple:
    """
    Perform k-fold cross validation and return confusion matrices and classification reports

        Args:
            X (np.array): data
            y (np.array): labels
            k (int): number of folds
            unilabel (bool): if True, will run single label classification instead of multilabel classification. Default: False
            n_estimators (int): number of trees in the forest. Default: 1000
            max_depth (int): maximum depth of the tree. Default: 100
            cutoff (float): cutoff for classification. Default: 0.5
            target_names (list): list of target names for multilabel classification
        Returns:
            tuple: confusion matrices and classification reports
    """
    # Perform k-fold cross validation.
    cms = []  # confusion matrices
    cl_reps = []  # classification reports
    importances = []
    for X_train, X_test, y_train, y_test in tqdm(kfold_yielder(X, y, k=k), total=k):
        classifier = train_classifier(
            X_train, y_train, n_estimators=n_estimators, max_depth=max_depth
        )
        importances.append(classifier.feature_importances_)
        if unilabel:
            y_pred = predict_one(X_test, classifier)
            cl_rep = classification_report(y_test, y_pred, zero_division=np.nan)
            cm = confusion_matrix(y_test, y_pred)
        else:
            # multilabel, consider cutoff
            y_proba = proba(X_test, classifier)
            y_pred = cutoffr(y_proba, cutoff=cutoff)
            # classification report with classes as target names for multilabel
            cl_rep = classification_report(
                y_test, y_pred, target_names=target_names, zero_division=np.nan
            )
            cm = multilabel_confusion_matrix(y_test, y_pred)
            # feature importance not possible for predict_proba,
            # therefore not included if multilabel (i.e. proba is possible)
        # append results
        cms.append(cm)
        cl_reps.append(cl_rep)
    return cms, cl_reps, importances


def get_wrong_predictions(
    y_test: np.array, y_pred: np.array, ids_test: np.array
) -> tuple:
    """
    Get wrong predictions

        Args:
            y_test (np.array): test labels
            y_pred (np.array): predicted labels
            ids_test (np.array): test ids
        Returns:
            tuple: array of ids, test labels and predicted labels
    """

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
    """
    Get wrong predictions in a numpy array

        Args:
            ids_test (np.array): test ids
            y_test (np.array): test labels
            y_pred (np.array): predicted labels
            class_index (dict): dictionary of classification types
        Returns:
            np.array: array of wrong predictions
    """
    # get wrong predictions.
    w_id, w_test, w_pred = get_wrong_predictions(y_test, y_pred, ids_test)
    # get the classification predictions in strings:
    w_test_str = i_to_cl(w_test, class_index)
    w_pred_str = i_to_cl(w_pred, class_index)
    wrong = np.array([w_id, w_test_str, w_pred_str]).T
    logging.info(f"{wrong.shape[0]} wrong predictions")
    return wrong


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


def write_unilab_confusion_matrices(
    cms: list, classification_types: list, title: str = "confusion_matrices.txt"
) -> None:
    """
    Write unilabel confusion matrices to file

        Args:
            cms (list): list of confusion matrices
            classification_types (list): list of classification types
            title (str): title of file. Default: "confusion_matrices.txt"
        Returns:
            None
    """
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
    """
    Write multilabel confusion matrices to file

        Args:
            cms (list): list of confusion matrices
            classification_types (list): list of classification types
            title (str): title of file. Default: "confusion_matrices.txt"
        Returns:
            None
    """

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
    """
    Combine confusion matrices

        Args:
            cms (list[np.ndarray]): list of confusion matrices
        Returns:
            np.ndarray: combined confusion matrix
    """
    combined = np.zeros(cms[0].shape)
    for cm in cms:
        combined += cm
    return combined.astype(int)


def write_classification_report(
    classification_results: list, classes: list, title="classification_report.txt"
) -> None:
    """
    Write classification report to file

        Args:
            classification_results (list): list of classification reports
            classes (list): list of classification types
            title (str): title of file. Default: "classification_report.txt"
        Returns:
            None
    """
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
    """
    Write arguments as yaml to file

        Args:
            args (argparse.Namespace): arguments
            title (str): title of file. Default: "arguments.yaml"
        Returns:
            None
    """
    args_dict = vars(args)
    # Write arguments as yaml.
    with open(title, "w") as fo:
        for key, value in args_dict.items():
            fo.write(f"{key}: {value}\n")
    return None


def handle_outdirs(db_name: str, fp_name: str, unilabel: bool) -> str:
    """
    Make directories for output

        Args:
            db_name (str): database name
            fp_name (str): fingerprint name
            unilabel (bool): if True, will make unilabel directory
        Returns:
            str: initial working directory
    """
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


def get_decision_trees(
    classifier: RandomForestClassifier, title: str = "decision_trees"
):
    """
    Get decision trees from random forest classifier

        Args:
            classifier (RandomForestClassifier): trained classifier
            title (str): title of file. Default: "decision_trees"
        Returns:
            None

    Should be updated to have updated feature names
    """
    # write decision tree.
    feature_names = [
        "coenzyme_a",
        "nadh",
        "nadph",
        "standard_amino_acids",
        "non-standard_amino_acids",
        "open_pyranose",
        "open_furanose",
        "pyranose",
        "furanose",
        "indoleC2N",
        "phenylC2N",
        "C5N",
        "C4N",
        "phenylC3",
        "phenylC2",
        "phenylC1",
        "isoprene",
        "acetyl",
        "methylmalonyl",
        "ethyl",
        "methyl",
        "phosphate",
        "sulfonate",
        "fluorine",
        "chlorine",
        "bromine",
        "iodine",
        "nitrate",
        "epoxy",
        "ether",
        "hydroxyl",
        "C3_ring",
        "C4_ring",
        "C5_ring",
        "C6_ring",
        "C7_ring",
        "C8_ring",
        "C9_ring",
        "C10_ring",
    ]
    export_graphviz(
        classifier.estimators_[0],
        out_file="tree.dot",
        feature_names=feature_names,
    )
    return None


def main() -> None:
    # set logging to debug
    logging.basicConfig(level=logging.DEBUG)
    n_estimators = 1000  # number of trees
    max_depth = 100
    args = cli()
    logging.info(args)
    include_none = args.include_none  # preserve value of args.include_none for later

    # get fingerprint information from path
    folder = os.path.basename(os.path.dirname(args.fingerprints))
    fp_name = args.fingerprints.split("/")[-1].split(".")[0]
    db_name = folder.split("/")[-1]

    # make directory for output, get initial working directory for later
    iwd = handle_outdirs(db_name, fp_name, args.unilabel)

    # Parse fingerprints from input file.
    delimiter = "\t" if args.fingerprints.endswith(".tsv") else ","
    X = np.loadtxt(args.fingerprints, delimiter=delimiter, dtype=int)
    # try to reduce dimensionality
    # # indices = np.array(
    # #     [3, 4, 5, 11, 14, 17, 18, 21, 22, 24, 27, 28, 29, 30, 31, 32, 33, 34, 36]
    # # )
    # indices = np.array([3, 4, 5, 11, 14, 17, 18, 21, 22, 24, 27, 28, 29, 30, 33, 34])
    # # make mask where indices are true
    # filtr = np.zeros(X.shape[1], dtype=bool)
    # filtr[indices] = True
    # X = X[:, filtr]

    h_labels, cl_idx, y = read_classifications(args.classifications, args.unilabel)
    ids = get_ids(args.names, X.shape[0])
    logging.info(f"X, y:{X.shape}, {y.shape}")

    if not include_none:
        h_labels, cl_idx, filtered, nones = separate_nones(
            h_labels, cl_idx, (X, y, ids)
        )
        X, y, ids = filtered
        nones_X, _, nones_ids = nones
    else:
        nones_X, nones_ids = np.array([]), np.array([])

    logging.debug(cl_idx)
    classes = list(cl_idx.keys())

    # Subsample randomly from X and y.
    if isinstance(args.subsample, int):
        X, y, ids = subsampler(
            (X, y, ids), size=args.subsample, random_seed=args.random_seed
        )

    # Split into train and test.
    train, test = train_test_split((X, y, ids), fract=0.8)
    X_train, y_train, ids_train = train
    X_test, y_test, ids_test = test

    # get initial predictions. ------------------------------------------------------
    logging.info("RF for getting wrong predictions")

    classifier = train_classifier(
        X_train,
        y_train,
        n_estimators=n_estimators,
        max_depth=max_depth,
        # class_weight="balanced",
    )
    importances = classifier.feature_importances_
    if args.no_proba:
        y_pred = predict_one(X_test, classifier)
    else:
        y_proba = proba(X_test, classifier)
        y_pred = cutoffr(y_proba, cutoff=args.cutoff)

    if args.unilabel:
        cm = confusion_matrix(y_test, y_pred)
    else:
        cm = multilabel_confusion_matrix(y_test, y_pred)

    wrongs = wrongs_array(ids_test, y_test, y_pred, cl_idx)
    np.savetxt("wrong_predictions.tsv", wrongs, delimiter="\t", fmt="%s")

    # print(cm, importances, sep="\n\n")
    # print(f'with {X.shape[1]} features and "balanced" class weights')
    # print(classification_report(y_test, y_pred, zero_division=np.nan))

    # none predictions. -------------------------------------------------------------
    if nones_X.shape[0] == 0:
        if not args.include_none:
            logging.info("no nones in dataset")
    else:
        logging.info('get prediction for the "Nones"')
        full_classifier = train_classifier(X, y, n_estimators=n_estimators)
        if args.no_proba:
            ny_pred = predict_one(nones_X, full_classifier)
        else:
            ny_proba = proba(nones_X, full_classifier)
            ny_pred = cutoffr(ny_proba, cutoff=args.cutoff)

        ny_pred = i_to_cl(ny_pred, cl_idx, unilabel=args.unilabel)
        nones_preds = np.concatenate([nones_ids[:, None], ny_pred[:, None]], axis=1)
        np.savetxt("nones_predictions.tsv", nones_preds, delimiter="\t", fmt="%s")
        if not args.no_proba:
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
    # write importances with 3 decimals
    np.savetxt("importances.tsv", importances, delimiter="\t", fmt="%.3f")

    # write arguments as yaml. ------------------------------------------------------
    args_yaml(args, title="arguments.yaml")

    # save model
    if args.export:
        # make directory for model
        if not os.path.exists("./model"):
            os.mkdir("./model")
        os.chdir("./model")

        logging.info("exporting full model")

        # measure the memory usage of training
        tracemalloc.start()
        # time the training precisely
        tic = time.perf_counter()
        full_classifier = train_classifier(X, y, n_estimators=n_estimators)
        toc = time.perf_counter()
        current, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        logging.info(
            f"Current memory usage is {current / 10**6}MB; Peak was {peak / 10**6}MB"
        )
        logging.info(f"training took {toc-tic} seconds")
        np.savetxt(
            "time_train_full.tsv", np.array([toc - tic]), delimiter="\t", fmt="%s"
        )
        np.savetxt(
            "memory_train_full.tsv",
            np.array([peak / 10**6]),
            delimiter="\t",
            fmt="%s",
        )

        pickle.dump(full_classifier, open("model.pkl", "wb"))
        # save labels for indexes
        np.savetxt("model_labels.tsv", classes, delimiter="\t", fmt="%s")
        # export decision tree
        if "full" in args.fingerprints:
            get_decision_trees(full_classifier)
        os.chdir("../")
    # return to initial working directory
    os.chdir(iwd)
    exit(0)
    return None


if __name__ == "__main__":
    main()
