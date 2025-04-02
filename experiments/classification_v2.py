# -*- coding: utf-8 -*-

"""This script tests the biosyfoni fingerprints for classifying compounds based on biosynthetic class."""

import argparse
import gzip
import logging
import os
import shutil
import typing as ty
from collections import Counter

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
from biosynfoni import Biosynfoni
from rdkit import Chem, RDLogger
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator
from sklearn.cluster import KMeans
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV, GroupKFold
from sklearn.metrics import f1_score, make_scorer, multilabel_confusion_matrix
from tqdm import tqdm


# Set random seed for reproducibility.
np.random.seed(42)


# Allow usage progress_apply for pandas DataFrame.
tqdm.pandas()


CHEBI_COMPLETE_3STAR_LINK = "https://ftp.ebi.ac.uk/pub/databases/chebi/SDF/ChEBI_complete_3star.sdf.gz"
CHEBI_3STAR_CLASSIFICATIONS_LINK = "https://zenodo.org/records/14791205/files/ChEBI_3star_classifications.tsv?download=1"


CHEBI_TARGET_CLASSES = [
    "alkaloid",
    "amino_acid",
    "carbohydrate",
    "fatty_acid",
    "isoprenoid",
    "phenylpropanoid",
    "polyketide",
]
CHEBI_CLASS_MAPPING = {
        "alkaloid": ["alkaloid"],
        "aminoacid": ["amino_acid"],
        "aminofattyacid": ["amino_acid", "fatty_acid"],
        "aminoglycan": ["amino_acid", "carbohydrate"],
        "aminoglycoside": ["amino_acid", "carbohydrate"],
        "carbohydrate": ["carbohydrate"],
        "fattyacid": ["fatty_acid"],
        "flavanglycoside": ["carbohydrate", "phenylpropanoid"],
        "flavones": ["phenylpropanoid"],
        "flavonoid": ["phenylpropanoid"],
        "flavonoids": ["phenylpropanoid"],
        "glycoalkaloid": ["alkaloid", "carbohydrate"],
        "glycolipid": ["carbohydrate", "fatty_acid"],
        "glycopeptide": ["amino_acid", "carbohydrate"],
        "glycosaminoglycan": ["amino_acid", "carbohydrate"],
        "glycosinolate": ["amino_acid", "carbohydrate"],
        "glycosyloxyflavone": ["carbohydrate", "phenylpropanoid"],
        "glycosyloxyisoflavone": ["carbohydrate", "phenylpropanoid"],
        "isoprenoid": ["isoprenoid"],
        "lipid": ["fatty_acid"],
        "lipopeptide": ["amino_acid", "fatty_acid"],
        "liposaccharide": ["carbohydrate", "fatty_acid"],
        "peptide": ["amino_acid"],
        "phenylpropanoid": ["phenylpropanoid"],
        "polyketide": ["polyketide"],
        "polysaccharide": ["carbohydrate"],
        "saccharolipid": ["carbohydrate", "fatty_acid"],
        "steroid": ["isoprenoid"],
        "steroidalkaloid": ["alkaloid", "isoprenoid"],
        "terpene": ["isoprenoid"],
        "terpenealkaloid": ["alkaloid", "isoprenoid"],
        "terpeneglycoside": ["carbohydrate", "isoprenoid"],
        "terpenoid": ["isoprenoid"],
}


def setup_logging(logger_level: str = "INFO", logger_file_path: ty.Optional[str] = None) -> None:
    """Set up logging configuration.

    :param logger_level: logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL), default is INFO.
    :param logger_file_path: path to the log file, default is None.
    """
    # Delete log file if it exists.
    if logger_file_path and os.path.exists(logger_file_path):
        os.remove(logger_file_path)

    # Create logger instance.
    logger = logging.getLogger()
    logger.setLevel(logger_level)

    # Specify logging format.
    formatter = logging.Formatter("[%(asctime)s - %(name)s - %(levelname)s] %(message)s")

    # Create console handler and set level.
    ch = logging.StreamHandler()
    ch.setLevel(logger_level)
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    # Create file handler and set level.
    fh = logging.FileHandler(logger_file_path)
    fh.setLevel(logger_level)
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    return logger


def cli() -> argparse.Namespace:
    """Command line interface.
    
    :return: Namespace: parsed arguments.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--out-dir", 
        type=str, 
        required=True, 
        help="Output directory to be used (and created if it doesn't exist)."
    )
    parser.add_argument(
        "--logger-level",
        type=str,
        default="INFO",
        help="Logger level (default: DEBUG).",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]
    )
    return parser.parse_args()


def download_file(url: str, out_path: str, logger: logging.Logger) -> None:
    """Download a file from a URL and save it to a specified path.

    :param url: URL of the file to download.
    :param out_path: Path to save the downloaded file.
    :param logger: Logger instance.
    """
    # Download the file.
    logger.info(f"Downloading {url} to {out_path}...")
    with requests.get(url, stream=True) as response:
        response.raise_for_status()
        total_size = int(response.headers.get("content-length", 0))
        block_size = 1024  # 1 KB
        with open(out_path, "wb") as file, tqdm(
            desc=out_path,
            total=total_size,
            unit="iB",
            unit_scale=True,
            unit_divisor=1024,
        ) as bar:
            for data in response.iter_content(block_size):
                file.write(data)
                bar.update(len(data))
    logger.info(f"Downloaded {url} to {out_path}.")


def download_data(out_dir: str, logger: logging.Logger) -> ty.Tuple[str, str]:
    """Download required data files.

    :param out_dir: Output directory to save the downloaded files.
    :param logger: Logger instance.
    :return: Tuple of paths to the downloaded files.
    """
    chebi_complete_3star_path = os.path.join(out_dir, "ChEBI_complete_3star.sdf.gz")
    unzipped_file_path_chebi_complete_3star = chebi_complete_3star_path.replace(".gz", "")
    chebi_3star_classifications_path = os.path.join(out_dir, "ChEBI_3star_classifications.tsv")

    # Check if the files already exist before downloading. Also skip if the unzipped file already exists.
    if not os.path.exists(chebi_complete_3star_path) and not os.path.exists(unzipped_file_path_chebi_complete_3star):
        download_file(CHEBI_COMPLETE_3STAR_LINK, chebi_complete_3star_path, logger) 
    else:
        logger.info(f"{chebi_complete_3star_path} already exists. Skipping download.")

    # Check if the classifications file already exists before downloading.
    if not os.path.exists(chebi_3star_classifications_path):
        download_file(CHEBI_3STAR_CLASSIFICATIONS_LINK, chebi_3star_classifications_path, logger)
    else:
        logger.info(f"{chebi_3star_classifications_path} already exists. Skipping download.")

    # If the unzipped sdf file not exists, unzip the gzipped file.
    if not os.path.exists(unzipped_file_path_chebi_complete_3star):
        with gzip.open(chebi_complete_3star_path, "rb") as f_in:
            with open(unzipped_file_path_chebi_complete_3star, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        logger.info(f"Unzipped {chebi_complete_3star_path} to {unzipped_file_path_chebi_complete_3star}.")
    else:
        logger.info(f"{unzipped_file_path_chebi_complete_3star} already exists. Skipping unzipping.")

    return unzipped_file_path_chebi_complete_3star, chebi_3star_classifications_path


def retrieve_annotated_data(out_dir: str, logger: logging.Logger) -> pd.DataFrame:
    """Retrieve annotated data from ChEBI SDF file and classifications.

    :param out_dir: Output directory to save the downloaded files.
    :param logger: Logger instance.
    :return: DataFrame containing ChEBI ID, SMILES, and target classes.
    """
    # Download required data if not already present in the output directory.
    chebi_sdf_path, chebi_classes_path = download_data(out_dir, logger)

    # Create ChEBI->classification mapping.
    chebi_classes_df = pd.read_csv(chebi_classes_path, sep="\t", header=0, index_col=0)  # Index is ChEBI ID, every column is a class with boolean values.
    chebi_classes_lookup = chebi_classes_df.to_dict(orient="index")
    
    # Read the SDF file and create a DataFrame with ChEBI ID and SMILES. Only retain ChEBI IDs that are in the classification file.
    succeeded_to_read, failed_to_read, wildcard_in_smiles = 0, 0, 0
    sdf_supplier = Chem.SDMolSupplier(chebi_sdf_path, removeHs=False)
    chebi_smiles_lookup = {}
    for mol in tqdm(sdf_supplier, desc="Processing ChEBI SDF file", unit="mol"):
        if mol is None:
            failed_to_read += 1
            continue
        chebi_id = mol.GetProp("ChEBI ID")
        if chebi_id in chebi_classes_lookup:
            smiles = Chem.MolToSmiles(mol)

            # Check for wildcard characters in SMILES.
            if "*" in smiles:
                wildcard_in_smiles += 1
                continue

            chebi_smiles_lookup[chebi_id] = smiles
            succeeded_to_read += 1
    logger.info(f"Successfully read {succeeded_to_read} compounds from the ChEBI SDF file.")
    logger.info(f"Failed to read {failed_to_read} compounds from the ChEBI SDF file.")
    logger.info(f"Found {wildcard_in_smiles} compounds with wildcard characters in SMILES.")

    # Create a DataFrame with ChEBI ID, SMILES, and ChEBI target classes.
    chebi_df = pd.DataFrame.from_dict(chebi_smiles_lookup, orient="index", columns=["SMILES"])
    chebi_df["ChEBI ID"] = chebi_df.index
    chebi_df = chebi_df.reset_index(drop=True)
    
    # Get all classes for ChEBI ID and convert them to target classes.
    # Check for which original classes ChEBI ID has True, and retrieve target classes for each True original class.
    # Then convert target classes to set to remove duplicates, and convert to list. Every target class will be a column in the DataFrame.
    chebi_df["ChEBI classes"] = chebi_df["ChEBI ID"].apply(lambda x: [k for k, v in chebi_classes_lookup[x].items() if v])
    chebi_df["ChEBI target classes"] = chebi_df["ChEBI classes"].apply(
        lambda x: list(set([target_class for original_class in x for target_class in CHEBI_CLASS_MAPPING.get(original_class, [])]))
    )
    # Convert target classes to boolean values. Every target class will be a column in the DataFrame.
    for target_class in CHEBI_TARGET_CLASSES:
        chebi_df[target_class] = chebi_df["ChEBI target classes"].apply(lambda x: target_class in x)

    # Remove original classes from the DataFrame.
    chebi_df = chebi_df.drop(columns=["ChEBI classes"]) 

    # Remove original target classes from the DataFrame.
    chebi_df = chebi_df.drop(columns=["ChEBI target classes"])

    return chebi_df


def smiles_to_biosynfoni_fingerprint(smiles: str) -> ty.List[int]:
    """Convert SMILES to a biosynfoni fingerprint.

    :param smiles: SMILES string to convert.
    :return: List of integers representing the biosynfoni fingerprint.
    """
    mol = Chem.MolFromSmiles(smiles)
    fingerprint = Biosynfoni(mol).fingerprint
    return fingerprint


def smiles_morgan_fingerprint(smiles: str, radius: int, num_bits: int) -> np.ndarray:
    """Convert SMILES to a Morgan fingerprint.
    
    :param smiles: SMILES string to convert.
    :param radius: Radius for the Morgan fingerprint.
    :param num_bits: Number of bits for the Morgan fingerprint.
    :return: Numpy array representing the Morgan fingerprint.
    """
    mol = Chem.MolFromSmiles(smiles)
    fingerprint_generator = GetMorganGenerator(radius=radius, fpSize=num_bits)
    return np.array(fingerprint_generator.GetFingerprint(mol))


def evaluate(
    split_name: str,
    analysis_name: str, 
    X_train: np.ndarray, 
    X_test: np.ndarray,
    y_train: np.ndarray, 
    y_test: np.ndarray,
    cluster_labels_train: np.ndarray,
    out_dir: str, 
    logger: logging.Logger,
) -> None:
    """Evaluate the model using RandomForestClassifier and GridSearchCV.

    :param split_name: Name of the split.
    :param analysis_name: Name of the analysis.
    :param X_train: Training data features.
    :param X_test: Test data features.
    :param y_train: Training data labels.
    :param y_test: Test data labels.
    :param cluster_labels_train: Cluster labels for the training data.
    :param out_dir: Output directory to save the results.
    :param logger: Logger instance.
    """
    logger.info(f"Evaluating {split_name}-{analysis_name}...")

    # (Hyper)parameter grid for RandomForest Classifier
    param_grid = {
        "n_estimators": [100, 500, 1000],
        "max_depth": [None, 10, 50, 100],
        "max_features": ["sqrt", "log2"],
    }

    # Calculate micro averaged F1 scores for each parameter combination, report on average score and its std over the folds.
    rf = RandomForestClassifier(random_state=42, n_jobs=-1)
    group_kfold = GroupKFold(n_splits=5)
    scorer = make_scorer(f1_score, average="macro")
    grid_search = GridSearchCV(rf, param_grid, cv=group_kfold, scoring=scorer, n_jobs=-1, verbose=3)
    grid_search.fit(X_train, y_train, groups=cluster_labels_train)  # Make sure cluster is either in train or test set.

    # Get scores for all parameter combinations, and its standard deviations over the folds.
    results = grid_search.cv_results_
    labeled_values = []
    for i in range(len(results["params"])):
        labeled_values.append({
            "params": results["params"][i],
            "mean_test_score": results["mean_test_score"][i],
            "std_test_score": results["std_test_score"][i],
            "rank_test_score": results["rank_test_score"][i],
        })

    # Create a DataFrame with the labeled values.
    grid_search_results_df = pd.DataFrame(labeled_values)

    # For every parameter combination, train the model and report on the test score.
    test_scores = []
    for i in tqdm(range(len(grid_search_results_df)), desc="Evaluating test scores", unit="param combination"):
        params = grid_search_results_df.iloc[i]["params"]
        rf = RandomForestClassifier(**params, random_state=42, n_jobs=-1)
        rf.fit(X_train, y_train)
        y_pred = rf.predict(X_test)
        test_score = f1_score(y_test, y_pred, average="macro")
        test_scores.append(test_score)

    # Add test scores to the DataFrame.
    grid_search_results_df["left_out_score"] = test_scores

    # Sort first on rank_test_score, then on test_score.
    grid_search_results_df = grid_search_results_df.sort_values(by=["rank_test_score", "left_out_score"], ascending=[True, False])

    # Save DataFrame to TSV file.
    grid_search_results_path = os.path.join(out_dir, f"grid_search_results_split_{split_name}_analysis_{analysis_name}.tsv")
    grid_search_results_df.to_csv(grid_search_results_path, sep="\t", index=False)
    logger.info(f"Saved grid search results for analysis {split_name}-{analysis_name} to {grid_search_results_path}.")

    # Also get best set of parameters, predict on the test set and write confusion matrices out to files.
    best_params = grid_search.best_params_
    logger.info(f"Best parameters for {split_name}-{analysis_name}: {best_params}")
    rf = RandomForestClassifier(**best_params, random_state=42, n_jobs=1)
    rf.fit(X_train, y_train)
    y_pred = rf.predict(X_test)
    cm = multilabel_confusion_matrix(y_test, y_pred)
    
    # Create single figure with all confusion matrices.
    fig, axes = plt.subplots(nrows=len(CHEBI_TARGET_CLASSES), ncols=1, figsize=(10, 5 * len(CHEBI_TARGET_CLASSES)))
    fig.suptitle(f"Confusion matrices for {split_name}-{analysis_name} with best parameters: {best_params}", fontsize=16)
    for i, target_class in enumerate(CHEBI_TARGET_CLASSES):
        axes[i].imshow(cm[i], interpolation="nearest", cmap=plt.cm.Blues)
        axes[i].set_title(f"Confusion matrix for {target_class}")
        axes[i].set_xlabel("Predicted label")
        axes[i].set_ylabel("True label")
        axes[i].set_xticks([0, 1])
        axes[i].set_yticks([0, 1])
        axes[i].set_xticklabels(["Not " + target_class, target_class])
        axes[i].set_yticklabels(["Not " + target_class, target_class])
        axes[i].text(0, 0, cm[i][0][0], ha="center", va="center", color="black")
        axes[i].text(0, 1, cm[i][0][1], ha="center", va="center", color="black")
        axes[i].text(1, 0, cm[i][1][0], ha="center", va="center", color="black")
        axes[i].text(1, 1, cm[i][1][1], ha="center", va="center", color="black")
    plt.tight_layout()
    plt.subplots_adjust(top=0.9)
    plt.savefig(os.path.join(out_dir, f"confusion_matrices_split_{split_name}_analysis_{analysis_name}.png"))
    plt.close(fig)
    logger.info(f"Saved confusion matrices for {split_name}-{analysis_name} to {out_dir}.")
    logger.info(f"Finished evaluating {split_name}-{analysis_name}.")
    

def main() -> None:
    """Entry point for the script."""
    # Parse command line arguments.
    args = cli()

    # Mute RDKit warnings if logger level is not DEBUG.
    if args.logger_level != "DEBUG":
        RDLogger.DisableLog("rdApp.*")

    # Create output directory if it doesn't exist.
    os.makedirs(args.out_dir, exist_ok=True)

    # Set up logging.
    logger = setup_logging(
        logger_level=args.logger_level, 
        logger_file_path=os.path.join(args.out_dir, f"{args.logger_level.lower()}_log.txt")
    )
    logger.info(f"command line arguments: {args}")

    # Define output path for parsed compound data.
    data_path = os.path.join(args.out_dir, "annotated_compounds.tsv")

    if not os.path.exists(data_path):
        chebi_df = retrieve_annotated_data(args.out_dir, logger)

        # Save the DataFrame to a CSV file.
        chebi_df.to_csv(data_path, index=False, header=True, sep="\t")
        logger.info(f"Saved annotated compounds to {data_path}.")
    else:
        # If the file already exists, read it.
        chebi_df = pd.read_csv(data_path, sep="\t", header=0)
        logger.info(f"Loaded annotated compounds from {data_path}.")

    # Report on number of compounds and number of classes.
    logger.info(f"Number of compounds: {len(chebi_df)}")
    logger.info(f"Number of classes: {len(chebi_df.columns) - 2}")  # Exclude ChEBI ID and SMILES columns.
    logger.info(f"Classes: {list(chebi_df.columns[2:])}")

    # Report on class balance for every target class.
    class_counts = Counter()
    for target_class in CHEBI_TARGET_CLASSES:
        class_counts[target_class] = chebi_df[target_class].sum()
    logger.info("Class balance:")
    for target_class, count in class_counts.items():
        logger.info(f"  * {target_class}: {count} / {len(chebi_df)} ({count / len(chebi_df) * 100:.2f}%)")

    # Featurize data in several ways and evaluate.
    split_names = ["biosynfoni", "morgan39", "morgan2048"]
    analysis_names = ["biosynfoni", "morgan39", "morgan2048"]

    # Featurize data in several ways.
    X_biosyfoni = np.array(chebi_df["SMILES"].progress_apply(smiles_to_biosynfoni_fingerprint).tolist())
    X_morgan39 = np.array(chebi_df["SMILES"].progress_apply(lambda x: smiles_morgan_fingerprint(x, radius=2, num_bits=39)).tolist())
    X_morgan2048 = np.array(chebi_df["SMILES"].progress_apply(lambda x: smiles_morgan_fingerprint(x, radius=2, num_bits=2048)).tolist())
    y = np.array(chebi_df[CHEBI_TARGET_CLASSES].values.tolist())

    for split_name in split_names:

        # Select X for splitting.
        if split_name == "biosynfoni":
            X_splitting = X_biosyfoni
        elif split_name == "morgan39":
            X_splitting = X_morgan39
        elif split_name == "morgan2048":
            X_splitting = X_morgan2048
        else:
            raise ValueError(f"Unknown split name: {split_name}")

        # Cluster data and select clusters for the left-out set.
        num_clusters_target = 1000
        cluster_labels = KMeans(n_clusters=num_clusters_target, random_state=42).fit_predict(X_splitting)
        # Randomly select ratio of the clusters for the test set.
        test_ratio = 0.3
        min_samples_per_class = 100
        while True:
            test_clusters = np.random.choice(np.unique(cluster_labels), size=int(num_clusters_target * test_ratio), replace=False)
            inds_test = np.where(np.isin(cluster_labels, test_clusters))[0]
            inds_train = np.where(~np.isin(cluster_labels, test_clusters))[0]
            y_test = y[inds_test]
            X_splitting_train = X_splitting[inds_train]  # Need those for clustering.
            y_train = y[inds_train]
            # Check if every target class is represented in the test and train sets with at least N samples.
            if (
                all(np.sum(y_test[:, i]) >= min_samples_per_class for i in range(y_test.shape[1]))
                and all(np.sum(y_train[:, i]) >= min_samples_per_class for i in range(y_train.shape[1]))
            ):
                break
            else:
                logger.warning(f"Test and/or train clusters for analysis {analysis_name} not valid, retrying...")

        # Cluster the train data by KMeans.
        num_clusters_target = 1000
        cluster_labels_train = KMeans(n_clusters=num_clusters_target, random_state=42).fit_predict(X_splitting_train)

        # Report on sizes of train and test sets.
        logger.info(f"Train set size {split_name}: {len(inds_train)}")
        logger.info(f"Test set size {split_name}: {len(inds_test)}")

        # Split the data into training and test sets.
        for analysis_name in analysis_names:
            
            # Based on the splitting indices, select for the correct data for training and testing.
            if analysis_name == "biosynfoni":
                X_train = X_biosyfoni[inds_train]
                X_test = X_biosyfoni[inds_test]
            elif analysis_name == "morgan39":
                X_train = X_morgan39[inds_train]
                X_test = X_morgan39[inds_test]
            elif analysis_name == "morgan2048":
                X_train = X_morgan2048[inds_train]
                X_test = X_morgan2048[inds_test]
            else:
                raise ValueError(f"Unknown analysis name: {analysis_name}")
            
            # Report on experiment:
            logger.info(f"Running experiment {split_name}-{analysis_name}...")
            logger.info(f"X_train shape: {X_train.shape}")
            logger.info(f"X_test shape: {X_test.shape}")
            logger.info(f"y_train shape: {y_train.shape}")
            logger.info(f"y_test shape: {y_test.shape}")
            logger.info(f"cluster_labels_train shape: {cluster_labels_train.shape}")

            # Evaluate the model.
            evaluate(
                split_name=split_name,
                analysis_name=analysis_name,
                X_train=X_train,
                X_test=X_test,
                y_train=y_train,
                y_test=y_test,
                cluster_labels_train=cluster_labels_train,
                out_dir=args.out_dir,
                logger=logger,
            )


if __name__ == "__main__":
    main()
