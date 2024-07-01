from enum import Enum

from rdkit import Chem

import numpy as np
import pandas as pd
from tqdm import tqdm

# turn off rdkit warnings
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")


def convert_classes(key: str) -> list:
    cd = {
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
    return ";".join(sorted(set(cd[key])))


def cleanup(query: str) -> str:
    """Clean up query for use as column name."""
    return "_".join(query.lower().split(" "))


def get_properties(sdf: Chem.SDMolSupplier, query: list) -> pd.DataFrame:
    """Get properties from sdf file."""
    all_props = dict()
    for ind, mol in tqdm(enumerate(sdf)):
        if mol:
            all_props[ind] = dict()
            for prop in query:
                try:
                    all_props[ind][prop] = mol.GetProp(prop)
                except:
                    continue
    df = pd.DataFrame.from_dict(all_props, orient="index")
    df["sdf_index"] = df.index
    df.columns = [cleanup(i) for i in df.columns.to_list()]
    return df


def classifications(df, classification_path: str) -> pd.DataFrame:
    """Returns dataframe with indexes of the sdf"""
    classification = pd.read_csv(classification_path, sep="\t", index_col=0)
    class_cols = classification.columns.to_list()
    classification["chebi_id"] = classification.index
    classification = classification.reset_index(drop=True)
    classification.index = classification.index + 50000  # to avoid index overlap

    # get idx, chebi and classification df
    cl_df = pd.merge(df, classification, on="chebi_id", how="inner")
    cl_df = cl_df.set_index("chebi_id")
    # since flavonoids are just subset of (only) phenylpropanoids
    # counts are contained in phenylpropanoids, so we drop flavonoid:
    cl_df = cl_df[class_cols].drop(columns=["flavonoid"])
    cl_df.columns = [convert_classes(i) for i in cl_df.columns.to_list()]

    # merge columns with same name, with OR membership
    return cl_df.groupby(level=0, axis=1).any()


def main():
    sdf_path = "raw_data/ChEBI_complete_3star.sdf"
    classification_path = "raw_data/ChEBI_3star_classifications.tsv"

    # get properties and classifications
    sdf = Chem.SDMolSupplier(sdf_path)
    df = get_properties(sdf, ["ChEBI ID", "SMILES", "InChI"])
    cl_df = classifications(df, classification_path)
    cl_df.to_csv("3star_converted_bool.tsv", sep="\t")  # membership table

    # write classifications to csv
    with open("chebi_classes.csv", "w") as fo:
        for chebi_id, data in cl_df.iterrows():
            true_only = ";".join(cl_df.columns[data.values])  # keeps only True
            true_only = ";".join(set(true_only.split(";")))  # removes duplicates
            fo.write(f"{chebi_id},{true_only}\n")

    # write sdf of classified mols
    _is_classified = df.chebi_id.isin(cl_df.index.to_list())
    sdf_indices = df[_is_classified].sdf_index.to_list()
    writer = Chem.SDWriter("chebi.sdf")
    for index in tqdm(sdf_indices, desc="writing sdf"):
        writer.write(sdf[index])  # keeps properties
    writer.close()

    return None


if __name__ == "__main__":
    main()
