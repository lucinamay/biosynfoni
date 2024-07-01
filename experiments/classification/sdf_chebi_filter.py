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
                    all_props[ind][cleanup(prop)] = mol.GetProp(prop)
                except:
                    continue
    return pd.DataFrame.from_dict(all_props, orient="index")


def main():
    sdf_path = "raw_data/ChEBI_complete_3star.sdf"
    # classification_path = "chebi_classifications.tsv"
    classification_path = "raw_data/ChEBI_3star_classifications.tsv"

    # get properties as dictionary
    prop_query = [
        "ChEBI ID",
        # "Secondary ChEBI ID",
        # "ChEBI Name",
        "SMILES",
        "InChI",
        # "InChIKey",
    ]

    # make dataframe
    sdf = Chem.SDMolSupplier(sdf_path)
    df = get_properties(sdf, prop_query)
    df["sdf_index"] = df.index

    classification = pd.read_csv(classification_path, sep="\t", index_col=0)
    classification["chebi_id"] = classification.index
    classification.reset_index(inplace=True, drop=True)
    classification.index = classification.index + 50000
    print(classification.head())
    classes = classification.columns.to_list()

    df = pd.merge(df, classification, on="chebi_id", how="inner")
    cl_df = df[classes].drop(columns="flavonoid")
    # since flavonoids are just subset of (only) phenylpropanoids, counts are contained in phenylpropanoids
    cl_df.columns = [convert_classes(i) for i in cl_df.columns.to_list()]
    print(cl_df.head())
    # merge columns with same name, with OR membership
    cl_df = cl_df.groupby(level=0).any()
    cl_df.to_csv("3star_converted_bool.tsv", sep="\t")  # membership table

    # loop over cl_df rows:
    with open("chebi_classes.csv", "w") as fo:
        for row, data in cl_df.iterrows():
            true_only = cl_df.columns[data.values]  # masks out False column values
            tabbed = ";".join(true_only)
            fo.write(f"{row},{tabbed}\n")

    sdf_indices = df.sdf_index.to_list()
    writer = Chem.SDWriter("chebi.sdf")
    for index in tqdm(sdf_indices, desc="writing sdf"):
        writer.write(sdf[index])
    writer.close()

    # sdf_df["sdf_index"] = sdf_df.index
    # sdf_df = sdf_df.set_index("chebi_id")

    # read classifications
    classification = pd.read_csv(classification_path, sep="\t", index_col=0)
    classes = classification.columns.to_list()

    # join mols with their classifications, keep only those that have classifications
    # df = pd.concat(
    #     [sdf_df, classification],
    #     axis=1,
    #     join="inner",
    #     sort=False,
    # )  # inner join on index
    # print(df.head())
    df = pd.merge(df, classification, on="chebi_id", how="inner", left_index=True)
    print(df.head())
    # df.to_csv("3star_class_props.tsv", sep="\t")
    # df[classes].to_csv("3star_classifications.tsv", sep="\t")
    # df = pd.read_csv("3star_classifications.tsv", sep="\t", index_col=0)

    # make sdf of mols with classifications only
    sdf_indices = df.index.to_list()
    assert len(sdf_indices) == len(df[classes])
    writer = Chem.SDWriter("chebi.sdf")
    for index in tqdm(sdf_indices, desc="writing sdf"):

        writer.write(sdf[index])
    writer.close()

    # make classification input ready
    cl_df = df[classes].drop(
        columns="flavonoid"
    )  # since flavonoids are just subset of (only) phenylpropanoids, counts are contained in phenylpropanoids
    cl_df.columns = [convert_classes(i) for i in cl_df.columns.to_list()]
    print(cl_df.head())
    # merge columns with same name, with OR membership
    cl_df = cl_df.groupby(level=0).any()
    cl_df.to_csv("3star_converted_bool.tsv", sep="\t")  # membership table

    # loop over cl_df rows:
    with open("chebi_classes.csv", "w") as fo:
        for row, data in cl_df.iterrows():
            true_only = cl_df.columns[data.values]  # masks out False column values
            tabbed = ";".join(true_only)
            fo.write(f"{row},{tabbed}\n")

    return None


if __name__ == "__main__":
    main()
