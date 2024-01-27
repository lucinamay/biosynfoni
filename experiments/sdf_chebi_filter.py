from enum import Enum

from rdkit import Chem

import numpy as np
import pandas as pd
from tqdm import tqdm

# turn off rdkit warnings
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")


def cleanup(query: str) -> str:
    return "_".join(query.lower().split(" "))


def get_properties(sdf: Chem.SDMolSupplier, query: list) -> dict:
    all_props = dict()
    for ind, mol in tqdm(enumerate(sdf)):
        if mol:
            all_props[ind] = dict()
            for prop in query:
                try:
                    all_props[ind][cleanup(prop)] = mol.GetProp(prop)
                except:
                    continue
    return all_props


def main():
    supplier_path = "/Users/lucina-may/thesis/chebiclasses/ChEBI_complete_3star.sdf"
    classification_path = (
        "/Users/lucina-may/thesis/chebiclasses/chebi_classifications.tsv"
    )

    # read sdf
    sdf = Chem.SDMolSupplier(supplier_path)
    # get properties as dictionary
    prop_query = [
        "ChEBI ID",
        "Secondary ChEBI ID",
        "ChEBI Name",
        "SMILES",
        "InChI",
        "InChIKey",
    ]
    all_props = get_properties(sdf, prop_query)

    # make dataframe
    sdf_df = pd.DataFrame.from_dict(all_props, orient="index")
    sdf_df["sdf_index"] = sdf_df.index
    sdf_df.set_index("chebi_id", inplace=True)

    # read classifications
    classification = pd.read_csv(classification_path, sep="\t", index_col=0)
    classes = classification.columns.to_list()

    # join dataframes on inner
    df = pd.concat([sdf_df, classification], axis=1, join="inner")

    # write to file
    df.to_csv("3star_class_props.tsv", sep="\t")
    df[classes].to_csv("3star_classifications.tsv", sep="\t")

    # read back in
    df = pd.read_csv("3star_classifications.tsv", sep="\t", index_col=0)

    # make reduced sdf
    sdf_indexes = df.sdf_index.to_list()
    assert len(sdf_indexes) == len(df[classes])
    writer = Chem.SDWriter("3star_class.sdf")
    for index in tqdm(sdf_indexes, desc="writing sdf"):
        writer.write(sdf[index])
    writer.close()

    # read back in to check
    sdf_reduced = Chem.SDMolSupplier("3star_class.sdf")

    # make classification input ready
    cl_df = df[classes].drop(columns="flavonoid")
    classes = cl_df.columns.to_list()
    conversion = np.loadtxt("conversion.tsv", dtype=str, delimiter="\t")
    conversion = {i[0]: i[1] for i in conversion}
    converted_columns = [conversion[i] for i in classes]
    cl_df.columns = converted_columns
    # merge columns with same name, give true if any of the column values was true
    cl_df = cl_df.groupby(level=0, axis=1).any()
    cl_df.to_csv("3star_converted_bool.tsv", sep="\t")

    # loop over cl_df rows:
    with open("3star_converted_onlytrue.tsv", "w") as fo:
        for row, data in cl_df.iterrows():
            true_only = cl_df.columns[data.values]
            tabbed = "\t".join(true_only)
            fo.write(f"{row}\t{tabbed}\n")

    with open("3star_converted_onlytrue_commad.tsv", "w") as fo:
        for row, data in cl_df.iterrows():
            true_only = cl_df.columns[data.values]
            commad = ",".join(true_only)
            fo.write(f"{commad}\n")

    # run biosynfoni

    return None


if __name__ == "__main__":
    main()
