import logging, os, sys, subprocess, re
from functools import partial
from pathlib import Path


import numpy as np
import pandas as pd
from tqdm import tqdm
from rdkit import Chem
from rdkit import RDLogger  # for muting warnings

from biosynfoni.inoutput import readr, entry_parser

from helper import ChangeDirectory

RDLogger.DisableLog("rdApp.*")


# ============================== COCONUT =======================================


def consolidate_npclassifier_pathway(key: str) -> list:
    cd = {
        "Polyketides": ["polyketide"],
        "Alkaloids": ["alkaloid"],
        "Amino acids and Peptides": ["amino_acid"],
        "Terpenoids": ["isoprenoid"],
        "Fatty acids": ["fatty_acid"],
        "Carbohydrates": ["carbohydrate"],
        "Shikimates and Phenylpropanoids": ["phenylpropanoid"],
    }
    return ";".join(sorted(set(cd[key]))) if key else ""


def all_natural_products(raw_sdf: Path) -> None:
    """raw_data_folder: path where the complete natural product set, coconut_complete-10-2024.sdf, is located"""
    writer = Chem.SDWriter("coconut.sdf")
    assert raw_sdf.exists(), f"{raw_sdf} not found"
    classification_outfile = "coconut_classes.csv"
    for mol in Chem.SDMolSupplier(raw_sdf):
        if mol:
            writer.write(mol)
            classification = consolidate_npclassifier_pathway(
                mol.GetProp("np_classifier_pathway")
            )
            with open(classification_outfile, "a") as fo:
                fo.write(f"{mol.GetProp("identifier")},{classification}\n")
    writer.close()

    return None


# ============================== SYNTHETICS ====================================


def _zinc_ids(zinc_file) -> list:
    assert os.path.exists(zinc_file), f"{zinc_file} not found"
    cmd = ["grep", "-Eo", "ZINC[0-9]{10,}", f"{zinc_file}"]
    _run = partial(subprocess.run, stdout=subprocess.PIPE, text=True, check=True)
    zinc_id = _run(cmd).stdout.splitlines()
    return [i for i in zinc_id if i]


def synthetic_mols(raw_sdf, np_file, biogenic_file):
    """
    This script is for filtering the synthetic molecules from the ZINC database as provided
    by NaPLeS zenodo files. It reads the synthetic zinc numbers and the zinc numbers and
    smiles from the NaPLeS zenodo files and writes the synthetic smiles to a file.
    """
    assert raw_sdf.exists()
    assert np_file.exists()
    assert biogenic_file.exists()

    writer = Chem.SDWriter("synthetic.sdf")

    with ChangeDirectory(raw_sdf.parent):
        nonsynthetic = list(set(_zinc_ids(np_file) + _zinc_ids(biogenic_file)))
    for mol in Chem.SDMolSupplier(raw_sdf):
        if not mol:
            continue
        if mol.GetProp("zinc_id") in nonsynthetic:
            continue
        mol.SetProp("identifier", mol.GetProp("zinc_id"))
        writer.write(mol)
    writer.close()
    return None


# ============================== CHEBI ==========================================


def consolidate_chebi_classes(key: str) -> list:
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


def get_chebi_properties(sdf: Chem.SDMolSupplier, query: list) -> pd.DataFrame:
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
    df.columns = ["_".join(i.lower().split(" ")) for i in df.columns.to_list()]
    return df


def get_chebi_classifications(df, classification_path: str) -> pd.DataFrame:
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
    cl_df.columns = [consolidate_chebi_classes(i) for i in cl_df.columns.to_list()]

    # merge columns with same name, with OR membership
    return cl_df.groupby(level=0, axis=1).any()


def chebi_and_classification(raw_sdf: Path, raw_classifications: Path):
    assert raw_sdf.exists(), f"{raw_sdf} not found"
    assert raw_classifications.exists(), f"{raw_classifications} not found"

    # get properties and classifications
    sdf = Chem.SDMolSupplier(raw_sdf)
    df = get_chebi_properties(sdf, ["ChEBI ID", "SMILES", "InChI"])
    cl_df = get_chebi_classifications(df, raw_classifications)
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
        mol = sdf[index]
        mol.SetProp("identifier", mol.GetProp("ChEBI ID"))
        writer.write(mol)  # keeps properties
    writer.close()
    return None


# ============================== METACYC ========================================


def extract(lines: list, starts: list = []) -> dict[list[str]]:
    """
    Extracts lines starting with terms in starts, returning them as a list to accomodate multiple values per entry
        Args:
            lines (list): list of strings
            starts (list): list of strings to search for at the start of each line
        Returns:
            dict: dictionary with starts as keys and lists of strings as values
    """
    extraction = {}
    for start in starts:
        extraction[start] = []
    for line in lines:
        for start in starts:
            if line.startswith(f"{start} "):  # include space to avoid partial matches
                extraction[start].append(
                    line.strip().replace(start, "").replace(" - ", "").strip()
                )
    return extraction


def parse_dat(
    info_loc: str,
    info2extract: list[str],
) -> pd.DataFrame:
    """
    Extracts information from entries in a file starting with terms in info2extract
        Args:
            info_loc (str): location of file
            info2extract (list): list of strings to search for at the start of each line
        Returns:
            pd.DataFrame: dataframe with info2extract as columns
    """
    entries = entry_parser(
        readr(info_loc, ignore_start="#", encoding="iso8859-1"), sep="//"
    )

    all_vals = [extract(entry, starts=info2extract) for entry in entries]
    return pd.DataFrame.from_records(all_vals)


def _switch_directions(
    df: pd.DataFrame, subset: tuple[str] = ("left", "right"), on: str = "direction"
) -> pd.DataFrame:
    """
    Switch left and right columns if direction is R2L

        Args:
            df (pd.DataFrame): dataframe with left, right, direction columns
            subset (tuple): tuple of strings for left and right column names
            on (str): name of the direction column
        Returns:
            pd.DataFrame: dataframe with left and right columns switched if direction is R2L

    """
    left, right = subset
    df[left], df[right] = np.where(
        df[on] == "R2L", (df[right], df[left]), (df[left], df[right])
    )
    df[on] = df[on].str.replace("R2L", "L2R")
    return df


def _lower_dunder(df: pd.DataFrame) -> pd.DataFrame:
    """Converts column names to lowercase and replaces - with _"""
    df.columns = df.columns.str.lower().str.replace("-", "_")
    return df


def _remove_empties(
    df: pd.DataFrame, subset: list[str] = ["left", "right"]
) -> pd.DataFrame:
    """Removes empty strings from subset columns"""
    df = df.replace(r"^\s*$", pd.NA, regex=True)
    return df.dropna(subset=subset)


def _remove_water(
    df: pd.DataFrame, subset: list[str] = ["left", "right"], water_str: str = "WATER"
) -> pd.DataFrame:
    """Removes water from subset columns as they are usually side products and
    if they are main products, the reaction is too low-level for natural products anyways
    """
    for colname in subset:
        df[colname] = df[colname].str.replace(water_str, "")
    return df


def _len1list_to_str(df: pd.DataFrame, colname: str) -> pd.DataFrame:
    """Converts list cell to string cell"""
    assert df[df[colname].str.len() > 1].empty, "error, multiple values in cell"
    df[colname] = df[colname].str[0]
    return df


def _parse_reaction(
    df: pd.DataFrame, column_name: str = "reaction_layout"
) -> pd.DataFrame:
    """
    Splits the reaction-layout column into reaction-id, left, direction, right
    """
    # first get separate entry for each REACTION-LAYOUT
    df = df.explode(column_name)
    df["reaction_id"] = df[column_name].str.split(" ", expand=True)[0].str.strip("(")
    df["left"] = (
        df[column_name]
        .str.split(":", expand=True)[1]
        .str.replace("LEFT-PRIMARIES", "")
        .str.strip(") (")
    )
    df["direction"] = df[column_name].str.split(":", expand=True)[3].str.strip(") (")
    df["right"] = (
        df[column_name]
        .str.split("RIGHT-PRIMARIES", expand=True)[1]
        .str.strip("))")
        .str.strip()
    )
    return df.drop(columns=[column_name])


def read_pathways(pathways_path: Path) -> pd.DataFrame:
    df = parse_dat(
        pathways_path,
        info2extract=[
            "UNIQUE-ID",
            "REACTION-LIST",
            "SPECIES",
            "TAXONOMIC-RANGE",
            "REACTION-LAYOUT",
        ],
    )

    df = _lower_dunder(df)
    df = _parse_reaction(df, column_name="reaction_layout")
    df = _switch_directions(df, subset=["left", "right"], on="direction")
    df = _remove_water(df, subset=["left", "right"])
    df = _remove_empties(df, subset=["left", "right"])
    df = _len1list_to_str(df, "unique_id")
    df.rename(columns={"unique_id": "pathway_id"}, inplace=True)
    return df[
        [
            "pathway_id",
            "reaction_id",
            "left",
            "direction",
            "right",
            "reaction_list",
            "species",
            "taxonomic_range",
        ]
    ]


def _str_to_mol(structure_string: str) -> Chem.Mol:
    if not isinstance(structure_string, str):
        return pd.NA
    if structure_string.startswith("InChI="):
        mol = Chem.MolFromInchi(structure_string)
    else:
        mol = Chem.MolFromSmiles(structure_string)
    if not mol:
        return pd.NA
    return mol


def get_compounds(compounds_path: Path) -> pd.DataFrame:
    df = parse_dat(
        compounds_path,
        info2extract=[
            "UNIQUE-ID",
            "SMILES",
            "INCHI",
            "NON-STANDARD-INCHI",
        ],
    )
    for col in df.columns:
        df = _len1list_to_str(df, col)
    df = _lower_dunder(df)
    df = df.set_index("unique_id")
    df.index.name = "compound_id"
    if "non_standard_inchi" not in df.columns:
        df["non_standard_inchi"] = pd.NA

    # inchi has less chance of having invalid parts such as [a arene] or sth
    df["mol"] = (
        df["inchi"].fillna(df["smiles"]).fillna(df["non_standard_inchi"]).fillna(pd.NA)
    )

    df["mol"] = df.mol.apply(_str_to_mol)
    df = df.dropna(subset=["mol"])
    df["inchi"] = df["mol"].apply(Chem.MolToInchi)
    df["smiles"] = df["mol"].apply(Chem.MolToSmiles)
    return df


def write_compound_tsv_sdf(df: pd.DataFrame) -> None:
    """
    Writes tsv and sdf files with recalculated INCHI and SMILES and mol objects, respectively
        Args:
            df (pd.DataFrame): dataframe
        Returns:
            None, writes sdf and tsv files
    """

    writer = Chem.SDWriter("metacyc.sdf")
    for index, row in tqdm(df.iterrows(), total=len(df), desc="writing sdf"):
        mol = row["mol"]
        if mol:
            mol.SetProp(
                "compound_id", index
            )  # @TODO: add in the actual MetaCyc ID as "identifier"
            mol.SetProp("smiles", str(row["smiles"]))
            mol.SetProp("inchi", str(row["inchi"]))
            writer.write(mol)
    df[["smiles", "inchi"]].to_csv(
        "metacyc_compounds.tsv", sep="\t", header=True, index=True
    )


def merge_merges(on_inchi: pd.DataFrame, on_smiles: pd.DataFrame) -> pd.DataFrame:
    """
    Merges the two dataframes resulting from merge_on_partial_inchi and merge_on_smiles
    """
    print(on_inchi.columns, on_smiles.columns, "\n\n")
    on_inchi.rename(columns={"smiles_meta": "smiles"}, inplace=True)
    # on_smiles.rename(columns={"inchi_meta": "inchi"}, inplace=True)
    print(on_inchi.columns, on_smiles.columns)
    print(on_inchi.head(), on_smiles.head())
    df = pd.merge(
        on_inchi,
        on_smiles,
        on=["identifier", "compound_id", "inchi_meta", "smiles"],
        how="outer",
        suffixes=("_inchi", "_smiles"),
    )
    df.rename(columns={"inchi_meta": "inchi"}, inplace=True)
    return df


def get_common_compounds(
    coco_mols: pd.DataFrame, meta_mols: pd.DataFrame
) -> pd.DataFrame:
    """
    Gets common compounds from coco and meta mols on partial inchi and smiles
        Args:
            coco_mols (pd.DataFrame): dataframe
            meta_mols (pd.DataFrame): dataframe
        Returns:
            pd.DataFrame: dataframe with common compounds' coconut and meta ids, inchi and smiles
    """
    # get the index as a column
    meta_mols = meta_mols.reset_index()
    coco_mols = coco_mols.reset_index()

    coco_mols = coco_mols.rename(
        columns={"canonical_smiles": "smiles", "standard_inchi": "inchi"}
    )

    _partial_inchi = lambda x: (
        "/".join(x.split("/")[:3]) if isinstance(x, str) else np.nan
    )
    meta_mols["partial_inchi"] = meta_mols["inchi"].apply(_partial_inchi)
    coco_mols["partial_inchi"] = coco_mols["inchi"].apply(_partial_inchi)

    logging.info(
        (
            f"{len(coco_mols[~coco_mols["partial_inchi"].isna()])} coco mols have inchi",
            f"{len(meta_mols[~meta_mols["partial_inchi"].isna()])} meta mols have inchi",
        )
    )  # 406529

    same_p_inchi = pd.merge(
        coco_mols,
        meta_mols,
        on="partial_inchi",
        how="inner",
        suffixes=("_coco", "_meta"),
    )
    same_smiles = pd.merge(
        coco_mols, meta_mols, on="smiles", how="inner", suffixes=("_coco", "_meta")
    )
    coco_meta = merge_merges(same_p_inchi, same_smiles)

    logging.info(f"{len(same_p_inchi)} compounds have a partial inchi in common")
    logging.info(f"{len(same_smiles)} compounds have a smiles in common")
    logging.info(f"{len(coco_meta)} compounds in common")  # 406529

    return coco_meta[["identifier", "compound_id", "inchi", "smiles"]]


# ------------------------ get relevant pathways -----------------------------


def get_pathways_of_interest(
    pathways: pd.DataFrame, containing: list[str]
) -> pd.DataFrame:
    """
    Returns pathways containing compounds of interest

        Args:
            pathways (pd.DataFrame): dataframe
            containing (list of str): list of compounds of interest
        Returns:
            pd.DataFrame: dataframe of pathways containing compounds of interest

    Remark:
        - filters pathways on common compounds by "compound_id"

    """
    # group by pathway
    pws = []
    for pathway_id, group in pathways.groupby("pathway_id"):
        for compound in containing:
            if compound in " ".join(group["left"].tolist()).split():
                pws.append(pathway_id)
                break
            if compound in " ".join(group["right"].tolist()).split():
                pws.append(pathway_id)
                break
    return pathways[pathways["pathway_id"].isin(pws)]


def metacyc_and_pathways(
    compounds_path: Path, pathways_path: Path, coconut_csv: str
) -> None:

    temp_dir = Path("intermediate_files")
    temp_dir.mkdir(parents=True, exist_ok=True)

    meta_mols = get_compounds(compounds_path)
    write_compound_tsv_sdf(meta_mols)

    common_compound_path = temp_dir / "metacyc_np_compounds.tsv"
    if not common_compound_path.exists():
        coco_mols = pd.read_csv(coconut_csv, index_col="identifier")
        compounds = get_common_compounds(coco_mols, meta_mols)
        compounds.to_csv(common_compound_path, sep="\t", index=False)
    compounds = pd.read_csv(common_compound_path, sep="\t")

    pathways = read_pathways(pathways_path)
    pws_oi = get_pathways_of_interest(pathways, compounds["compound_id"].unique())
    pws_oi.to_csv("metacyc_pathways.tsv", sep="\t", header=True, index=True)


def main():
    raw_data_folder = Path(sys.argv[1]).resolve(strict=True)

    raw_coconut_sdf = raw_data_folder / "coconut_complete-10-2024.sdf"

    raw_zinc_sdf = raw_data_folder / "zinc-all-for-sale.sdf"
    raw_zinc_natural_products = raw_data_folder / "13062018.natural-products.sdf"
    raw_zinc_biogenic = raw_data_folder / "12.07.2018.biogenic.smi"

    raw_chebi_sdf = raw_data_folder / "ChEBI_complete_3star.sdf"
    raw_chebi_classifications = raw_data_folder / "ChEBI_3star_classifications.tsv"

    raw_metacyc_compounds = raw_data_folder / "compounds.dat"
    raw_metacyc_pathways = raw_data_folder / "pathways.dat"
    raw_coconut_csv = raw_data_folder / "coconut_complete-10-2024.csv"

    with ChangeDirectory(raw_data_folder.parent / "input"):
        # all_natural_products(raw_coconut_sdf)
        # synthetic_mols(raw_zinc_sdf, raw_zinc_natural_products, raw_zinc_biogenic)
        # chebi_and_classification(raw_chebi_sdf, raw_chebi_classifications)
        metacyc_and_pathways(
            raw_metacyc_compounds, raw_metacyc_pathways, raw_coconut_csv
        )

    return None


if __name__ == "__main__":
    main()
