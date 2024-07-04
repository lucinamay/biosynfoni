import sys
from pathlib import Path
from functools import partial

import numpy as np
import pandas as pd
from rdkit import Chem
from tqdm import tqdm

from biosynfoni.inoutput import readr, entry_parser


# ============================== input handling ===============================


def extract(
    lines: list, starts: list = [], strip_extra: str = " - "
) -> dict[list[str]]:
    """
    Extracts lines starting with terms in starts, returning them as a list to accomodate multiple values per entry

        Args:
            lines (list): list of strings
            starts (list): list of strings to search for at the start of each line
            strip_extra (str): string to strip from the start of the line
        Returns:
            dict: dictionary with starts as keys and lists of strings as values

    Remark:
        within a collection of lines, only extracts the lines starting with
        terms in starts, returning them as a list to accomodate multiple values
        per entry
    """
    extraction = {}
    for start in starts:
        extraction[start] = []
    for line in lines:
        for start in starts:
            if line.startswith(f"{start} "):  # include space to avoid partial matches
                extraction[start].append(
                    line.strip().replace(start, "").replace(strip_extra, "").strip()
                )
    return extraction


def info_from_entries(
    info_loc: str,
    info2extract: list[str],
) -> pd.DataFrame:
    """
    Extracts information from entries in a file

        Args:
            info_loc (str): location of file
            info2extract (list): list of strings to search for at the start of each line
        Returns:
            pd.DataFrame: dataframe with info2extract as columns

    Remark:
        within a collection of lines, only extracts the lines starting with
        terms in starts, returning them as a list to accomodate multiple values
        per entry
    """
    entries = entry_parser(
        readr(info_loc, ignore_start="#", encoding="iso8859-1"), sep="//"
    )

    all_vals = [
        extract(entry, starts=info2extract, strip_extra=" - ") for entry in entries
    ]
    return pd.DataFrame.from_records(all_vals)


get_pathways = partial(
    info_from_entries,
    info2extract=[
        "UNIQUE-ID",
        "REACTION-LIST",
        "SPECIES",
        "TAXONOMIC-RANGE",
        "REACTION-LAYOUT",
    ],
)
get_compounds = partial(
    info_from_entries,
    info2extract=["UNIQUE-ID", "SMILES", "INCHI", "NON-STANDARD-INCHI"],
)


def _str_to_mol(structure_string: str) -> Chem.Mol:
    if structure_string.startswith("InChI="):
        mol = Chem.MolFromInchi(structure_string)
    else:
        mol = Chem.MolFromSmiles(structure_string)
    if not mol:
        return pd.NA
    return mol


def write_compounds(df: pd.DataFrame) -> pd.DataFrame:
    """
    Writes csv and sdf files with recalculated INCHI and SMILES and mol objects, respectively

        Args:
            df (pd.DataFrame): dataframe
        Returns:
            pd.DataFrame: dataframe with recalculated INCHI and SMILES
    """
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

    df["mol"] = df.mol.apply(lambda x: _str_to_mol(x) if isinstance(x, str) else pd.NA)
    df = df.dropna(subset=["mol"])
    df["inchi"] = df["mol"].apply(Chem.MolToInchi)
    df["smiles"] = df["mol"].apply(Chem.MolToSmiles)

    writer = Chem.SDWriter("metacyc.sdf")
    for index, row in tqdm(df.iterrows(), total=len(df), desc="writing sdf"):
        mol = row["mol"]
        if mol:
            mol.SetProp("compound_id", index)
            mol.SetProp("smiles", str(row["smiles"]))
            mol.SetProp("inchi", str(row["inchi"]))
            writer.write(mol)
    df[["smiles", "inchi"]].to_csv(
        "metacyc_compounds.tsv", sep="\t", header=True, index=True
    )
    return df


# ============================== formatting df =====================================


def _process_reaction_layout(
    df: pd.DataFrame, column_name: str = "reaction_layout"
) -> pd.DataFrame:
    """
    Splits the reaction-layout column into reaction-id, left, direction, right

        Args:
            df (pd.DataFrame): dataframe with reaction-layout column
            column_name (str): name of the reaction-layout column
        Returns:
            pd.DataFrame: dataframe with reaction-layout column split into reaction-id, left, direction, right
    Remark:
        - specific to pathway-tools dataframe format (metacyc pathway file)
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
    """Removes water from subset columns"""
    for colname in subset:
        df[colname] = df[colname].str.replace(water_str, "")
    return df


def _filter_by_num_reactions(
    df: pd.DataFrame, min_num: int = 4, reaction_list_col="reaction_list"
) -> pd.DataFrame:
    """Filters dataframe by minimum number of reactions in reaction_list_col"""
    df = df[df[reaction_list_col].str.len() >= min_num]
    return df


def _len1list_to_str(df: pd.DataFrame, colname: str) -> pd.DataFrame:
    """Converts list cell to string cell"""
    assert df[df[colname].str.len() > 1].empty, "error, multiple values in cell"
    df[colname] = df[colname].str[0]
    return df


def _rename_col(df: pd.DataFrame, old: str, new: str) -> pd.DataFrame:
    """Renames column in dataframe"""
    df.rename(columns={old: new}, inplace=True)
    return df


def clean_pw_df(
    df: pd.DataFrame, min_num: int = 4, reaction_list_col="reaction_list"
) -> pd.DataFrame:
    """Cleans pathway dataframe"""
    df = _lower_dunder(df)
    df = _process_reaction_layout(df, column_name="reaction_layout")
    df = _switch_directions(df, subset=["left", "right"], on="direction")
    df = _remove_water(df, subset=["left", "right"])
    df = _remove_empties(df, subset=["left", "right"])
    df = _filter_by_num_reactions(
        df, min_num=min_num, reaction_list_col=reaction_list_col
    )
    df = _len1list_to_str(df, "unique_id")
    df = _rename_col(df, "unique_id", "pathway_id")
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


def main():
    raw_data = Path(sys.argv[1]).resolve()
    compounds_path = raw_data / "compounds.dat"  # linked-compounds does not contain all
    pathways_path = raw_data / "pathways.dat"

    meta_mols = get_compounds(compounds_path)
    write_compounds(meta_mols)

    pathways = get_pathways(pathways_path)
    pathways = clean_pw_df(pathways)
    pathways.to_csv("metacyc_pathways.tsv", sep="\t", header=True, index=True)

    return None


if __name__ == "__main__":
    main()
