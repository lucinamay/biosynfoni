#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 16:15:34 2023

@author: lucina-may
"""
import sys, subprocess, os, logging, argparse
from sys import argv
from functools import partial


import numpy as np
import pandas as pd
from tqdm import tqdm
from rdkit import Chem

from biosynfoni.inoutput import *


def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-p",
        "--pathways_path",
        # metavar = "pathways.dat",
        type=str,
        help="path to MetaCyc pathways.dat file",
        
    )
    parser.add_argument(
        "-c",
        "--compounds_path",
        # metavar="compound-links.dat",
        type=str,
        help="path to MetaCyc compound-links.dat file containing info per compound",
    )
    parser.add_argument(
        "-f",
        "--filter",
        # metavar = "filter_for_these.sdf",
        required=False,
        default=None,
        help = ("[under_construction] (optional) path to .sdf file with compounds",
                " to filter pathways on (for inclusion). Default: no filter",)
    )
    parser.add_argument(
        "-o",
        "--output_path",
        # metavar="output",
        type=str,
        help="path to output .tsv file. Default: metacyc_chains.tsv",
        default="metacyc_chains.tsv",
    )
    parser.add_argument(
        "-s",
        "--smiles",
        action="store_true",
        default=False,
        help=(
            "[under construction]: output as smiles instead of MetaCyc IDs",
            "(first column will still be MetaCyc pathway IDs)",
        ),
    )
    return parser.parse_args()
    

# ============================== input handling ===============================


def extract_linestartswith(
    lines: list, starts: list = [], remove_starts: bool = True, strip_extra: str = " - "
) -> dict[list[str]]:
    """
    Extracts lines starting with terms in starts, returning them as a list to accomodate multiple values per entry

        Args:
            lines (list): list of strings
            starts (list): list of strings to search for at the start of each line
            remove_starts (bool): whether to remove the start of the line
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
            if line.startswith(start):
                extraction[start].append(
                    line.strip().replace(start, "").replace(strip_extra, "").strip()
                )
    return extraction


def info_from_entries(
    info_loc: str,
    info2extract: list[str],
    ignore_start: str = "#",
    entry_sep="//",
    encoding="iso8859-1",
    remove_starts=True,
    strip_extra: str = " - ",
) -> pd.DataFrame:
    """
    Extracts information from entries in a file

        Args:
            info_loc (str): location of file
            info2extract (list): list of strings to search for at the start of each line
            ignore_start (str): ignore lines starting with this string
            entry_sep (str): separator for entries
            encoding (str): file encoding
            remove_starts (bool): whether to remove the start of the line
            strip_extra (str): string to strip from the start of the line
        Returns:
            pd.DataFrame: dataframe with info2extract as columns

    Remark:
        within a collection of lines, only extracts the lines starting with
        terms in starts, returning them as a list to accomodate multiple values
        per entry
    """
    entries = entry_parser(
        readr(info_loc, ignore_start=ignore_start, encoding=encoding), sep=entry_sep
    )

    all_vals = per_entry(
        entries,
        extract_linestartswith,
        starts=info2extract,
        remove_starts=remove_starts,
        strip_extra=strip_extra,
    )
    df = pd.DataFrame.from_records(all_vals)

    return df


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


# ============================== formatting df =====================================


def _df_get_all_rxn_layouts(
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
    return df


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
    df.columns = df.columns.str.lower()
    df.columns = df.columns.str.replace("-", "_")
    return df


def _remove_empties(
    df: pd.DataFrame, subset: list[str] = ["left", "right"]
) -> pd.DataFrame:
    """Removes empty strings from subset columns"""
    df = df.replace(r"^\s*$", np.nan, regex=True)
    df = df.dropna(subset=subset)
    return df


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


def _listcell_to_strcell(df: pd.DataFrame, colname: str) -> pd.DataFrame:
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
    df = _df_get_all_rxn_layouts(df, column_name="reaction_layout")
    df = _switch_directions(df, subset=["left", "right"], on="direction")
    df = _remove_water(df, subset=["left", "right"])
    df = _remove_empties(df, subset=["left", "right"])
    df = _filter_by_num_reactions(
        df, min_num=min_num, reaction_list_col=reaction_list_col
    )
    df = _listcell_to_strcell(df, "unique_id")
    df = _rename_col(df, "unique_id", "pathway_id")
    return df


# ============================== additional pathway filtering/adjusting =====================================
def mols_per_pathway(pathways: pd.DataFrame) -> pd.DataFrame:
    """
    Returns dataframe with all mols per pathway

        Args:
            pathways (pd.DataFrame): dataframe with pathways
        Returns:
            pd.DataFrame: dataframe with all mols per pathway

    Remark:
        - specific to pathway-tools dataframe format but then lowercased and with _ instead of -
    """
    pw = pathways  # for brevity

    def paste(x):
        return ",".join(x)

    pw["all_right"] = pw.groupby(pw.index)["right"].transform(paste).str.split(",")
    pw["all_left"] = pw.groupby(pw.index)["left"].transform(paste).str.split(",")
    pw["all_mols"] = pw["all_left"] + pw["all_right"]
    pw.drop(columns=["all_left", "all_right"], inplace=True)

    # condense rows on unique_id, then explode for each reagent one row
    pw_mols = pw.groupby(pw.index).agg({"pathway_id": "first", "all_mols": "first"})
    pw_mols = pw_mols.explode("all_mols")
    # get rid of duplicate rows of same pathway_id and mol
    pw_before = pw_mols["pathway_id"].unique().tolist()  # for checking
    pw_mols = pw_mols.drop_duplicates(subset=["pathway_id", "all_mols"])
    pw_after = pw_mols["pathway_id"].unique().tolist()  # for checking
    assert pw_after == pw_before, "pathways lost after unduplicating"
    return pw_mols


def filter_pathwayid(
    df: pd.DataFrame, is_in: list[str], col: str = "pathway_id"
) -> pd.DataFrame:
    """Filters dataframe by pathway_id"""
    return df[df[col].isin(is_in)]


def explode_on_products(df: pd.DataFrame) -> pd.DataFrame:
    """
    Explodes dataframe downwards, to get one product per row

        Args:
            df (pd.DataFrame): dataframe with pathways
        Returns:
            pd.DataFrame: dataframe with one product per row
    Remark:
        - if more than one product, explode the dataframe downwards, to get one product per row
        - this makes same precursors/reactions appear more than once
    """
    # unfold dataframe downwards if there are more than one right or left primary
    exploded = df.copy()
    # first make into list, then explode
    exploded["right"] = exploded["right"].str.split(" ")
    exploded = exploded.explode("right")
    return exploded


# ============================ obtain annotations ============================
# ---------------------------- classes annotation ----------------------------


def cmd_clean_classes() -> None:  # under construction
    """Linux command to clean classes.dat (obsolete)"""
    cmd = 'grep -E "^[A-Za-z]|/{2}" classes.dat | grep -E -v "SYNONYMS" | grep -E -v "COMMENT" | grep -E -v "TYPES" > cleaner_classes.dat '
    subprocess.run(cmd, shell=True)
    return None


def get_idandname(dat_entry: list) -> tuple[str]:
    """Returns idnum and name from dat_entry"""
    idnum, name = "", ""
    for line in dat_entry:
        if line.startswith("UNIQUE-ID"):
            idnum = line.split(" - ")[-1]
        elif line.startswith("COMMON-NAME"):
            name = line.split(" - ")[-1]
    return idnum, name


def get_classes_annotation(filename: str = "../metacyc/cleaner_classes.dat") -> dict:
    """Returns dictionary of classes annotation from filename"""
    annot_entries = entry_parser(
        readr("../metacyc/cleaner_classes.dat", encoding="iso8859-1"), sep="//"
    )
    annot = dictify(per_entry(annot_entries, get_idandname))
    return annot


# ----------------------------------------------------------------------------
# ---------------------------- reaction annotation (obsolete)-----------------
def cmd_clean_rxns() -> None:
    """Linux command to clean reactions.dat (obsolete)"""
    cmd = 'grep -E "^[A-Za-z]|/{2}" reactions.dat | grep -E "^UNI|^LEFT|^RIGHT|^REACTION\-DIRECTION|^/{2}" > cleaner_rxns.dat'
    subprocess.run(cmd, shell=True)
    return None


def get_preandpost(entry: list[str]) -> tuple[str]:
    """Returns pre and post from entry"""
    idnum, name = "", ""
    for line in entry:
        if line.startswith("UNIQUE-ID"):
            idnum = line.split(" - ")[-1]
        elif line.startswith("COMMON-NAME"):
            name = line.split(" - ")[-1]
    return idnum, name


# ----------------------------------------------------------------------------
# ---------------------------- compound annotation ---------------------------


def file_to_compounds(
    filename: str, column_names: list = ["unique_id", "inchi", "SMILES"]
) -> pd.DataFrame:
    """Returns dataframe from filename

    Args:
        filename (str): location of file
        column_names (list): list of column names
    Returns:
        pd.DataFrame: dataframe from filename

    """
    usecols = (i for i in range(len(column_names)))
    array = np.loadtxt(filename, dtype=str, delimiter="\t", usecols=usecols)
    df = pd.DataFrame(array, columns=column_names)
    return df


def _sdf_to_records(
    sdf_path: str,
    column_names: list = ["coconut_id", "inchi", "SMILES", "rdk_inchi", "rdk_SMILES"],
) -> pd.DataFrame:
    """
    Returns records from sdf_path

        Args:
            sdf_path (str): location of file
            column_names (list): list of column names
        Returns:
            pd.DataFrame: dataframe from sdf_path

    """
    sdf_props = []

    supplier = Chem.SDMolSupplier(sdf_path)
    for i, mol in tqdm(enumerate(supplier), total=len(supplier)):
        this_mol_props = {q: "" for q in column_names}
        if mol is None:
            logging.warning("molecule {} could not be read".format(i))
            continue
        else:
            for prop in column_names:
                if mol.HasProp(prop):
                    this_mol_props[prop] = mol.GetProp(prop)
        sdf_props.append(this_mol_props)

    return sdf_props


def sdf_to_compounds(*args, **kwargs) -> pd.DataFrame:
    """Returns dataframe from sdf file with properties as columns"""
    records = _sdf_to_records(*args, **kwargs)
    return pd.DataFrame(records)


def add_partial_inchi(
    df, inchi_column: str = "inchi", inchi_parts: int = 3
) -> pd.DataFrame:
    """
    Adds partial inchi to dataframe

        Args:
            df (pd.DataFrame): dataframe with inchi_column
            inchi_column (str): name of the inchi column
            inchi_parts (int): number of parts to keep
        Returns:
            pd.DataFrame: dataframe with added partial inchi

    Remark:
        - partial inchi refers to the first inchi_parts parts of the inchi divided by /
    """
    newcol = f"{inchi_parts}_{inchi_column}"
    # # create mask for rows with inchi that are strings
    # inchi_present = ~df[inchi_column].isna()
    # # create mask for rows with inchi that are strings
    inchi_present = df[inchi_column].apply(lambda x: isinstance(x, str) and x.startswith("InChI="))
    df[newcol] = df[inchi_present][inchi_column].apply(lambda x: x.split("/"))
    # df[newcol] = df[~df[inchi_column].isna()][inchi_column].str.split("/") # gives error with coco_mols, saying they are not always string values
    df[newcol] = df[inchi_present][newcol].apply(
        lambda x: "/".join(x[:inchi_parts])
    ).astype(str)
    df[newcol] = df[newcol].fillna("")
    return df


def takeover_rdk_partials(
    df: pd.DataFrame, columns: dict = {"3_rdk_inchi": "3_inchi"}
) -> pd.DataFrame:
    """
    Takes over partial inchi from one column to another

        Args:
            df (pd.DataFrame): dataframe with columns
            columns (dict): dictionary of columns
        Returns:
            pd.DataFrame: dataframe with partial inchi taken over

    """
    for key, val in columns.items():
        df[val] = df[val].fillna(df[key])
        df.loc[df[key] == df[val], val] = np.nan
        assert df[df[key] == df[val]].empty, "error in takeover_rdk_partials merging"
    return df


def merge_on_partial_inchi(
    coco: pd.DataFrame, meta: pd.DataFrame, on: str = "3_inchi"
) -> pd.DataFrame:
    """
    Merges dfs on partial inchi

        Args:
            coco (pd.DataFrame): dataframe
            meta (pd.DataFrame): dataframe
            on (str): column name
        Returns:
            pd.DataFrame: dataframe
    """
    df = pd.merge(coco, meta, on=on, how="inner", suffixes=("_coco", "_meta"))
    return df


def merge_on_smiles(
    coco: pd.DataFrame, meta: pd.DataFrame, on: str = "SMILES"
) -> pd.DataFrame:
    """
    Merges dfs on smiles

        Args:
            coco (pd.DataFrame): dataframe
            meta (pd.DataFrame): dataframe
            on (str): column name
        Returns:
            pd.DataFrame: dataframe
    """

    df = pd.merge(coco, meta, on="SMILES", how="inner", suffixes=("_coco", "_meta"))
    return df


def merge_merges(on_inchi: pd.DataFrame, on_smiles: pd.DataFrame) -> pd.DataFrame:
    """
    Merges the two dataframes resulting from merge_on_partial_inchi and merge_on_smiles
    """
    on_smiles.rename(columns={"SMILES": "SMILES_meta"}, inplace=True)
    df = pd.merge(
        on_inchi,
        on_smiles,
        on=["coconut_id", "unique_id", "inchi_meta", "SMILES_meta"],
        how="outer",
        suffixes=("_inchi", "_SMILES"),
    )
    df.rename(columns={"SMILES_meta": "SMILES", "inchi_meta": "inchi"}, inplace=True)
    return df


def get_common_compounds(
    coco_mols: pd.DataFrame, meta_mols: pd.DataFrame,
    inchi_parts: int = 3
) -> pd.DataFrame:
    """
    Gets common compounds from coco and meta mols

        Args:
            coco_mols (pd.DataFrame): dataframe
            meta_mols (pd.DataFrame): dataframe
        Returns:
            pd.DataFrame: dataframe with common compounds' coconut and meta ids, inchi and smiles

    Remark:
        - gets compounds that are in both coco and meta mols
        - filters on partial inchi OR smiles in common
    """
    meta_mols = add_partial_inchi(
        meta_mols, inchi_column="inchi", inchi_parts=inchi_parts
    )
    coco_mols = add_partial_inchi(
        coco_mols, inchi_column="inchi", inchi_parts=inchi_parts
    ) # gives error!
    coco_mols = add_partial_inchi(
        coco_mols, inchi_column="rdk_inchi", inchi_parts=inchi_parts
    )
    coco_mols = takeover_rdk_partials(coco_mols, columns={"3_rdk_inchi": "3_inchi"})

    log_col = f"{inchi_parts}_inchi"
    logging.info(
        f"{len(coco_mols[~coco_mols[log_col].isna()])} coco mols have inchi"
    )  # 406529
    logging.info(
        f"{len(meta_mols[~meta_mols[log_col].isna()])} meta mols have inchi"
    )  # 19173

    on_inchi = merge_on_partial_inchi(coco_mols, meta_mols)
    on_smiles = merge_on_smiles(coco_mols, meta_mols)
    coco_meta = merge_merges(on_inchi, on_smiles)

    logging.info(f"{len(on_inchi)} compounds have a partial inchi in common")
    logging.info(f"{len(on_smiles)} compounds have a smiles in common")
    logging.info(f"{len(coco_meta)} compounds in common")  # 406529
    return coco_meta[["coconut_id", "unique_id", "inchi", "SMILES"]]


# ------------------------ get relevant pathways -----------------------------


def get_pathways_of_interest(
    pathways: pd.DataFrame, compounds_of_interest: pd.DataFrame
) -> pd.DataFrame:
    """
    Returns pathways containing compounds of interest

        Args:
            pathways (pd.DataFrame): dataframe
            common_compounds (pd.DataFrame): dataframe
        Returns:
            pd.DataFrame: dataframe of pathways containing compounds of interest

    Remark:
        - filters pathways on common compounds by "compound_id"

    """
    common_compound_ids = compounds_of_interest["unique_id"].unique() # filters on MetaCyc ID
    pw_reagents = mols_per_pathway(pathways)
    pw_reagents_oi = df_filter(pw_reagents, col="all_mols", is_in=common_compound_ids)
    pw_names_oi = pw_reagents_oi["pathway_id"].unique()
    pw_oi = df_filter(pathways, col="pathway_id", is_in=pw_names_oi)
    return pw_oi


# ----------------------------------------------------------------------------
# ============================================================================


# ============================= traverse pathway =============================
def _row_precursors(row: pd.Series) -> list[str]:
    """returns list of precursors for a given row"""
    return row["left"].split(" ")


def _row_products(row: pd.Series) -> list[str]:
    """returns list of products for a given row"""
    return row["right"].split(" ")


def _is_next_rxn(parent_row: pd.Series, child_row: pd.Series) -> bool:
    """checks if the parent_row is a precursor of the rxn2_prec"""
    for query in _row_products(parent_row):
        if query in _row_precursors(child_row):
            return True
    return False


def get_overlap_mol(parent_row: pd.Series, child_row: pd.Series) -> str:
    """returns overlap molecule between two rows of reactions"""
    overlap_mols = []
    for par_mol in _row_products(parent_row):
        for child_mol in _row_precursors(child_row):
            if par_mol == child_mol:
                overlap_mols.append(par_mol)
    return overlap_mols


def get_beginnings(one_pw: pd.DataFrame):
    """returns list of reactions that are not products of any other reaction"""
    beginning_rxns = []
    all_rights = []
    for right in one_pw["right"].tolist():
        all_rights = all_rights + right.split(" ")

    for _, rxn_row in one_pw.iterrows():
        for precursor in _row_precursors(rxn_row):
            if precursor not in all_rights:
                # add rxn_row index to list
                beginning_rxns.append(rxn_row["reaction_id"])
    return beginning_rxns


def get_child_rxns(parent: str, pw: pd.DataFrame, unused_rxns: list[str]) -> list[str]:
    """for parent compound, returns list of the children reactions"""

    parent = pw[pw["reaction_id"] == parent]
    # get single row of parent as pd.Series
    parent = parent.iloc[0]
    children = []

    unused_pw = pw[pw["reaction_id"].isin(unused_rxns)]
    for _, rxn_row in unused_pw.iterrows():
        # if the parent appears in the reaction's products
        if _is_next_rxn(parent, rxn_row):
            children.append(rxn_row["reaction_id"])
    logging.debug(
        f"searching for{parent['reaction_id']}from\n{unused_pw['reaction_id']}\nfound{children}\n\n"
    )
    return children


# remove for child in children: continue, remove from unused_rxns


def make_tree(parent: str, pw: pd.DataFrame, unused: list[str]) -> list[str]:
    """
    Returns tree of reactions from parent reaction

        Args:
            parent (str): reaction id
            pw (pd.DataFrame): dataframe
            unused (list): list of unused reactions
        Returns:
            list: tree of reactions from parent reaction

    Remark:
        - recursive function
        - uses get_child_rxns to get children of parent
        - uses make_tree to get children of children
    """

    tree = []

    if parent in unused:
        unused.remove(parent)

    # get children of parent
    children = get_child_rxns(parent, pw, unused)
    if not children:
        return tree
        # return parent
    # if reactions preceding this parent, recurse:
    for child in children:
        this_tree = [parent]
        this_unused_rxn = unused.copy()
        this_unused_rxn.remove(child)
        child_tree = make_tree(child, pw, this_unused_rxn)
        if len(child_tree) > 0:
            this_tree.append(child_tree)
        tree.append(this_tree)

    # get any duplicate lists in tree
    for i, item in enumerate(tree):
        if item in tree[:i]:
            tree[i] = None
    tree = [item for item in tree if item is not None]
    if len(tree) == 1:
        tree = tree[0]
    return tree


def _contains_list(lst: list) -> bool:
    """
    Returns True if lst contains a list
    """
    for item in lst:
        if isinstance(item, list):
            return True
    return False


def tree_to_chains(tree):
    """
    Extracts rxn chains from tree of reactions
    """
    chains = []
    if not _contains_list(tree):
        return tree
    elif isinstance(tree[0], str) and _contains_list(tree):
        parent = tree[0]
        children = tree[1]
        for chain in tree_to_chains(children):
            if isinstance(chain, str):
                chain = [chain]
            chains.append([parent] + chain)
    elif isinstance(tree[0], list):
        for lst in tree:
            if isinstance(lst, str):
                lst = [lst]
            for chain in tree_to_chains(lst):
                chains.append(chain)
    return chains


def remove_duplicate_chains(chains: list[list[str]]) -> list[list[str]]:
    """Removes duplicate chains from chains"""
    for i, chain in enumerate(chains):
        if chain in chains[:i]:
            chains[i] = None
    chains = [chain for chain in chains if chain is not None]
    return chains


def remove_subchains(chains: list[list[str]]) -> list[list[str]]:
    """Removes subchains from chains"""
    # check if a chain is a subchain of another chain
    for i, chain in enumerate(chains):
        for j, other_chain in enumerate(chains):
            if i == j:
                continue
            if chain is None or other_chain is None:
                continue
            if ",".join(chain) in ",".join(other_chain):
                chains[i] = None
    chains = [chain for chain in chains if chain is not None]
    return chains


def __get_chain_mols(pw: pd.DataFrame, rxn_chain: list[str]) -> list[str]:
    """returns list of mols for a chain of reactions"""
    chain_of_mols = []

    for i, rxn_id in enumerate(rxn_chain):
        rxn_row = pw.loc[pw["reaction_id"] == rxn_id]
        if rxn_row.empty:
            continue
        rxn_row = rxn_row.iloc[0]
        if i == 0:
            chain_of_mols.append(_row_precursors(rxn_row))
            continue
        prev_rxn = pw.loc[pw["reaction_id"] == rxn_chain[i - 1]].iloc[0]
        chain_of_mols.append(get_overlap_mol(prev_rxn, rxn_row))
    return chain_of_mols


def __choose_first_mols(chain_of_mols: list[list[str]]) -> list[str]:
    """returns list of mols for a chain of reactions"""
    return [mols[0] for mols in chain_of_mols]


def _pw_to_chains_mols(one_pw: pd.DataFrame) -> list[list[list[str]]]:
    """
    Returns list of chains of mols for a pathway

        Args:
            one_pw (pd.DataFrame): dataframe containing one pathway
        Returns:
            list: list of chains of mols for a pathway with [mols for rxn1, mols for rxn2, ...]

    Remark:
        - uses make_tree to get chains of reactions
        - uses __get_chain_mols to get mols for a chain of reactions
        - uses __choose_first_mols to choose first mols for a chain of reactions
        - sorts chains by length
    """
    chains = []
    chains_mols = []
    all_rxns = one_pw["reaction_id"].tolist()
    beginnings = get_beginnings(one_pw)
    for beginning in beginnings:
        tree = make_tree(beginning, one_pw, all_rxns)
        chains = chains + tree_to_chains(tree)
    for chain in chains:
        chain_of_mols = __get_chain_mols(one_pw, chain)
        chain_of_mols = __choose_first_mols(chain_of_mols)
        chains_mols.append(chain_of_mols)
    chains_mols = remove_duplicate_chains(chains_mols)
    chains_mols = remove_subchains(chains_mols)
    return sorted(chains_mols, key=len, reverse=True)


def chains_per_pathway(pw: pd.DataFrame) -> pd.DataFrame:
    """Returns dataframe with chains of mols per pathway"""
    pw_chains = {}
    for pw_id in pw["pathway_id"].unique():
        chains_mols = _pw_to_chains_mols(pw[pw["pathway_id"] == pw_id])
        if len(chains_mols) > 10:
            chains_mols = chains_mols[:10]
        pw_chains[pw_id] = chains_mols

    pw_chains = pd.DataFrame.from_dict(pw_chains, orient="index")
    pw_chains["chains"] = pw_chains.values.tolist()
    pw_chains = pw_chains[["chains"]]
    # explode on chains
    pw_chains = pw_chains.explode("chains")
    # drop na
    pw_chains = pw_chains.dropna()

    # get length of longest chain
    pw_chains["chain_len"] = pw_chains["chains"].apply(lambda x: len(x))
    max_len = pw_chains["chain_len"].max()
    # filter out chains that are too short
    pw_chains = pw_chains[pw_chains["chain_len"] > 4]
    # drop chain_len column
    pw_chains.drop(columns=["chain_len"], inplace=True)

    # split chains into columns
    pw_chains = pw_chains["chains"].apply(pd.Series)

    # keep only first 15 columns
    pw_chains = pw_chains.iloc[:, :15]

    return pw_chains


def compound_id_to_structure(
    df: pd.DataFrame, dict1: dict, dict2: dict, dict3: dict
):  # under construction
    """
    Returns dataframe with structure for compound_id from dict1, dict2, dict3
    """
    newdf = df.copy()
    # fill all cells as nan
    newdf = newdf.applymap(lambda x: np.nan)

    for col in df.columns:
        newdf[col] = df[col].map(dict1)
        newdf[col] = newdf[col].fillna(df[col].map(dict2))
        newdf[col] = newdf[col].fillna(df[col].map(dict3))
    return newdf


# ============================ df transformations =============================


def df_filter(df: pd.DataFrame, col: str, is_in: list) -> pd.DataFrame:
    """Filters dataframe by col"""
    mask = df[col].isin(is_in)
    df_of_interest = df[mask]
    return df_of_interest


def get_indexes(df: pd.DataFrame) -> list[int]:
    """Returns list of indexes from dataframe"""
    unique_index = df.index.unique()
    return unique_index


def remove_cols(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    """Removes columns from dataframe"""
    new = df
    for i in cols:
        new = new.drop(columns=i).drop_duplicates()
    return new


def annotate(df: pd.DataFrame, replace_dict: dict) -> pd.DataFrame:
    """Annotates dataframe with replace_dict"""
    return df.replace(to_replace=replace_dict)


# def filter_by_rxn_len(df: pd.DataFrame, length: int = 4) -> pd.DataFrame:
#     mask = np.array([len(x) > 3 for x in colval_per_index(df)])
#     all_ind = np.array(df.index.unique())
#     long_enough = all_ind[mask]
#     new = df[df.index.isin(long_enough)]
#     return new


def mergeinchis(df: pd.DataFrame) -> pd.DataFrame:
    """Merges inchi columns in dataframe (obsolete)"""
    newdf = df
    newdf["InChi"] = newdf["INCHI"].fillna(newdf["NON-STANDARD-INCHI"])
    return newdf


def to_conversion_dict(df: pd.DataFrame, allcapskeys: bool = True) -> dict:
    """
    Makes dictionary from dataframe for conversion

        Args:
            df (pd.DataFrame): dataframe
            allcapskeys (bool): whether to make keys all caps
        Returns:
            dict: dictionary for conversion

    Remark:
        - fills inchi, smiles and non-standard-inchi columns
    """
    ndf = df
    ndf.columns = ndf.columns.str.upper()
    ndf.columns = ndf.columns.str.replace("_", "-")

    # ndf["filled"] = ndf["INCHI"].fillna(ndf["smiles"])#.fillna(ndf["non-standard-inchi"])
    ndf["filled"] = ndf["SMILES"]#.fillna(ndf["INCHI"]).fillna(ndf["NON-STANDARD-INCHI"])
    ndf["filled"].fillna("")
    dic = pd.Series(ndf["filled"].values, index=ndf["UNIQUE-ID"]).to_dict()
    full_dic = {}
    for key, val in dic.items():
        if isinstance(val, str) and val:
            if allcapskeys:
                full_dic[key.upper()] = val
            else:
                full_dic[key] = val
    return full_dic



# ============================================================================


def main():
    args = cli()

    # pathways_path = "/Users/lucina-may/thesis/metacyc/pathways.dat"
    
    pathways = get_pathways(args.pathways_path)
    pathways = clean_pw_df(pathways)

    meta_mols = file_to_compounds(args.compounds_path)

    if isinstance(args.filter, str):
        # create temporary save directory to not have to redo the filtering
        tmp_dir = "metacyc_extract_temp"
        if not os.path.exists(tmp_dir):
            os.makedirs(tmp_dir)
        # common_compound_path = "/Users/lucina-may/thesis/metacyc/coco_meta.tsv"
        common_compound_path = os.path.join(tmp_dir, "coco_meta.tsv")
        if not os.path.exists(common_compound_path):
            # coco_mols_path = "/Users/lucina-may/thesis/metacyc/coconut-links.tsv"
            # meta_mols_path = "/Users/lucina-may/thesis/metacyc/compound-links.dat"
            # meta_mols_path = args.compounds_path
            coco_mols_path = os.path.join(tmp_dir, "coconut-links.tsv")

            if not os.path.exists(coco_mols_path):
                # coco_sdf = "/Users/lucina-may/thesis/input/coconut.sdf"
                coco_sdf = args.filter
                coco_mols = sdf_to_compounds(sdf_path=coco_sdf)  # standard query
                coco_mols.to_csv(coco_mols_path, sep="\t", index=False)

            coco_mols = pd.read_csv(coco_mols_path, sep="\t")
            # meta_mols = file_to_compounds(meta_mols_path) # moved out of if
            compounds = get_common_compounds(coco_mols, meta_mols, inchi_parts=3)

            compounds.to_csv(common_compound_path, sep="\t", index=False)
        compounds = pd.read_csv(common_compound_path, sep="\t")

        pathways_oi = get_pathways_of_interest(pathways, compounds)
    else:
        pathways_oi = pathways

    pw_chains = chains_per_pathway(pathways_oi)
    # change unique_ids to smiles with moldict
    if args.smiles:
        mol_dict = to_conversion_dict(meta_mols)
        pw_chains_smiles = pw_chains.replace(mol_dict)
        pw_chains_smiles.to_csv("metacyc_chains_smiles.tsv", sep="\t", index=True)


    pw_chains.to_csv("metacyc_chains.tsv", sep="\t", index=True)
    exit()


if __name__ == "__main__":
    main()
