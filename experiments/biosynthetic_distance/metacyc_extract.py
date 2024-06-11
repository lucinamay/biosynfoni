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
        # metavar="compounds.dat",
        type=str,
        help="path to MetaCyc compounds.dat file containing info per compound",
    )
    parser.add_argument(
        "-f",
        "--filter",
        # metavar = "compounds_to_filter_for.sdf",
        required=False,
        default=None,
        help = ("(optional) path to .sdf file with compounds",
                " to filter pathways on (for inclusion). Default: no filter",)
    )
    parser.add_argument(
        "-o",
        "--output_path",
        # metavar="output",
        type=str,
        help="path to output .tsv file. Default: metacyc_chains.tsv",
        default="metacyc_chains",
    )
    parser.add_argument(
        "-s",
        "--smiles",
        action="store_true",
        default=False,
        help=(
            "also outputs as smiles instead of only as MetaCyc IDs",
            "(first column will still be MetaCyc pathway IDs)",
        ),
    )
    args = parser.parse_args()
    args.pathways_path = os.path.abspath(args.pathways_path)
    args.compounds_path = os.path.abspath(args.compounds_path)
    if args.filter:
        args.filter = os.path.abspath(args.filter)
    args.output_path = os.path.abspath(args.output_path)
    return args
    

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
            if line.startswith(f"{start} "): # to avoid ones not exactly starting with start
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

def _str_to_valid_mol(structure_string: str, clean: bool = False) -> Chem.Mol:
    mol = ""
    if structure_string.startswith("InChI="):
        mol = Chem.MolFromInchi(structure_string)
    elif clean:
        mol = Chem.MolFromSmiles(
            structure_string.split(
                "[a ",
            )[0]
        ) # remove any non-specific parts of the smiles string
    else:
        mol = Chem.MolFromSmiles(structure_string)
    if mol:
        return mol
    else:
        return ""

def write_compounds(df: pd.DataFrame) -> dict:
    """
    Makes valid dataframe for conversion of UNIQUE-ID to INCHI and SMILES

        Args:
            df (pd.DataFrame): dataframe
            allcapskeys (bool): whether to make keys all caps
        Returns:
            dict: dictionary for conversion

    Remark:
        - fills inchi, smiles and non-standard-inchi columns
    """
    df.columns = df.columns.str.upper()
    df.columns = df.columns.str.replace("_", "-")
    df = df.set_index("UNIQUE-ID")
    if "NON-STANDARD-INCHI" not in df.columns:
        df["NON-STANDARD-INCHI"] = np.nan


    df["mol"] = df["INCHI"].fillna(df["SMILES"]).fillna(df["NON-STANDARD-INCHI"]).fillna(np.nan)
    #print mol values if isinstance list
    # df["mol"].apply(lambda x: print(x) if isinstance(x, list) else x)

    df["mol"] = df["mol"].apply(lambda x: _str_to_valid_mol(x) if isinstance(x, str) else np.nan)
    df["INCHI"] = df["mol"].apply(lambda x: Chem.MolToInchi(x) if isinstance(x, Chem.Mol) else np.nan)
    df["SMILES"] = df["mol"].apply(lambda x: Chem.MolToSmiles(x) if isinstance(x, Chem.Mol) else np.nan)
    df.drop(columns=["NON-STANDARD-INCHI", "mol"], inplace=True)

    df.columns = df.columns.str.lower()
    df.columns = df.columns.str.replace("-", "_")
    df.index.name = "unique_id"
    
    df.to_csv("meta_compounds.tsv", sep="\t", index=True)
    return df


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


# ----------------------------------------------------------------------------
def sdf_to_compounds(sdf_path: str) -> pd.DataFrame:
    """
    Returns dataframe from sdf file with properties as columns

        Args:
            sdf_path (str): location of file
    """
    supplier = Chem.SDMolSupplier(sdf_path)
    # df = pd.DataFrame(len(supplier), columns=["coconut_id", "inchi", "smiles"])
    sdf_props = []
    for mol in tqdm(supplier, total=len(supplier), desc="reading sdf"):
        if not mol:
            continue

        this_mol_props = {"coconut_id": np.nan}
        if mol.HasProp("coconut_id"):
            coco_id = mol.GetProp("coconut_id")
        this_mol_props["coconut_id"] = coco_id
        this_mol_props["inchi"] = Chem.MolToInchi(mol)
        this_mol_props["smiles"] = Chem.MolToSmiles(mol)
        sdf_props.append(this_mol_props)

    df = pd.DataFrame(sdf_props)
    df.index = df["coconut_id"]
    return df

def _partial_inchi(inchi:str, parts:int = 3)->str:
    """returns partial inchi, where parts is the number of parts to keep. 3 is up to and including atom connectivity"""
    parts = inchi.split("/")
    return "/".join(parts[:3])


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
        on=["coconut_id", "unique_id", "inchi_meta", "smiles"],
        how="outer",
        suffixes=("_inchi", "_smiles"),
    )
    df.rename(columns={"inchi_meta": "inchi"}, inplace=True)
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
    #get the index as a column
    meta_mols = meta_mols.reset_index()
    coco_mols = coco_mols.reset_index()

    p_i = "partial_inchi"
    meta_mols[p_i] = meta_mols["inchi"].apply(lambda x: _partial_inchi(x, parts=inchi_parts) if isinstance(x, str) else np.nan)
    coco_mols[p_i] = coco_mols["inchi"].apply(lambda x: _partial_inchi(x, parts=inchi_parts) if isinstance(x, str) else np.nan)


    logging.info(
        (f"{len(coco_mols[~coco_mols[p_i].isna()])} coco mols have inchi",
        f"{len(meta_mols[~meta_mols[p_i].isna()])} meta mols have inchi"))  # 406529
  
    same_p_inchi = pd.merge(coco_mols, meta_mols, on=p_i, how="inner", suffixes=("_coco", "_meta"))
    same_smiles = pd.merge(coco_mols, meta_mols, on="smiles", how="inner", suffixes=("_coco", "_meta"))
    coco_meta = merge_merges(same_p_inchi, same_smiles)

    meta_mols = meta_mols.set_index("unique_id")
    coco_mols = coco_mols.set_index("coconut_id")

    logging.info(f"{len(same_p_inchi)} compounds have a partial inchi in common")
    logging.info(f"{len(same_smiles)} compounds have a smiles in common")
    logging.info(f"{len(coco_meta)} compounds in common")  # 406529

    print(f"{len(same_p_inchi)} compounds have a partial inchi in common")
    print(f"{len(same_smiles)} compounds have a smiles in common")
    print(f"{len(coco_meta)} compounds in common")  # 406529

    return coco_meta[["coconut_id", "unique_id", "inchi", "smiles"]]


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
    # return df[df[col].isin(is_in)] @todo: check if this is correct
    mask = df[col].isin(is_in)
    df_of_interest = df[mask]
    return df_of_interest

# ============================================================================


def main():
    args = cli()

    output_dir = os.path.abspath(os.path.dirname(args.output_path))
    temp_dir = os.path.join(output_dir, "intermediate_files")
    iwd = os.getcwd()
    os.makedirs(temp_dir, exist_ok=True)
    os.chdir(temp_dir)

    pathways = get_pathways(args.pathways_path)
    pathways = clean_pw_df(pathways)

    meta_mols = get_compounds(args.compounds_path)
    for col in meta_mols.columns:
        meta_mols = _listcell_to_strcell(meta_mols, col)
    write_compounds(meta_mols)
    meta_mols = pd.read_csv("meta_compounds.tsv", sep="\t", index_col="unique_id")

    if isinstance(args.filter, str):
        common_compound_path = os.path.join(temp_dir, "coco_meta.tsv")
        if not os.path.exists(common_compound_path):
            coco_mols_path = os.path.join(temp_dir, "coconut_compounds.tsv")
            if not os.path.exists(coco_mols_path):
                coco_sdf = args.filter
                coco_mols = sdf_to_compounds(sdf_path=coco_sdf)  # standard query
                coco_mols.to_csv(coco_mols_path, sep="\t", index=True)

            coco_mols = pd.read_csv(coco_mols_path, sep="\t", index_col="coconut_id")
            compounds = get_common_compounds(coco_mols, meta_mols, inchi_parts=3)
            compounds.to_csv(common_compound_path, sep="\t", index=False)
        compounds = pd.read_csv(common_compound_path, sep="\t")

        pathways_oi = get_pathways_of_interest(pathways, compounds)
    else:
        pathways_oi = pathways

    pw_chains = chains_per_pathway(pathways_oi)
    # change unique_ids to smiles with moldict
    
    # go to parent directory
    os.chdir(output_dir)
    if args.smiles or args.filter:
        mol_dict = pd.Series(meta_mols["smiles"].values, index=meta_mols.index).to_dict()
        #@todo: add in all unique_ids that are in pathways.dat but not in compounds.dat
        # mol_dict = to_conversion_dict(meta_mols)
        # pd.DataFrame.from_dict(mol_dict, orient="index").to_csv(f"{args.output_path}_moldict.tsv", sep="\t", index=True)
        pw_chains_structs = pw_chains.replace(mol_dict)
        pw_chains_structs.to_csv(f"{args.output_path}_smiles.tsv", sep="\t", index=True)


    pw_chains.to_csv(f'{args.output_path}.tsv', sep="\t", index=True)
    
    os.chdir(iwd)
    exit(0)


if __name__ == "__main__":
    main()
