#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 16:15:34 2023

@author: lucina-may
"""
import sys, subprocess, os, logging
from sys import argv
from functools import partial

import numpy as np
import pandas as pd
from tqdm import tqdm
from rdkit import Chem

sys.path.append("../src/")
from biosynfoni.inoutput import *

# ============================== input handling ===============================


def extract_linestartswith(
    lines: list, starts: list = [], remove_starts: bool = True, strip_extra: str = " - "
) -> dict[list[str]]:
    """within a collection of lines, only extracts the lines starting with
    terms in starts, returning them as a list to accomodate multiple values
    per entry"""
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
    """non-general function to split the reaction-layout column into reaction-id, left, direction, right"""
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
    """switches left and right columns if direction is R2L"""
    left, right = subset
    df[left], df[right] = np.where(
        df[on] == "R2L", (df[right], df[left]), (df[left], df[right])
    )
    df[on] = df[on].str.replace("R2L", "L2R")
    return df


def _lower_dunder(df: pd.DataFrame) -> pd.DataFrame:
    df.columns = df.columns.str.lower()
    df.columns = df.columns.str.replace("-", "_")
    return df


def _remove_empties(
    df: pd.DataFrame, subset: list[str] = ["left", "right"]
) -> pd.DataFrame:
    df = df.replace(r"^\s*$", np.nan, regex=True)
    df = df.dropna(subset=subset)
    return df


def _remove_water(
    df: pd.DataFrame, subset: list[str] = ["left", "right"], water_str: str = "WATER"
) -> pd.DataFrame:
    for colname in subset:
        df[colname] = df[colname].str.replace(water_str, "")
    return df


def _filter_by_num_reactions(
    df: pd.DataFrame, min_num: int = 4, reaction_list_col="reaction_list"
) -> pd.DataFrame:
    df = df[df[reaction_list_col].str.len() >= min_num]
    return df


def _listcell_to_strcell(df: pd.DataFrame, colname: str) -> pd.DataFrame:
    assert df[df[colname].str.len() > 1].empty, "error, multiple values in cell"
    df[colname] = df[colname].str[0]
    return df


def _rename_col(df: pd.DataFrame, old: str, new: str) -> pd.DataFrame:
    df.rename(columns={old: new}, inplace=True)
    return df


def clean_pw_df(
    df: pd.DataFrame, min_num: int = 4, reaction_list_col="reaction_list"
) -> pd.DataFrame:
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
    return df[df[col].isin(is_in)]


def explode_on_products(df: pd.DataFrame) -> pd.DataFrame:
    """if more than one product, explode the dataframe downwards, to get one product per row.
    this makes same precursors/reactions appear more than once"""
    # unfold dataframe downwards if there are more than one right or left primary
    exploded = df.copy()
    # first make into list, then explode
    exploded["right"] = exploded["right"].str.split(" ")
    exploded = exploded.explode("right")
    return exploded


# ============================ obtain annotations ============================
# ---------------------------- classes annotation ----------------------------


def cmd_clean_classes() -> None:  # under construction
    cmd = 'grep -E "^[A-Za-z]|/{2}" classes.dat | grep -E -v "SYNONYMS" | grep -E -v "COMMENT" | grep -E -v "TYPES" > cleaner_classes.dat '
    subprocess.run(cmd, shell=True)
    return None


def get_idandname(dat_entry: list) -> tuple[str]:
    idnum, name = "", ""
    for line in dat_entry:
        if line.startswith("UNIQUE-ID"):
            idnum = line.split(" - ")[-1]
        elif line.startswith("COMMON-NAME"):
            name = line.split(" - ")[-1]
    return idnum, name


def get_classes_annotation(filename: str = "../metacyc/cleaner_classes.dat") -> dict:
    annot_entries = entry_parser(
        readr("../metacyc/cleaner_classes.dat", encoding="iso8859-1"), sep="//"
    )
    annot = dictify(per_entry(annot_entries, get_idandname))
    return annot


# ----------------------------------------------------------------------------
# ---------------------------- reaction annotation (obsolete)-----------------
def cmd_clean_rxns() -> None:
    cmd = 'grep -E "^[A-Za-z]|/{2}" reactions.dat | grep -E "^UNI|^LEFT|^RIGHT|^REACTION\-DIRECTION|^/{2}" > cleaner_rxns.dat'
    subprocess.run(cmd, shell=True)
    return None


def get_preandpost(entry: list[str]) -> tuple[str]:
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
    usecols = (i for i in range(len(column_names)))
    array = np.loadtxt(filename, dtype=str, delimiter="\t", usecols=usecols)
    df = pd.DataFrame(array, columns=column_names)
    return df


def _sdf_to_records(
    sdf_path: str,
    column_names: list = ["coconut_id", "inchi", "SMILES", "rdk_inchi", "rdk_SMILES"],
) -> pd.DataFrame:
    sdf_path = "/Users/lucina-may/thesis/input/coconut.sdf"
    coco_props = []

    supplier = Chem.SDMolSupplier(sdf_path)
    for i, mol in tqdm(enumerate(supplier), total=len(supplier)):
        this_mol_props = {q: "" for q in column_names}
        if mol is None:
            print("molecule {} could not be read".format(i))
            continue
        else:
            for prop in column_names:
                if mol.HasProp(prop):
                    this_mol_props[prop] = mol.GetProp(prop)
        coco_props.append(this_mol_props)
    return coco_props


def sdf_to_compounds(*args, **kwargs) -> pd.DataFrame:
    records = _sdf_to_records(*args, **kwargs)
    return pd.DataFrame(records)


def add_partial_inchi(
    df, inchi_column: str = "inchi", inchi_parts: int = 3
) -> pd.DataFrame:
    newcol = f"{inchi_parts}_{inchi_column}"
    df[newcol] = df[~df[inchi_column].isna()][inchi_column].str.split("/")
    df[newcol] = df[~df[inchi_column].isna()][newcol].apply(
        lambda x: "/".join(x[:inchi_parts])
    )
    return df


def takeover_rdk_partials(
    df: pd.DataFrame, columns: dict = {"3_rdk_inchi": "3_inchi"}
) -> pd.DataFrame:
    for key, val in columns.items():
        df[val].fillna(df[key], inplace=True)
        df.loc[df[key] == df[val], val] = np.nan
        assert df[df[key] == df[val]].empty, "error in takeover_rdk_partials merging"
    return df


def merge_on_partial_inchi(
    coco: pd.DataFrame, meta: pd.DataFrame, on: str = "3_inchi"
) -> pd.DataFrame:
    df = pd.merge(coco, meta, on=on, how="inner", suffixes=("_coco", "_meta"))
    return df


def merge_on_smiles(
    coco: pd.DataFrame, meta: pd.DataFrame, on: str = "SMILES"
) -> pd.DataFrame:
    df = pd.merge(coco, meta, on="SMILES", how="inner", suffixes=("_coco", "_meta"))
    return df


def merge_merges(on_inchi: pd.DataFrame, on_smiles: pd.DataFrame) -> pd.DataFrame:
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
    coco_mols: pd.DataFrame, meta_mols: pd.DataFrame
) -> pd.DataFrame:
    inchi_parts = 3
    meta_mols = add_partial_inchi(
        meta_mols, inchi_column="inchi", inchi_parts=inchi_parts
    )
    coco_mols = add_partial_inchi(
        coco_mols, inchi_column="inchi", inchi_parts=inchi_parts
    )
    coco_mols = add_partial_inchi(
        coco_mols, inchi_column="rdk_inchi", inchi_parts=inchi_parts
    )
    coco_mols = takeover_rdk_partials(coco_mols, columns={"3_rdk_inchi": "3_inchi"})

    log_col = "3_inchi"
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
    pathways: pd.DataFrame, common_compounds: pd.DataFrame
) -> pd.DataFrame:
    common_compound_ids = common_compounds["compound_id"].unique()
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
    overlap_mols = []
    for par_mol in _row_products(parent_row):
        for child_mol in _row_precursors(child_row):
            if par_mol == child_mol:
                overlap_mols.append(par_mol)
    return overlap_mols


def get_beginnings(one_pw: pd.DataFrame):
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


def _list_present(lst: list) -> bool:
    for item in lst:
        if isinstance(item, list):
            return True
    return False


def tree_to_chains(tree):
    chains = []
    if not _list_present(tree):
        return tree
    elif isinstance(tree[0], str) and _list_present(tree):
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
    for i, chain in enumerate(chains):
        if chain in chains[:i]:
            chains[i] = None
    chains = [chain for chain in chains if chain is not None]
    return chains


def remove_subchains(chains: list[list[str]]) -> list[list[str]]:
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


def __choose_first_mols(chain_of_mols) -> list[str]:
    return [mols[0] for mols in chain_of_mols]


def _pw_to_chains_mols(one_pw: pd.DataFrame) -> list[list[list[str]]]:
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
    mask = df[col].isin(is_in)
    df_of_interest = df[mask]
    return df_of_interest


def get_indexes(df: pd.DataFrame) -> list[int]:
    unique_index = df.index.unique()
    return unique_index


def remove_cols(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    new = df
    for i in cols:
        new = new.drop(columns=i).drop_duplicates()
    return new


def annotate(df: pd.DataFrame, replace_dict: dict) -> pd.DataFrame:
    return df.replace(to_replace=replace_dict)


# def filter_by_rxn_len(df: pd.DataFrame, length: int = 4) -> pd.DataFrame:
#     mask = np.array([len(x) > 3 for x in colval_per_index(df)])
#     all_ind = np.array(df.index.unique())
#     long_enough = all_ind[mask]
#     new = df[df.index.isin(long_enough)]
#     return new


def mergeinchis(df: pd.DataFrame) -> pd.DataFrame:
    newdf = df
    newdf["InChi"] = newdf["INCHI"].fillna(newdf["NON-STANDARD-INCHI"])
    return newdf


def to_conversion_dict(df: pd.DataFrame, allcapskeys: bool = True) -> dict:
    ndf = df
    ndf["filled"] = ndf["INCHI"].fillna(ndf["SMILES"]).fillna(ndf["NON-STANDARD-INCHI"])
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
    pathways_path = "/Users/lucina-may/thesis/metacyc/pathways.dat"
    pathways = get_pathways(pathways_path)
    pathways = clean_pw_df(pathways)

    common_compound_path = "/Users/lucina-may/thesis/metacyc/coco_meta.tsv"
    if not os.path.exists(common_compound_path):
        coco_mols_path = "/Users/lucina-may/thesis/metacyc/coconut-links.tsv"
        meta_mols_path = "/Users/lucina-may/thesis/metacyc/compound-links.dat"

        if not os.path.exists(coco_mols_path):
            coco_sdf = "/Users/lucina-may/thesis/input/coconut.sdf"
            coco_mols = sdf_to_compounds(sdf_path=coco_sdf)  # standard query
            coco_mols.to_csv(coco_mols_path, sep="\t", index=False)

        coco_mols = pd.read_csv(coco_mols_path, sep="\t")
        meta_mols = file_to_compounds(meta_mols_path)
        compounds = get_common_compounds(coco_mols, meta_mols)

        compounds.to_csv(common_compound_path, sep="\t", index=False)
    compounds = pd.read_csv(common_compound_path, sep="\t")

    pathways_oi = get_pathways_of_interest(pathways, compounds)

    pw_chains = chains_per_pathway(pathways_oi)

    pw_chains.to_csv("metacyc_chains.tsv", sep="\t", index=True)

    exit()


if __name__ == "__main__":
    main()
