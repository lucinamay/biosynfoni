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
    df[colname] = df[colname].str.join[0]
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
    parent = pw["reaction_id" == parent]
    children = []

    for rxn_id in unused_rxns:
        current_row = pw[pw["reaction_id"] == rxn_id]
        # if the parent appears in the reaction's products
        if _is_next_rxn(parent, current_row):
            children.append(rxn_id)
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


# later, remove unfit molecules within chain
# def get_chain_mols(pw: pd.DataFrame, rxn_chain: list[str]) -> list[str]:
#     chain_of_mols = []
#     for i, rxn_id in enumerate(rxn_chain):
#         rxn = pw[pw["reaction_id"] == rxn_id]
#         if i == 0:
#             chain_of_mols.append(_row_precursors(rxn))
#             continue
#         prev_rxn = pw[pw["reaction_id"] == rxn_chain[i - 1]]
#         chain_of_mols.append(get_overlap_mol(prev_rxn, rxn))
#     return chain_of_mols


def get_chain_mols(pw: pd.DataFrame, rxn_chain: list[str]) -> list[str]:
    chain_of_mols = []

    for i, (_, rxn_row) in enumerate(pw[pw["reaction_id"].isin(rxn_chain)].iterrows()):
        # for i, rxn_id in enumerate(rxn_chain):
        if i == 0:
            chain_of_mols.append(_row_precursors(rxn_row))
            continue
        prev_rxn = pw[pw["reaction_id"] == rxn_chain[i - 1]].iloc[0]
        chain_of_mols.append(get_overlap_mol(prev_rxn, rxn_row))
    return chain_of_mols


def choose_first_mols(chain_of_mols) -> list[str]:
    return [mols[0] for mols in chain_of_mols]


def _pw_to_chains_mols(one_pw: pd.DataFrame) -> list[list[list[str]]]:
    chains = []
    chains_mols = []
    all_rxns = one_pw["reaction_id"].tolist()
    beginnings = get_beginnings(one_pw)
    for beginning in beginnings:
        tree = make_tree(beginning, one_pw, all_rxns)
        chains.append(tree_to_chains(tree))
    for chain in chains:
        chain_of_mols = get_chain_mols(one_pw, chain)
        chain_of_mols = choose_first_mols(chain_of_mols)
        chains_mols.append(chain_of_mols)
    return chains_mols


def chains_per_pathway(pw: pd.DataFrame) -> pd.DataFrame:
    pw_chains = {}
    for pw_id in pw["pathway_id"].unique():
        chains_mols = _pw_to_chains_mols(pw[pw["pathway_id"] == pw_id])
        pw_chains[pw_id] = chains_mols
    pw_chains = pd.DataFrame.from_dict(pw_chains, orient="index", columns=["chains"])

    # explode on chains
    pw_chains = pw_chains.explode("chains")

    # split chains list into individual columns for each list item
    pw_chains = pw_chains.join(pd.DataFrame(pw_chains["chains"].tolist()))  # test!

    return pw_chains


# ----------------------------------------------------------------------------
# ------------------------------ find start point ----------------------------


# def find_ancestor_rxn(split_rxns: list[list[str]]) -> list[list[str]]:
#     """with a list [pre, post]"""
#     # find oldest ancestor (parent) reaction
#     parents = []
#     all_pre = [item[0] for item in split_rxns]
#     all_post = [item[1] for item in split_rxns]

#     split_pre = []
#     for pres in all_pre:
#         splitted = pres.split(" ")
#         if isinstance(splitted, list):
#             split_pre = split_pre + splitted
#         elif isinstance(splitted, str):
#             split_pre.append(splitted)

#     for i in range(len(all_post)):
#         if len(all_post[i].split(" ")) > 1:
#             continue
#         elif all_post[i] in split_pre:
#             continue
#         else:
#             parents.append(split_rxns[i])
#     return parents


# def check_ancestor_res(ancestor_rxns_res: list[list[str]]) -> str:
#     parent = ""

#     if ancestor_rxns_res:
#         # select the reaction from list
#         parent_rxn = ancestor_rxns_res[0]
#         assert len(parent_rxn) == 2, "error in the parent rxn result"
#         # check if the result reaction exists and
#         if parent_rxn[1]:
#             # print(parent_rxn[1])
#             if len(parent_rxn[1].split(" ")) == 1:
#                 parent = parent_rxn[1]  # selected product
#     return parent


# def get_checked_ancestor(prec_prods: list[list[str]]) -> tuple[str, str]:
#     """gets parent compound of pathway, and returns if its a synthesis
#     or a degradation type pathway"""
#     parent_type = "synthesis"

#     synthesis_par_rxns = find_ancestor_rxn(prec_prods)
#     # check if there is a result and if there is only one result rxn
#     parent = check_ancestor_res(synthesis_par_rxns)
#     if not parent:
#         parent_type = "degradation"
#         # try to find the degradation parent
#         degradation_par_rxns = find_ancestor_rxn(inverse_rxns(prec_prods))
#         parent = check_ancestor_res(degradation_par_rxns)
#         if not parent:
#             parent_type = "first"
#             print(
#                 "Could not find parent",
#                 "defaulted to first product of first rxn:",
#                 prec_prods[0][1],
#             )
#             try:
#                 parent = prec_prods[0][1].split(" ")[0]
#             except:
#                 print("Failed getting 1st product of 1st rxn, return empty")
#                 parent = ""
#                 parent_type = ""
#     return parent, parent_type


# # ----------------------------------------------------------------------------
# # --------------------------------- tree-making ------------------------------


# def get_child_rxn(parent: list[str], rxn_layout: list[list[str]]):  # obsolete?
#     """par_child = [parent]
#     child_found = False
#     for child_i in range(len(rxn_layout)):
#         print(rxn_layout[child_i][1], parent[0])
#         if rxn_layout[child_i][1] == parent[0]:
#             child_found = True
#             print('True')
#             grandchildren = get_child_rxn(rxn_layout[child_i],rxn_layout)
#             par_child.append(rxn_layout[child_i],[grandchildren])
#     if not child_found:
#         return [parent]
#     return par_child"""
#     return None


# def getchild_rxns(parent: str, prec_prods: list[list[str]]) -> list[list[str]]:
#     """for parent compound, returns list of the children reactions"""
#     children = []
#     for reaction in range(len(prec_prods)):
#         # if the parent appears in the reaction's products
#         if parent in prec_prods[reaction][1].split(" "):
#             # append the reactions [precursors, products]
#             children.append(prec_prods[reaction])
#     return children


# def make_tree_dict(
#     parent: str, pw_rxns: list[list[str]], parent_lvl: int = 0, max_lvl: int = 5
# ) -> dict[int, str, dict]:
#     tree_dict = {}  # total dictionary
#     tree_dict["level"] = parent_lvl
#     tree_dict["compound"] = parent

#     # get children of parent
#     children_rxns = getchild_rxns(parent, pw_rxns)
#     child_lvl = parent_lvl + 1

#     # check level to prevent overrecursion
#     if child_lvl > max_lvl:
#         # print('stop at lvl{}: {}'.format(max_lvl,parent))
#         return tree_dict

#     # if no reactions leading to parent's compound, return the tree dict
#     if not children_rxns:
#         return tree_dict

#     # if reactions preceding this parent, recurse:
#     else:
#         # list of all reactions producing parent
#         tree_dict["reactions"] = []

#         for children_rxn in children_rxns:
#             # check if valid result, otherwise will stop walking this branch
#             if len(children_rxn) != 2:
#                 continue
#             # for each child reaction collect all precursors, side-products,
#             # then, each child will get own tree as well
#             # rxn_num = 'rxn1{}'.format(i)
#             # tree_dict[rxn_num]={}

#             reaction_dict = {}  # make dictionary for reaction
#             # later, add in reaction name?
#             reaction_dict["byproducts"] = []
#             reaction_dict["precursors"] = []

#             # get any side-products
#             child_prods = children_rxn[1].split(" ")  # list of child's products
#             for product in child_prods:
#                 if product != parent:
#                     reaction_dict["byproducts"].append(product)

#             # get precursors, for each precursor a dictionary:
#             child_precs = children_rxn[0].split(" ")
#             for prec in child_precs:
#                 if prec:
#                     prec_dict = make_tree_dict(
#                         parent=prec,
#                         pw_rxns=pw_rxns,
#                         parent_lvl=child_lvl,
#                         max_lvl=max_lvl,
#                     )
#                     reaction_dict["precursors"].append(prec_dict)

#             # append dictionary to the reactions list
#             tree_dict["reactions"].append(reaction_dict)
#         return tree_dict


# # ----------------------------------------------------------------------------
# # ------------------------------- combine all --------------------------------


# def get_pw_tree(
#     pw_leftrights: list[list[str, str]], max_lvl: int = 5
# ) -> list[list[str]]:
#     """creates a list where each index corresponds to the amount of reactions
#     it takes to get to the final product (i.e. first item is the final product
#     !!does not take into accound double left-side primaries"""
#     reaction_tree = []
#     # will get parent either in synthesis or reaction version
#     parent, parent_type = get_checked_ancestor(pw_leftrights)
#     if parent_type == "degradation":
#         pw_leftrights = inverse_rxns(pw_leftrights)

#     reaction_tree = make_tree_dict(parent, pw_leftrights, parent_lvl=0, max_lvl=max_lvl)
#     reaction_tree["type"] = parent_type
#     return reaction_tree


# # ----------------------------------------------------------------------------


# def get_reactions(entry: dict):  # obsolete?
#     """split=split_rxns(entry['REACTION-LAYOUT'],pre_post = True)
#     reactions= get_pw_tree(split)
#     return reactions"""
#     return None


# # ============================= getting reactions =============================


# def traverse(  # under construction
#     subdict: dict[int, str, dict],
#     reaction_list: list = [],
#     level_list: list = [],
# ) -> list[dict]:  # under constr.
#     reaction_list.append(subdict["compound"])
#     reactions = reaction_list  # ???? -> check later
#     for i in range(len(reaction_list)):  # ??? added without checking
#         if len(reactions) > 1:
#             reaction_list = [reaction_list.append(traverse(i)) for i in reactions]
#         elif len(reactions) == 1:
#             reaction_list.append(traverse(reactions[i]))

#     return


# def get_long_chain(pathway_tree: dict, max_level: int = 4) -> dict:
#     """for a given pathway_tree, traverses to get the reaction. gives first
#     resulting chain, often the only chain"""
#     # chain with minimum length of min_len
#     current_dict = pathway_tree
#     level = 0
#     chain = {level: pathway_tree["compound"]}
#     while level < max_level:
#         try:
#             current_dict = current_dict["reactions"][0]["precursors"][0]
#             level = current_dict["level"]
#             compound = current_dict["compound"]
#             if compound:
#                 chain[level] = compound
#         except:
#             print("stopped")
#             break
#     return chain


# def get_struct_rxn(
#     long_chain: dict,
#     compound_struct: dict,
#     second_dict: dict = {},
#     third_dict: dict = {},
# ) -> list[str]:
#     """for given chain of reactions, translates compounds to their inchies"""
#     structures = []
#     compounds = list(long_chain.values())
#     for i in compounds:
#         structure = ""
#         if i in compound_struct.keys():
#             structure = compound_struct[i]
#         else:
#             if second_dict and (i in second_dict) and second_dict[i]:
#                 structure = second_dict[i]
#             else:
#                 if third_dict and (i in third_dict) and third_dict[i]:
#                     structure = third_dict[i]
#         structures.append(structure)
#     return structures

# first, need to get


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


# def get_normalised_db(list_of_dict: list[dict[list]]) -> pd.DataFrame:
#     df = pd.DataFrame.from_records(list_of_dict)
#     for col in df.columns:
#         df = df.explode(col)
#     return df


# def split_columns(
#     df: pd.DataFrame, col: str, new_headers: list[str]
# ) -> pd.DataFrame:  # obsolete?
#     df_new = df
#     df_new[new_headers] = df_new[col].str.split("\(:", expand=True)
#     for header in new_headers:
#         df_new[header] = df_new[header].str.replace(header + " ", "")
#         df_new[header] = df_new[header].str.replace(")", "")
#         df_new[header] = df_new[header].str.replace("(", "")
#         df_new[header] = df_new[header].str.replace(" :", "")
#     # drop the original 'Name' column
#     df_new.drop(col, axis=1, inplace=True)
#     return df


# def lost_function(df, column):  # where did this come from?
#     """for unique_val in df[column].unique():
#     row = df.loc[df[column] == unique_val]"""
#     return None


# def colval_per_index(df: pd.DataFrame, colname: str = "REACTION-LAYOUT") -> list[str]:
#     "yields pw_reactions"
#     # unique_ind = get_indexes(df).tolist()
#     # for i in unique_ind:
#     #     colval = df[df.index.isin([i])][colname].tolist()
#     #     yield colval


def filter_by_rxn_len(df: pd.DataFrame, length: int = 4) -> pd.DataFrame:
    mask = np.array([len(x) > 3 for x in colval_per_index(df)])
    all_ind = np.array(df.index.unique())
    long_enough = all_ind[mask]
    new = df[df.index.isin(long_enough)]
    return new


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
# ================================ output ====================================


# obsolete
def pickle_mols(
    pathways, dict1, max_level=7, dict2={}, dict3={}, title="metacyc_reactions"
):
    return None


def write_reactions(
    pathways,
    dict1,
    max_level=7,
    dict2={},
    dict3={},
    title="metacyc_reactions",
    annot: dict = {},
    fulldf=None,
):
    """annotation from classes annotation, fulldf needed for if you want to
    obtain all taxonomy for pathway"""

    reactions
    all_pw_ids = [x for x in colval_per_index(pathways, colname="UNIQUE-ID")]
    assert len(all_pw_ids) == len(
        all_reaction_layouts
    ), "the number of ids and layouts does not match up"
    outfile_name = outfile_namer(title) + ".tsv"
    annotation_name = outfile_namer(title, "taxonomy") + ".tsv"
    with open(outfile_name, "w") as nf:
        nf.write("")
    if annot:
        with open(annotation_name, "w") as an:
            an.write("")

    all_pathway_trees = []
    for i in range(len(all_reaction_layouts)):
        pathway_tree = get_pw_tree(all_reaction_layouts[i], max_lvl=max_level)
        all_pathway_trees.append(pathway_tree)
        if not pathway_tree:
            continue
        long_chain = get_long_chain(pathway_tree, max_level=max_level)
        structures = get_struct_rxn(
            long_chain, compound_struct=dict1, second_dict=dict2, third_dict=dict3
        )

        with open(outfile_name, "a") as nf:
            nf.write("{}\t".format(all_pw_ids[i][0]))
            nf.write("\t".join(structures))
            nf.write("\n")

        # optional annotation file
        if annot:
            annotation = []
            if fulldf is None:
                continue
            taxonomic_range = (
                fulldf[fulldf["UNIQUE-ID"] == all_pw_ids[i][0]]["TAXONOMIC-RANGE"]
                .unique()
                .tolist()
            )
            for tax in taxonomic_range:
                if tax in annot.keys():
                    annotation.append(annot[tax])
                else:
                    annotation.append(tax)
            with open(annotation_name, "a") as an:
                an.write("{}\t".format(all_pw_ids[i][0]))
                an.write("\t".join(annotation))
                an.write("\n")

    print("finished writing to", outfile_name)
    return outfile_name


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

    # map compound ids to structures with meta_mols

    # get class annotation
    class_ann_loc = "/Users/lucina-may/thesis/metacyc/cleaner_classes.dat"
    annot = get_classes_annotation(class_ann_loc)

    # write_reactions(
    #     pathways_oi,
    #     main_dict,
    #     max_level=7,
    #     dict2=compound_inchi,
    #     dict3=compound_smiles,
    #     title="metacyc_reactions",
    #     annot=annot,
    #     fulldf=pathways,
    # )

    """#enough_split = split_rxns(reaction_layouts[0]) #UNIQUE-ID:PWY-6527
    for i in range(len(all_reaction_layouts)):
        print(i)
        pathway_tree = get_pw_tree(all_reaction_layouts[i], max_lvl=7)    
        long_chain = get_long_chain(pathway_tree, max_level=7)
        levels.append(len(long_chain)-1)
    fil = [x>=7 for x in levels]
    sum(fil)
    pathway_tree = get_pw_tree(all_reaction_layouts[537],max_lvl=150)
    long_chain = get_long_chain(pathway_tree,max_level=7)
    structs = get_struct_rxn(long_chain, main_dict,compound_inchi,compound_smiles)
    
    
    
    
    inchies = get_struct_rxn(long_chain, compound_inchi) #list with inchies in order
    smiles=get_struct_rxn(long_chain,compound_smiles)

    
    rxn_headers = ['name','LEFT-PRIMARIES','DIRECTION','RIGHT-PRIMARIES']
    split = split_rxn_columns(annot_pw,'REACTION-LAYOUT',rxn_headers )
    picklr(annot_pw)"""


if __name__ == "__main__":
    main()
