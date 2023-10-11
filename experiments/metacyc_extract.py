#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 16:15:34 2023

@author: lucina-may
"""
import sys
import subprocess
from sys import argv

import numpy as np
import pandas as pd

sys.path.append("../src/")
from biosynfoni.inoutput import *

# ============================== input handling ===============================


def get_df(
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
        extractr,
        starts=info2extract,
        remove_starts=remove_starts,
        strip_extra=strip_extra,
    )

    # get data frame, will have doubles if rxnslists (because of 'expanding')
    df = get_normalised_db(all_vals)
    return df


def get_relevant_pws(
    df: pd.DataFrame, reactions_of_interest: list, min_pw_len: int = 4
) -> pd.DataFrame:
    """from DataFrame with pathways, for each pathway with at least min_pw_len,
    returns the Dataframe containing them and their pathway unique-id's"""
    # get pathways containing reactions----------------------------------------
    # pathway-df row containing rxn of interest
    pw_rxns_oi = df_filter(df, "REACTION-LAYOUT", reactions_of_interest)
    # filter out only those with index of interest
    ind_oi = get_indexes(pw_rxns_oi).tolist()  # index of interests
    pw_oi = df[df.index.isin(ind_oi)]

    # check if same indexes
    check_ind = get_indexes(pw_oi).tolist()
    assert check_ind == ind_oi, "error with COCONUT-df filtering"

    # remove taxonomic range and species for the rxn annotation + undouble:
    columns_to_rem = pw_oi.columns.values.tolist()
    columns_to_rem.remove("UNIQUE-ID")
    columns_to_rem.remove("REACTION-LAYOUT")
    clean_pw_oi = remove_cols(pw_oi, columns_to_rem)

    # get df with pathways that have at least four reactions in them
    long_enough_pws = filter_by_rxn_len(clean_pw_oi, length=min_pw_len)
    # all_reaction_layouts = [x for x in colval_per_index(enough_pw_oi)]
    return long_enough_pws


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


def id_repr(
    filename: str = "../metacyc/compound-links.dat", representation: str = "inchi"
) -> dict:
    reprs = ["inchi", "smiles"]
    assert representation in reprs, "representation type unavailable"
    column = reprs.index(representation) + 1
    id_repr = {}
    with open(filename) as comp:
        for line in comp:
            if line:
                row = line.strip().split("\t")
                if len(row) > 1:
                    if row[column] and row[0]:
                        try:
                            id_repr[row[0].upper()] = row[column]
                        except:
                            continue
    cleaned = clean_dict(id_repr)
    return cleaned


# fix capitalisation


# ----------------------------------------------------------------------------
# ============================================================================


# ============================= traverse pathway =============================
# ------------------------- handle the reaction-layout------------------------
def split_rxns(reaction_layouts: list[str], pre_post: bool = True) -> list[list[str]]:
    """splits reaction terms into left, direction, right. if pre_post == True
    the reaction is automatically written as [precursor,product]"""

    rxns = []
    for item in reaction_layouts:
        left = " ".join(item.split("(:")[1].strip(") ").split(" ")[1:])
        drxn = item.split("(:")[2].strip("\) ").split(":")[-1]
        right = " ".join(item.split("(:")[3].strip(") ").split(" ")[1:])
        if not (left and drxn and right):
            continue
        if left == right:
            continue
        if not pre_post:
            rxns.append([left, drxn, right])
        elif pre_post:
            if "R2L" in drxn or "2L" in drxn or "R2" in drxn:
                pre = right
                post = left
            else:
                pre, post = left, right
            rxns.append([pre, post])
    return rxns


def inverse_rxns(pre_posts: list[list[str]]) -> list[list[str]]:
    """mainly for degradation pathways, where the 'parent' is actually
    the terminal precursor in the pathway"""
    inversed = []
    for pair in pre_posts:
        inversed.append([pair[1], pair[0]])
    return inversed


# ----------------------------------------------------------------------------
# ------------------------------ find start point ----------------------------


def find_ancestor_rxn(split_rxns: list[list[str]]) -> list[list[str]]:
    """with a list [pre, post]"""
    # find oldest ancestor (parent) reaction
    parents = []
    all_pre = [item[0] for item in split_rxns]
    all_post = [item[1] for item in split_rxns]

    split_pre = []
    for pres in all_pre:
        splitted = pres.split(" ")
        if isinstance(splitted, list):
            split_pre = split_pre + splitted
        elif isinstance(splitted, str):
            split_pre.append(splitted)

    for i in range(len(all_post)):
        if len(all_post[i].split(" ")) > 1:
            continue
        elif all_post[i] in split_pre:
            continue
        else:
            parents.append(split_rxns[i])
    return parents


def check_ancestor_res(ancestor_rxns_res: list[list[str]]) -> str:
    parent = ""

    if ancestor_rxns_res:
        # select the reaction from list
        parent_rxn = ancestor_rxns_res[0]
        assert len(parent_rxn) == 2, "error in the parent rxn result"
        # check if the result reaction exists and
        if parent_rxn[1]:
            # print(parent_rxn[1])
            if len(parent_rxn[1].split(" ")) == 1:
                parent = parent_rxn[1]  # selected product
    return parent


def get_checked_ancestor(prec_prods: list[list[str]]) -> tuple[str, str]:
    """gets parent compound of pathway, and returns if its a synthesis
    or a degradation type pathway"""
    parent_type = "synthesis"

    synthesis_par_rxns = find_ancestor_rxn(prec_prods)
    # check if there is a result and if there is only one result rxn
    parent = check_ancestor_res(synthesis_par_rxns)
    if not parent:
        parent_type = "degradation"
        # try to find the degradation parent
        degradation_par_rxns = find_ancestor_rxn(inverse_rxns(prec_prods))
        parent = check_ancestor_res(degradation_par_rxns)
        if not parent:
            parent_type = "first"
            print(
                "Could not find parent",
                "defaulted to first product of first rxn:",
                prec_prods[0][1],
            )
            try:
                parent = prec_prods[0][1].split(" ")[0]
            except:
                print("Failed getting 1st product of 1st rxn, return empty")
                parent = ""
                parent_type = ""
    return parent, parent_type


# ----------------------------------------------------------------------------
# --------------------------------- tree-making ------------------------------


def get_child_rxn(parent: list[str], rxn_layout: list[list[str]]):  # obsolete?
    """par_child = [parent]
    child_found = False
    for child_i in range(len(rxn_layout)):
        print(rxn_layout[child_i][1], parent[0])
        if rxn_layout[child_i][1] == parent[0]:
            child_found = True
            print('True')
            grandchildren = get_child_rxn(rxn_layout[child_i],rxn_layout)
            par_child.append(rxn_layout[child_i],[grandchildren])
    if not child_found:
        return [parent]
    return par_child"""
    return None


def getchild_rxns(parent: str, prec_prods: list[list[str]]) -> list[list[str]]:
    """for parent compound, returns list of the children reactions"""
    children = []
    for reaction in range(len(prec_prods)):
        # if the parent appears in the reaction's products
        if parent in prec_prods[reaction][1].split(" "):
            # append the reactions [precursors, products]
            children.append(prec_prods[reaction])
    return children


def make_tree_dict(
    parent: str, pw_rxns: list[list[str]], parent_lvl: int = 0, max_lvl: int = 5
) -> dict[int, str, dict]:
    tree_dict = {}  # total dictionary
    tree_dict["level"] = parent_lvl
    tree_dict["compound"] = parent

    # get children of parent
    children_rxns = getchild_rxns(parent, pw_rxns)
    child_lvl = parent_lvl + 1

    # check level to prevent overrecursion
    if child_lvl > max_lvl:
        # print('stop at lvl{}: {}'.format(max_lvl,parent))
        return tree_dict

    # if no reactions leading to parent's compound, return the tree dict
    if not children_rxns:
        return tree_dict

    # if reactions preceding this parent, recurse:
    else:
        # list of all reactions producing parent
        tree_dict["reactions"] = []

        for children_rxn in children_rxns:
            # check if valid result, otherwise will stop walking this branch
            if len(children_rxn) != 2:
                continue
            # for each child reaction collect all precursors, side-products,
            # then, each child will get own tree as well
            # rxn_num = 'rxn1{}'.format(i)
            # tree_dict[rxn_num]={}

            reaction_dict = {}  # make dictionary for reaction
            # later, add in reaction name?
            reaction_dict["byproducts"] = []
            reaction_dict["precursors"] = []

            # get any side-products
            child_prods = children_rxn[1].split(" ")  # list of child's products
            for product in child_prods:
                if product != parent:
                    reaction_dict["byproducts"].append(product)

            # get precursors, for each precursor a dictionary:
            child_precs = children_rxn[0].split(" ")
            for prec in child_precs:
                if prec:
                    prec_dict = make_tree_dict(
                        parent=prec,
                        pw_rxns=pw_rxns,
                        parent_lvl=child_lvl,
                        max_lvl=max_lvl,
                    )
                    reaction_dict["precursors"].append(prec_dict)

            # append dictionary to the reactions list
            tree_dict["reactions"].append(reaction_dict)
        return tree_dict


# ----------------------------------------------------------------------------
# ------------------------------- combine all --------------------------------


def get_pw_tree(reaction_layouts: list[str], max_lvl: int = 5) -> list[list[str]]:
    """creates a list where each index corresponds to the amount of reactions
    it takes to get to the final product (i.e. first item is the final product
    !!does not take into accound double left-side primaries"""
    prec_prods = split_rxns(reaction_layouts)
    if not prec_prods:
        return []
    # will get parent either in synthesis or reaction version
    parent, parent_type = get_checked_ancestor(prec_prods)
    if parent_type == "degradation":
        prec_prods = inverse_rxns(prec_prods)

    reaction_tree = make_tree_dict(parent, prec_prods, parent_lvl=0, max_lvl=max_lvl)
    reaction_tree["type"] = parent_type
    return reaction_tree


# ----------------------------------------------------------------------------


def get_reactions(entry: dict):  # obsolete?
    """split=split_rxns(entry['REACTION-LAYOUT'],pre_post = True)
    reactions= get_pw_tree(split)
    return reactions"""
    return None


# ============================= getting reactions =============================


def traverse(  # under construction
    subdict: dict[int, str, dict],
    reaction_list: list = [],
    level_list: list = [],
) -> list[dict]:  # under constr.
    reaction_list.append(subdict["compound"])
    reactions = reaction_list  # ???? -> check later
    for i in range(len(reaction_list)):  # ??? added without checking
        if len(reactions) > 1:
            reaction_list = [reaction_list.append(traverse(i)) for i in reactions]
        elif len(reactions) == 1:
            reaction_list.append(traverse(reactions[i]))

    return


def get_long_chain(pathway_tree: dict, max_level: int = 4) -> dict:
    """for a given pathway_tree, traverses to get the reaction. gives first
    resulting chain, often the only chain"""
    # chain with minimum length of min_len
    current_dict = pathway_tree
    level = 0
    chain = {level: pathway_tree["compound"]}
    while level < max_level:
        try:
            current_dict = current_dict["reactions"][0]["precursors"][0]
            level = current_dict["level"]
            compound = current_dict["compound"]
            if compound:
                chain[level] = compound
        except:
            print("stopped")
            break
    return chain


def get_struct_rxn(
    long_chain: dict,
    compound_struct: dict,
    second_dict: dict = {},
    third_dict: dict = {},
) -> list[str]:
    """for given chain of reactions, translates compounds to their inchies"""
    structures = []
    compounds = list(long_chain.values())
    for i in compounds:
        structure = ""
        if i in compound_struct.keys():
            structure = compound_struct[i]
        else:
            if second_dict and (i in second_dict) and second_dict[i]:
                structure = second_dict[i]
            else:
                if third_dict and (i in third_dict) and third_dict[i]:
                    structure = third_dict[i]
        structures.append(structure)
    return structures


def compound_id_to_structure():  # under construction
    return None


# ============================ df transformations =============================


def df_filter(df: pd.DataFrame, col: str, lst_yes: list) -> pd.DataFrame:
    mask = df[col].isin(lst_yes)
    df_of_interest = df[mask]
    return df_of_interest


def get_indexes(df: pd.DataFrame) -> list[int]:
    unique_index = df.index.unique()
    return unique_index


def get_unique_vals(df: pd.DataFrame, col: str) -> list:
    unique_col_vals = df.index.get_level_values(col).unique()
    return unique_col_vals


def remove_cols(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    new = df
    for i in cols:
        new = new.drop(columns=i).drop_duplicates()
    return new


def annotate(df: pd.DataFrame, replace_dict: dict) -> pd.DataFrame:
    return df.replace(to_replace=replace_dict)


def get_normalised_db(list_of_dict: list[dict[list]]) -> pd.DataFrame:
    df = pd.DataFrame.from_records(list_of_dict)
    for col in df.columns:
        df = df.explode(col)
    return df


def split_columns(
    df: pd.DataFrame, col: str, new_headers: list[str]
) -> pd.DataFrame:  # obsolete?
    df_new = df
    df_new[new_headers] = df_new[col].str.split("\(:", expand=True)
    for header in new_headers:
        df_new[header] = df_new[header].str.replace(header + " ", "")
        df_new[header] = df_new[header].str.replace(")", "")
        df_new[header] = df_new[header].str.replace("(", "")
        df_new[header] = df_new[header].str.replace(" :", "")
    # drop the original 'Name' column
    df_new.drop(col, axis=1, inplace=True)
    return df


def lost_function(df, column):  # where did this come from?
    """for unique_val in df[column].unique():
    row = df.loc[df[column] == unique_val]"""
    return None


def colval_per_index(df: pd.DataFrame, colname: str = "REACTION-LAYOUT") -> list[str]:
    "yields pw_reactions"
    unique_ind = get_indexes(df).tolist()
    for i in unique_ind:
        colval = df[df.index.isin([i])][colname].tolist()
        yield colval


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
    all_reaction_layouts = [x for x in colval_per_index(pathways)]
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
    # pathway_loc = argv[1]
    pathway_loc = "/Users/lucina-may/thesis/metacyc/pathways.dat"
    info2extract = [
        "UNIQUE-ID",
        "REACTION-LIST",
        "SPECIES",
        "TAXONOMIC-RANGE",
        "REACTION-LAYOUT",
    ]
    pathways = get_df(pathway_loc, info2extract)

    # get reactions containing COCONUT
    pathways_ofinterest_loc = "/Users/lucina-may/thesis/metacyc/interest_pws.tmp"
    reactions_of_interest = extractr(
        readr(pathways_ofinterest_loc), ["REACTION-LAYOUT"], strip_extra=" - "
    )["REACTION-LAYOUT"]
    # get pathways containing COCONUT compounds of defined minimal length
    pathways_oi = get_relevant_pws(pathways, reactions_of_interest, min_pw_len=4)
    all_reaction_layouts = [x for x in colval_per_index(pathways_oi)]
    pw_num = [x for x in colval_per_index(pathways_oi, "UNIQUE-ID")]

    # get annotation of compounds
    compounds_loc = "/Users/lucina-may/thesis/metacyc/cleaner_compounds.dat"
    second_db = "/Users/lucina-may/thesis/metacyc/compound-links.dat"

    compoundinfo = ["UNIQUE-ID", "SMILES", "INCHI", "NON-STANDARD-INCHI"]
    compounds = mergeinchis(get_df(compounds_loc, compoundinfo))
    main_dict = to_conversion_dict(compounds)
    compound_inchi = id_repr(filename=second_db)  # backup dictionaries
    compound_smiles = id_repr(filename=second_db, representation="smiles")

    # get class annotation
    class_ann_loc = "/Users/lucina-may/thesis/metacyc/cleaner_classes.dat"
    annot = get_classes_annotation(class_ann_loc)

    write_reactions(
        pathways_oi,
        main_dict,
        max_level=7,
        dict2=compound_inchi,
        dict3=compound_smiles,
        title="metacyc_reactions",
        annot=annot,
        fulldf=pathways,
    )

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
