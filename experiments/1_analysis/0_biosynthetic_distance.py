#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 16:15:34 2023

@author: lucina-may
"""
import os, logging, argparse

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


def compound_graph(reaction_info: pd.DataFrame, with_compounds: list) -> pd.DataFrame:
    """Returns graph of reagents per MetaCyc pathway ID"""

    compound_graph = nx.DiGraph()
    for _, row in reaction_info.iterrows():
        if not with_compounds or (
            row["left"] in with_compounds and row["right"] in with_compounds
        ):
            compound_graph.add_edge(
                row["left"], row["right"], reaction_id=row["reaction_id"]
            )
    return compound_graph


def longest_path(pathway_graph):
    if not nx.is_directed_acyclic_graph(pathway_graph):
        longest_path = []
        for node in pathway_graph.nodes:
            for path in nx.all_simple_paths(
                pathway_graph, source=node, target=pathway_graph.nodes
            ):
                if len(path) > len(longest_path):
                    longest_path = path
        return longest_path
    else:
        return nx.dag_longest_path(pathway_graph)


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


def prepare_reaction_chains(
    pathways: pd.DataFrame, compounds: pd.DataFrame
) -> pd.DataFrame:
    available_compounds = compounds["unique_id"].tolist()

    compound_graphs = {
        pathway_id: compound_graph(reaction_info, available_compounds)
        for pathway_id, reaction_info in pathways.groupby("pathway_id")
    }

    longest_paths = {
        pathway_id: longest_path(compound_graph)
        for pathway_id, compound_graph in compound_graphs.items()
    }

    longest_paths = pd.DataFrame.from_dict(longest_paths, orient="index")
    longest_paths = longest_paths.explode(0)  # explode to get one row per compound
    longest_paths = compound_id_to_structure(
        longest_paths,
        dict1=compounds["inchi"].to_dict(),
        dict2=compounds["smiles"].to_dict(),
        dict3=compounds["non_standard_inchi"].to_dict(),
    )
    return longest_paths


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

    pathways = pd.read_csv(
        "metacyc_pathways.tsv", sep="\t", index_col="pathway_id", header=True
    )
    meta_mols = pd.read_csv(
        "meta_compounds.tsv", sep="\t", index_col="unique_id", header=True
    )

    if args.filter:
        common_compound_path = os.path.join(temp_dir, "coco_meta.tsv")
        if not os.path.exists(common_compound_path):
            coco_mols_path = os.path.join(temp_dir, "coconut_compounds.tsv")
            if not os.path.exists(coco_mols_path):
                coco_mols = sdf_to_compounds(sdf_path=args.filter)  # standard query
                coco_mols.to_csv(coco_mols_path, sep="\t", index=True)

            coco_mols = pd.read_csv(coco_mols_path, sep="\t", index_col="coconut_id")
            compounds = get_common_compounds(coco_mols, meta_mols, inchi_parts=3)
            compounds.to_csv(common_compound_path, sep="\t", index=False)
        compounds = pd.read_csv(common_compound_path, sep="\t")

        pathways_oi = get_pathways_of_interest(pathways, compounds)
    else:
        pathways_oi = pathways

    pw_chains = prepare_reaction_chains(pathways_oi, meta_mols["unique_id"].to_list())

    # go to parent directory
    os.chdir(output_dir)
    if args.smiles or args.filter:
        mol_dict = pd.Series(
            meta_mols["smiles"].values, index=meta_mols.index
        ).to_dict()
        pw_chains_structs = pw_chains.replace(mol_dict)
        pw_chains_structs.to_csv(f"{args.output_path}_smiles.tsv", sep="\t", index=True)

    pw_chains.to_csv(f'{args.output_path}.tsv', sep="\t", index=True)

    os.chdir(iwd)
    exit(0)


if __name__ == "__main__":
    main()
