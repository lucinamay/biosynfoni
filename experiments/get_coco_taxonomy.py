#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 11:17:39 2023

@author: lucina-may
"""
import os
import sys
import numpy as np
import pandas as pd
from rdkit import Chem

sys.path.append(os.path.abspath(os.path.join(sys.path[0], os.pardir, "src")))
# for intra-biosynfoni-code running
sys.path.append(
    os.path.abspath(os.path.join(sys.path[0], os.pardir, "src", "biosynfoni"))
)
from biosynfoni.inoutput import jaropener
from biosynfoni.subkeys.get_version import get_biosynfoni
from tsne import tsner
from utils import figuremaking as fm

# from inoutput import jaropener
# from concerto_fp import get_biosynfoni
# from tsne import tsner
# import figuremaking as fm

cwd = os.getcwd()
taxonomy_files = "../wikidata_compounds_taxonomy/kingdom_annotated_mols"
os.chdir(taxonomy_files)
file_list = os.listdir(".")

all_taxonomy = pd.DataFrame()
for picklefile in file_list:
    a = jaropener(picklefile)
    all_taxonomy = pd.concat([all_taxonomy, a])

os.chdir(cwd)

all_taxonomy["molecular_structure_from_inchi"] = (
    all_taxonomy["inchi"].str.split("/").str[1:3]
)
all_taxonomy["shortinchi"] = all_taxonomy["molecular_structure_from_inchi"].apply(
    lambda x: "/".join(x)
)
all_taxonomy.drop("molecular_structure_from_inchi", axis=1)
taxonomy_as_list = all_taxonomy.groupby(all_taxonomy["structure"]).aggregate(list)

coco_info = ["CNP-ID", "InChI", "smiles", "molecular_formula", "NPLS"]
coco = pd.read_csv(
    "/Users/lucina-may/thesis/arch/coconut_full_info(incl_inchismiles)/output/0907_COCONUT_DB_info.tsv",
    sep="\t",
    header=None,
)

coco.columns = coco_info
coco["molecular_structure_from_inchi"] = coco["InChI"].str.split("/").str[1:3]
coco["shortinchi"] = coco["molecular_structure_from_inchi"].apply(lambda x: "/".join(x))

taxonomied_inchis = all_taxonomy["shortinchi"].tolist()
# taxonomied_mol_struct = ["/".join(x.split("/")[1:3]) for x in taxonomied_inchis]
# coco_inchi_short = [x.split('/')[1:3] for x in taxonomied_inchis]]
# coco["split"] = coco["InChI"].str.split("/")
# coco["short2"] = coco["InChI"].str.split("/")[2]

# coco["inchi_short"] = "/".join([coco["split"].str[1], coco["split"].Series[2]])

# coco["taxonomy_present"] = coco["shortinchi"].isin(taxonomied_inchis)
# coco_taxonomy = coco[coco["taxonomy_present"] == True]
# coco_taxonomy.drop('taxonomy_present', axis=1)

# get dataframe with available taxonomy labeling:


coco_taxonomy = pd.merge(coco, all_taxonomy, how="inner", on="shortinchi")
# get biosynfoni for all inchis (using the inchis from the coconut database)
coco_taxonomy["biosynfoni_3"] = coco_taxonomy["InChI"].apply(
    lambda x: get_biosynfoni(Chem.MolFromInchi(x), "fps_full_3")
)
array = np.array(coco_taxonomy["biosynfoni_3"].to_list())
tsne_df = tsner(array)
components = tsne_df.columns.to_list()
tsne_df["kingdom"] = coco_taxonomy["kingdom"]
fm.plot_two_cols(
    tsne_df,
    components[0],
    components[1],
    figtitle="BioSynFoNi t-SNE per taxonomic kingdom",
    colour_label="kingdom",
    colour_dict=fm.TAXONOMY_COLOURS,
)
