#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
||||||||||||  ㅇㅅㅇ  |||||||||||||||
____________________________________

title: inoutput                     ||
creaetd: 2023.09.28 09:49           ||
language: python                    ||
author: lucina-may nollen           || 
institute: WUR Bioinformatics       ||
student no: 1197142                 ||
____________________________________

||||||||||||  ()()()  |||||||||||||||

description:    general functions needed in the other file

- readr(filename:list)
- entry_parser(lines:list, sep:str = "$$$$")

"""

import rdkit
from rdkit import Chem
from rdkit.DataManip import Metric
from rdkit.Chem import AllChem
import pandas as pd
import numpy as np
import fingerprints as fp

#=========================== GLOBALS =========================
FP_FUNCTIONS={
    'biosynfoni':fp.binosynfoni_getter,
    'maccs':fp.maccs_getter,
    'morgan':fp.morgan_getter,
    'rdkit':fp.rdk_fp_getter,
    }

#============================================================


def tanimoto(expl_bitvectors:list)->float:
    array = Metric.rdMetricMatrixCalc.GetTanimotoDistMat(\
        expl_bitvectors)
    return array.tolist()[0]

def euclidean(expl_bitvectors:list):
    array =  Metric.rdMetricMatrixCalc.GetEuclideanDistMat(\
        expl_bitvectors)
    return array.tolist()[0]
    
DIST_FUNCTIONS={
    'tanimoto':tanimoto,
    'euclidean':euclidean,
    }

def fp_to_distance(fp1,fp2, metric='tanimoto'):
    """for given molecule pair fp1 and fp2, gives distance"""
    #input handling & checking
    metric = metric.lower()
    assert metric in DIST_FUNCTIONS.keys(),\
        'unsupported metric: {}'.format(metric)
        
    distance=DIST_FUNCTIONS[metric]([fp1,fp2])
    return distance

def molpair_to_distance(mol1,mol2, 
                        fingerprint:str='morgan',
                        metric:str='tanimoto')->float:
    """uses default settings, check fingerprints.py for
    default settings per fingerprint type"""
    
    fingerprint = fingerprint.lower()
    assert fingerprint in FP_FUNCTIONS.keys(),\
        "unsupported fingerprint: {}".format(fingerprint)
    
    fp1 = FP_FUNCTIONS[fingerprint](mol1)
    fp2 = FP_FUNCTIONS[fingerprint](mol2)
    distance = fp_to_distance(fp1,fp2, metric=metric)
    return distance

def get_step_pairs(chain_mols:list, degree_of_sep:int=1):
    """for given list of molecules originating from one reaction,
    gets the available(non-null) molecule pairs separated by 
    degree of separation"""
    pairs = []
    for i in range(len(chain_mols)-degree_of_sep):
        if chain_mols[i] and chain_mols[i+degree_of_sep]:
            pairs.append([chain_mols[i],chain_mols[i+degree_of_sep]])
    return pairs

def yield_row(struct_loc:str):
    struct = pd.read_csv(struct_loc,sep="\t", header=None)
    #get headers
    nums = ['{}'.format(i-1) for i in range(len(struct.columns))]
    nums.pop(0)
    headers = ['pathway']+nums
    struct.columns=headers
    #clean up
    struct.replace('', np.nan, inplace=True) # just in case one was missed
    filtered = struct.dropna(thresh=2) #drop completely empty ones
    filtered.replace(np.nan, '', inplace=True) # replace back for later fnxs
    #return pathways
    for i in range(len(filtered.index)):
        yield filtered.iloc[i].tolist()

def distances_dictify(molpair, distances_dict={}):
    if not distances_dict:
        #init.
        for metric in DIST_FUNCTIONS.keys():
            distances_dict[metric]={}
            for fingerprint in FP_FUNCTIONS.keys():
                distances_dict[metric][fingerprint]=[]
    
    for metric in distances_dict.keys():
        for fingerprint in distances_dict[metric].keys():
            mol1,mol2 = molpair[0],molpair[1]
            current_distance = None
            if mol1 and mol2:
                current_distance = molpair_to_distance(
                    mol1,
                    mol2, 
                    fingerprint=fingerprint, 
                    metric=metric
                    )
            distances_dict[metric][fingerprint].append(current_distance)
    return distances_dict


def main():
    struct_loc = '/Users/lucina-may/thesis/scripts/0927_metacyc_reactions.tsv'
    reactions = [x[1:] for x in yield_row(struct_loc)]
    pw_ids = [x[0] for x in yield_row(struct_loc)]
    pairs_per_chain=[]
    for chain in reactions[:10]:
        pairs_per_chain.append(get_step_pairs(chain))
    
    distances_per_chain = []
    for chain in pairs_per_chain:
        chain_dict = {}
        for pairs in chain:
            molpairs = [Chem.MolFromInchi(x) for x in pairs]
            chain_dict=distances_dictify(molpairs, chain_dict)
        distances_per_chain.append(chain_dict)
        

if __name__ == "__main__":
    main()