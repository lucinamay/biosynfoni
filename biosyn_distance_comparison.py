#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
||||||||||||  ㅇㅅㅇ  |||||||||||||||
____________________________________

title: inoutput                     ||
created: 2023.09.28 09:49           ||
author: lucina-may nollen           || 
institute: WUR Bioinformatics       ||
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
from inoutput import picklr 

#=========================== GLOBALS =========================
FP_FUNCTIONS={
    'biosynfoni':fp.binosynfoni_getter,
    'maccs':fp.maccs_getter,
    'morgan':fp.morgan_getter,
    'rdkit':fp.rdk_fp_getter,
    }

DIST_FUNCTIONS={
    'tanimoto':fp.tanimoto,
    'euclidean':fp.euclidean,
    }
#============================================================

    
def fp_to_distance(fp1, fp2, metric='tanimoto'):
    """for given molecule pair fp1 and fp2, gives distance"""
    #input handling & checking
    metric = metric.lower()
    assert metric in DIST_FUNCTIONS.keys(), \
        'unsupported metric: {}'.format(metric)
    assert fp1 and fp2,\
        ''
        
    distance=DIST_FUNCTIONS[metric]([fp1,fp2])
    return distance


from enum import Enum, auto

class FingerprintType(Enum):
    Morgan = "morgan"


def molpair_to_distance(
        molpair:list[Chem.Mol],
        fingerprint:str='morgan',
        #fingerprint:FingerprintType = FingerprintType.Morgan,
        metric:str='tanimoto'
        )->float:
    """uses default settings, check fingerprints.py for
    default settings per fingerprint type"""
    """if not isinstance(fingerprint, FingerprintType):
        raise TypeError(f"wrong signature for fingerprint: {type(fingerprint)} != FingerprintType" )
    """
    
    #fingerprint = fingerprint.value.lower()
    fingerprint = fingerprint.lower()
    assert fingerprint in FP_FUNCTIONS.keys(),\
        "unsupported fingerprint: {}".format(fingerprint)
    assert len(molpair)>1 and \
        isinstance(molpair[0],Chem.Mol) and\
        isinstance(molpair[1],Chem.Mol),\
        "please provide two Chem.Mol-type molecules"
    
    fp1 = FP_FUNCTIONS[fingerprint](molpair[0])
    fp2 = FP_FUNCTIONS[fingerprint](molpair[1])
    distance = fp_to_distance(fp1,fp2, metric=metric)
    return distance


def get_step_pairs(chain_mols: list, degree_of_sep: int = 1)->list:
    """for given list of molecules originating from one reaction,
    gets the available(non-null) molecule pairs separated by 
    degree of separation"""
    
    pairs = []
    for i in range(len(chain_mols)-degree_of_sep):
        if chain_mols[i] and chain_mols[i+degree_of_sep]:
            pairs.append([chain_mols[i],chain_mols[i+degree_of_sep]])
    return pairs


def yield_row(struct_loc: str):
    
    struct = pd.read_csv(struct_loc, sep = "\t", header = None)
    #get headers
    nums = ['{}'.format(i-1) for i in range(len(struct.columns))]
    nums.pop(0)
    headers = ['pathway']+nums
    struct.columns=headers
    #clean up
    struct.replace('', np.nan, inplace=True) # just in case one was missed
    filtered = struct.dropna(thresh=2) #drop completely empty ones
    final = filtered.replace(np.nan, '') # replace back for later fnxs
    #return pathways
    for i in range(len(final.index)):
        yield final.iloc[i].tolist()


def handle_dis_dict_kargs(
        metrics,
        fingerprints
        ) -> tuple[list[str]]:
    """handles non-list input, as well as giving the right options for 'all'"""
    
    #fix metrics
    if isinstance(metrics, str):
        if metrics == 'all':
            metrics = list(DIST_FUNCTIONS.keys())
        else:
            metrics = [metrics]
    elif isinstance(metrics, list) or isinstance(metrics, tuple):
        if 'all' in metrics:
            metrics = list(DIST_FUNCTIONS.keys())
    #fix fingerprints
    if isinstance(fingerprints,str):
        if fingerprints == 'all':
            fingerprints = list(FP_FUNCTIONS.keys())
        else:
            fingerprints = [fingerprints]
    elif isinstance(fingerprints, list) or isinstance(fingerprints, tuple):
        if 'all' in fingerprints:
            fingerprints = list(FP_FUNCTIONS.keys())
    #check validity        
    for metric in metrics:
        assert metric in DIST_FUNCTIONS.keys(),\
            "unsupported distance type: {}".format(metric)
    for fingerprint in fingerprints:
        assert fingerprint in FP_FUNCTIONS.keys(),\
            "unsupported fingerprint type: {}".format(fingerprint)
    return metrics, fingerprints 


def distances_dictify(
        molpair:list[Chem.Mol],
        distances_dict = {},
        metrics:list = ['all'],
        fingerprints:list = ['all']
        ) -> dict[dict[dict[list[float]]]]:
    
    metrics, fingerprints = handle_dis_dict_kargs(metrics, fingerprints)
    
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
            #if the mol is a Nonetype, distance will return None too
            #this is to not mess up further calculations (eg. with -1 or 0)
            if mol1 and mol2: 
                current_distance = molpair_to_distance(
                    mol1,
                    mol2, 
                    fingerprint=fingerprint, 
                    metric=metric
                    )
            distances_dict[metric][fingerprint].append(current_distance)
    return distances_dict


def yield_pairs(
        struct_loc:str,
        degree_of_sep:int=1
        ) -> list[list[list[str]]]:
    """returns as [pathway][pairnumber][0]product[1]precursor"""
    for i in yield_row(struct_loc):
        #remove the pathway annotation
        yield get_step_pairs(i[1:], degree_of_sep=degree_of_sep)


def mean_distances(chain_distances:list[float]) -> float:
    #remove the None values resulting from invalid Inchi's
    valid_list = [x for x in chain_distances if x] 
    return float(sum(valid_list))/float(len(valid_list))


def main():
    struct_loc = '/Users/lucina-may/thesis/scripts/0927_metacyc_reactions.tsv'
    degree_of_sep = 1
    pw_ids = [x[0] for x in yield_row(struct_loc)]
    pairs_per_pw =[x for x in yield_pairs(struct_loc, 
                                          degree_of_sep=degree_of_sep)]
    
    #getting a list of dictionaries per pathway, of the various distance values
    distances_per_pw = []
    for chain in pairs_per_pw:
        chain_dict = {}
        for pairs in chain:
            #later, add functionality to take the SMILES if it is not InChI
            molpairs = [Chem.MolFromInchi(x) for x in pairs]
            chain_dict=distances_dictify(molpairs, chain_dict)
        distances_per_pw.append(chain_dict)
    picklr(distances_per_pw, 'distances_per_pw_sep{}'.format(degree_of_sep))
    distances_per_pw['tanimoto']
        

if __name__ == "__main__":
    main()
