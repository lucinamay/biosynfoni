#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
||||||||||||  ㅇㅅㅇ  |||||||||||||||
____________________________________

title: biosyn distance comparison   ||
creaetd: 2023.09.28 09:49           ||
author: lucina-may nollen           || 
institute: WUR Bioinformatics       ||
____________________________________

||||||||||||  ()()()  |||||||||||||||

description: uses molecule pairs for comparison of distance estimation
between various fingerprints, including biosynfoni and binosynfoni(binary
biosynfoni)

"""

from rdkit import Chem
from rdkit.DataManip import Metric
from rdkit.Chem import AllChem
import pandas as pd
import numpy as np
import fingerprints as fp
from inoutput import picklr, jaropener, outfile_namer
from inoutput import entryfile_dictify as ann
import plotly
import plotly.express as px
from plotly.offline import download_plotlyjs, init_notebook_mode, iplot
import os
import matplotlib

#=========================== GLOBALS =========================
FP_FUNCTIONS={
    'biosynfoni':fp.biosynfoni_getter,
    'binosynfoni':fp.binosynfoni_getter,
    'maccs':fp.maccs_getter,
    'morgan':fp.morgan_getter,
    'rdkit':fp.rdk_fp_getter,
    }

SIM_FUNCTIONS={
    'c_tanimoto':fp.countanimoto,
    'tanimoto_dist':fp.tanimoto,
    'euclidean_dist':fp.euclidean,
    }

#for later:
from enum import Enum, auto

class FingerprintType(Enum):
    Morgan = "morgan"
    Daylight = "daylight"
"""if not isinstance(fingerprint, FingerprintType):
    raise TypeError(f"wrong signature for fingerprint: {type(fingerprint)} != FingerprintType" )
"""
#fingerprint.values.lower()


#=============================================================================

#=============================== distance fnx ================================
    
def fp_to_distance(fp1, fp2, metric='c_tanimoto'):
    """for given molecule pair fp1 and fp2, gives distance"""
    #input handling & checking
    metric = metric.lower()
    assert metric in SIM_FUNCTIONS.keys(), \
        'unsupported metric: {}'.format(metric)
    assert fp1 and fp2,\
        'please provide two fingerprints'
    distance = None #init
    distance=SIM_FUNCTIONS[metric]([fp1,fp2])
    return distance

def _twomols_to_distance(
        mol1:Chem.Mol,
        mol2:Chem.Mol,
        fingerprint:str = 'morgan',
        #fingerprint:FingerprintType = FingerprintType.Morgan,
        metric:str = 'c_tanimoto'
        ) -> float:
    fingerprint = fingerprint.lower()
    assert fingerprint in FP_FUNCTIONS.keys(),\
        "unsupported fingerprint: {}".format(fingerprint)
    assert isinstance(mol1,Chem.Mol) and isinstance(mol2,Chem.Mol),\
        "please provide two Chem.Mol-type molecules"
    fp1 = FP_FUNCTIONS[fingerprint](mol1)
    fp2 = FP_FUNCTIONS[fingerprint](mol2)
    distance = fp_to_distance(fp1,fp2, metric=metric)
    return distance
            

def molpair_to_distance(
        molpair:list[Chem.Mol],
        fingerprint:str = 'morgan',
        #fingerprint:FingerprintType = FingerprintType.Morgan,
        metric:str = 'c_tanimoto_similarity'
        ) -> float:
    """uses default settings, check fingerprints.py for
    default settings per fingerprint type"""
    
    assert len(molpair)>1, "please provide two molecules"
    distance = _twomols_to_distance(
        molpair[0],
        molpair[1],
        fingerprint=fingerprint,
        metric=metric)
    return distance

def inchi_pair_to_distance(
        pair:list[str],
        metric:str = 'c_tanimoto',
        fingerprint:str = 'binosynfoni',
        ) -> float:
    distance = None
    assert len(pair) > 1, "only one molecule given"
    
    mol1 = Chem.MolFromInchi(pair[0])
    mol2 = Chem.MolFromInchi(pair[1])
    if isinstance(mol1, Chem.Mol) and isinstance(mol2, Chem.Mol):
        distance = _twomols_to_distance(
            mol1,
            mol2,
            metric=metric,
            fingerprint=fingerprint
            )
    return distance

#=============================================================================

#=============================== getting pairs ===============================

def get_step_pairs(chain_mols: list, degree_of_sep: int = 1)->list:
    """for given list of molecules originating from one reaction,
    gets the available(non-null) molecule pairs separated by 
    degree of separation"""
    
    pairs = []
    for i in range(len(chain_mols)-degree_of_sep):
        if chain_mols[i] and chain_mols[i+degree_of_sep]:
            pairs.append([chain_mols[i],chain_mols[i+degree_of_sep]])
    return pairs

def yield_row_from_file(struct_loc: str):
    """yields rows from csv, cleaning up the rows in the process"""
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

def yield_pairs_from_file(
        struct_loc:str,
        degree_of_sep:int = 1
        ) -> list[list[list[str]]]:
    """returns as [pathway][pairnumber][0]product[1]precursor"""
    for i in yield_row_from_file(struct_loc):
        #remove the pathway annotation
        yield get_step_pairs(i[1:], degree_of_sep=degree_of_sep)

def dictify_pw_pairs(
        struct_loc:str, 
        degree_of_sep:int = 1
        ) -> dict[list[list[str]]]:
    """returns {'pathway-id':[[product,precursor],[precursor,preprecursor]]}"""
    dictionary = {}
    for row in yield_row_from_file(struct_loc):
        pathway_id = row[0]
        pairs = get_step_pairs(row[1:], degree_of_sep=degree_of_sep)
        dictionary[pathway_id] = pairs
    return dictionary

def _label_pairs(
        distances_per_pathway:dict[list[float]]
        ) -> tuple[str,float]:
    """returns list of (pathway,distances)"""
    labeled_distances = []
    for pw_id in distances_per_pathway.keys():
        for distance in distances_per_pathway[pw_id]:
            labeled_distances.append((pw_id,distance))
    return labeled_distances

def get_pairs_dist_df(
        labeled_distances:list[tuple[str,float]],
        metric = 'c_tanimoto_similarity'
        ) -> pd.DataFrame:
    df = pd.DataFrame(labeled_distances)
    df.columns = ['pathway', 'inchis_sep']
    #find way to convert to mol only once instead of continuously
    for fp_type in FP_FUNCTIONS.keys():
        df[fp_type] = df['inchis_sep'].apply(
            lambda x: inchi_pair_to_distance(
                x,
                fingerprint=fp_type,
                metric='c_tanimoto'))
    return df

def annotate_pathways(df, pw_tax_file:str, tax_text_file:str) -> pd.DataFrame:
    ndf = df.copy()
    start_val_sep = ' - '
    entry_sep="//"
    encoding='iso8859-1'
    pw_tax = ann(
        pw_tax_file,
        keyvals=('UNIQUE-ID','TAXONOMIC-RANGE'),
        start_val_sep=start_val_sep,
        encoding=encoding,
        entry_sep=entry_sep)
    
    class_text = ann(
        tax_text_file,
        keyvals=('UNIQUE-ID','COMMON-NAME'),
        start_val_sep=start_val_sep,
        encoding=encoding,
        entry_sep=entry_sep)
    print(len(pw_tax),len(class_text))
    ndf['tax_id']=ndf['pathway'].replace(to_replace=pw_tax)
    ndf['taxonomic_range']=ndf['tax_id'].replace(to_replace=class_text)
    ndf['pathway_name']=ndf['pathway'].replace(to_replace=class_text)
    #ndf['pathway_name']=ndf['pathway_name'].fillna(ndf['pathway'])
    return ndf

def plot_two_cols(
        df,
        col1,
        col2,
        colour_label='taxonomic_range',
        shape_label='',
        hover_data=['pathway_name','taxonomic_range']):
    #fig = px.scatter(df,x=col1, y=col2, hover_data=['class1,class2'])
    fig = px.scatter(
        df,
        x=col1, 
        y=col2,
        color=colour_label,
        hover_data=hover_data)
    fig.update_layout(
        title='Tanimoto similarity for biosynthetic (x,x+1) pairs'
        )
    outfilename = outfile_namer(f'{col1}_{col2}',colour_label)
    fig.write_html(f'{outfilename}.html', auto_open=True)
    return None

def main():
    struct_loc = '../scripts/0927_metacyc_reactions.tsv'
    degree_of_sep = 1
    
    """#add function for folder moving

    #pickle_loc = f"{outfile_namer(f'distances_per_pw_sep{degree_of_sep}')}.pickle"
    pickle_loc = '../scripts/0929_distances_per_pw.pickle'
    
    
    
    if not os.path.exists(pickle_loc):
        distances_per_pw =_run_distancer(struct_loc, degree_of_sep)
    else:
        distances_per_pw = jaropener(pickle_loc)
    #test = [list(x['tanimoto']) for x in distances_per_pw]
    """
    
if __name__ == "__main__":
    main()