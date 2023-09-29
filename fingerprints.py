#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 12:20:39 2023

@author: lucina-may
"""
import numpy as np
from rdkit import Chem, DataStructs
from rdkit.DataManip import Metric
from rdkit.Chem import AllChem
from concerto_fp import get_biosynfoni

DEFAULT_BIOSYNFONI_VERSION = 'fps_full_3'

#circular fingerprint -------------------------------------------------
def morgan_getter(
        mol:Chem.Mol, 
        useChirality:bool = False, 
        radius:int = 2, 
        nBits:int = 2048
        ) -> DataStructs.ExplicitBitVect:
    """returns explicit bit vector fp"""
    
    fingerprint = AllChem.GetMorganFingerprintAsBitVect(
        mol,
        useChirality=useChirality,
        radius=radius,
        nBits=nBits)
    return fingerprint

#topology path fingerprint ---------------------------------------------
def rdk_fp_getter(
        mol:Chem.Mol, 
        nbits:str = 2048
        ) -> DataStructs.ExplicitBitVect:
    """returns explicit bit vector fp"""
    
    fingerprint = Chem.RDKFingerprint(mol, fpSize=nbits)
    return fingerprint

#substructure key fingerprints ------------------------------------------
def maccs_getter(mol:Chem.Mol) -> DataStructs.ExplicitBitVect:
    """returns explicit bit vector fp """
    
    fingerprint  = AllChem.GetMACCSKeysFingerprint(mol)     
    return fingerprint


def biosynfoni_getter(
        mol:Chem.Mol, 
        version:str = DEFAULT_BIOSYNFONI_VERSION):
    """returns counted fingerprint list"""
    counted_fingerprint = get_biosynfoni(mol, version=version, 
                                         return_matches=False)
    return counted_fingerprint


def binosynfoni_getter(
        mol:Chem.Mol, 
        version:str = DEFAULT_BIOSYNFONI_VERSION
        ) -> DataStructs.ExplicitBitVect:
    """returns explicit bit vector"""
    counted = get_biosynfoni(mol, version=version, return_matches=False)
    binary=[]
    for i in counted:
        if i > 0:
            binary.append(1)
        elif i == 0:
            binary.append(0)
    assert (len(binary)==len(counted)), "error in obtaining binosynfoni"
    
    binosynfoni = list_to_bitvect(binary)
    return binosynfoni

#============================= distance, similarity =========================

def tanimoto(
        expl_bitvectors:list[DataStructs.ExplicitBitVect]
        ) -> float:
    
    array = Metric.rdMetricMatrixCalc.GetTanimotoDistMat(\
        expl_bitvectors)
    return array.tolist()[0]


def euclidean(
        expl_bitvectors:list[DataStructs.ExplicitBitVect]
        ) -> float:
    
    array =  Metric.rdMetricMatrixCalc.GetEuclideanDistMat(\
        expl_bitvectors)
    return array.tolist()[0]


def fingerprint_similarity(fp1: np.ndarray, fp2: np.ndarray) -> float:
    return (
        np.logical_and(fp1, fp2).sum() / (
            float(np.logical_or(fp1, fp2).sum())
        )
    )


#============================= formatting ===================================

def list_to_bitvect(lst:list):
    """returns DataStructs.cDataStructs.ExplicitBitVect type"""
    binarystring = ''.join([str(x) for x in lst])
    return DataStructs.cDataStructs.CreateFromBitString(binarystring)


def array_to_bitvect(numpy_array):
    return list_to_bitvect(numpy_array.tolist())

