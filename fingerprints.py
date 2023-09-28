#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 12:20:39 2023

@author: lucina-may
"""
import numpy as np
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
import def_biosynfoni
from concerto_fp import get_biosynfoni

DEFAULT_BIOSYNFONI_VERSION = 'fps3_full'

def morgan_getter(mol, useChirality=False, radius=2, nBits=2048):
    fingerprint = AllChem.GetMorganFingerprintAsBitVect(mol,
                                                   useChirality=useChirality,
                                                   radius=radius,
                                                   nBits=nBits)
    morgan = np.array(fingerprint)
    return morgan

def maccs_getter(mol):
    fingerprint  = AllChem.GetMACCSKeysFingerprint(mol)
    maccs = np.array(fingerprint)
    return maccs

def rdk_fp_getter(mol, nbits=2048):
    fingerprint = Chem.RDKFingerprint(mol, fpSize=nbits)
    daylight_like = np.array(fingerprint)
    return daylight_like


def biosynfoni_getter(mol, version=DEFAULT_BIOSYNFONI_VERSION):
    fingerprint = get_biosynfoni(mol, version=version,matches=False)
    biosynfoni = np.array(fingerprint)
    return biosynfoni

def binosynfoni_getter(mol, version=DEFAULT_BIOSYNFONI_VERSION):
    counted = get_biosynfoni(mol, version=version, matches=False)
    binary=[]
    for i in counted:
        if i > 0:
            binary.append(1)
        elif i == 0:
            binary.append(0)
    assert (len(binary)==len(counted)), "error in obtaining binosynfoni"
    binosynfoni = np.array(binary)
    return binosynfoni
            