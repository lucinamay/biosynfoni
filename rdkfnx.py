#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 14:59:28 2023

@author: lucina-may
"""
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit import RDLogger  #for muting warnings

def sdf_writr(mols:list, outfile:str) -> None:
    """writes sdf of mols, 
    * needs RDKit Chem
    """
    writer = Chem.SDWriter(outfile)
    for mol in mols:
        writer.write(mol)
    return None