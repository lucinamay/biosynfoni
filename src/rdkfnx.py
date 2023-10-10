#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 14:59:28 2023
||||||||||||  ㅇㅅㅇ  |||||||||||||||
____________________________________

title: RDKit functions              ||
created: 2023-09-13                 ||
author: Lucina-May Nollen           || 
institute: WUR Bioinformatics       ||
____________________________________
 
||||||||||||  ()()()  |||||||||||||||

contains general reusable functionalities depending on RDKit packages,
mainly supplier-handling
"""
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit import RDLogger  # for muting warnings


def sdf_writr(mols: list, outfile: str) -> None:
    """writes sdf of mols,
    * needs RDKit Chem
    """
    writer = Chem.SDWriter(outfile)
    for mol in mols:
        writer.write(mol)
    return None


def get_supplier(sdf_file: str, supplier_only: bool = True) -> list:
    suppl = Chem.SDMolSupplier(sdf_file)

    print("reading sdf, number of entries:", len(suppl))

    if supplier_only:
        return suppl
    else:
        mols = [x for x in suppl]
        return mols
