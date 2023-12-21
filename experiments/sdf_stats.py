import argparse
import os, sys
from enum import Enum

from tqdm import tqdm
import numpy as np
from rdkit import Chem


def cli():
    parser = argparse.ArgumentParser()

    parser.add_argument("sdf", type=str, help="input sdf file")
    args = parser.parse_args()
    
    return args


def main():
    args = cli()

    sdf_path = os.abspath(args.sdf)
    folder = os.path.dirname(sdf_path)
    name = os.path.basename(sdf_path).split(".")[0]
    # out_path = os.path.join(folder, name + '.smi')

    iwd = os.getcwd()
    os.chdir(folder)

    suppl = Chem.SDMolSupplier(sdf_path)

    properties = {}

    nones = np.zeros(len(suppl))
    nones.fill(np.nan)

    properties['molecular_weight'] = np.copy(nones)
    properties['n_atoms'] = np.copy(nones)
    properties['carbons'] = np.copy(nones)
    properties['nitrogens'] = np.copy(nones)
    properties['oxygens'] = np.copy(nones)
    properties['halogens'] = np.copy(nones)
    properties['heteroatoms'] = np.copy(nones) 

    properties['double_bonds'] = np.copy(nones)
    properties['triple_bonds'] = np.copy(nones)
    properties['aromatic_bonds'] = np.copy(nones)
    properties['rings'] = np.copy(nones)
    properties['ring_atoms'] = np.copy(nones)
    properties['ring_bonds'] = np.copy(nones)
    properties['rotatable_bonds'] = np.copy(nones)
    
    for i, mol in tqdm(enumerate(suppl), total = len(suppl)):
        if mol is None:
            continue
        properties['molecular_weight'][i] = Chem.Descriptors.MolWt(mol)
        properties['n_atoms'][i] = mol.GetNumAtoms()
        properties['carbons'][i] = mol.GetNumAtoms(Chem.rdchem.Element.C)
        properties['nitrogens'][i] = mol.GetNumAtoms(Chem.rdchem.Element.N)
        properties['oxygens'][i] = mol.GetNumAtoms(Chem.rdchem.Element.O)
        properties['halogens'][i] = mol.GetNumAtoms(Chem.rdchem.Element.F) + mol.GetNumAtoms(Chem.rdchem.Element.Cl) + mol.GetNumAtoms(Chem.rdchem.Element.Br) + mol.GetNumAtoms(Chem.rdchem.Element.I)
        properties['heteroatoms'][i] = mol.GetNumHeavyAtoms() - mol.GetNumAtoms(Chem.rdchem.Element.C)

        properties['double_bonds'][i] = mol.GetNumBonds(Chem.rdchem.BondType.DOUBLE)
        properties['triple_bonds'][i] = mol.GetNumBonds(Chem.rdchem.BondType.TRIPLE)
        properties['aromatic_bonds'][i] = mol.GetNumBonds(Chem.rdchem.BondType.AROMATIC)
        properties['rings'][i] = Chem.rdMolDescriptors.CalcNumRings(mol)
        properties['ring_atoms'][i] = Chem.rdMolDescriptors.CalcNumAromaticRings(mol)
        properties['ring_bonds'][i] = Chem.rdMolDescriptors.CalcNumAromaticRings(mol)
        properties['rotatable_bonds'][i] = Chem.rdMolDescriptors.CalcNumRotatableBonds(mol)


if __name__ == "__main__":
    main()
