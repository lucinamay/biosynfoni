"""Checks if biosynfoni gives same results for molecules and chirality-removed mols"""
import os, sys
import subprocess
from sys import argv

from rdkit import Chem

# os.chdir(os.path.abspath(os.path.join(sys.path[0], os.pardir)))

from biosynfoni.inoutput import outfile_namer
from biosynfoni.rdkfnx import get_mol_supplier, sdf_writr
from biosynfoni.concerto_fp import loop_over_supplier

runwd = os.getcwd()
original_suppl_loc = argv[1]
original_suppl_dir = "/".join(original_suppl_loc.split("/")[:-1])
original_suppl_name = original_suppl_loc.split("/")[-1].split(".")[0]
# run biosynfoni on all molecules in sdf
original_suppl = get_mol_supplier(original_suppl_loc)
chir_biosynf = loop_over_supplier(original_suppl)

# run biosynfoni on all molecules in sdf after removing chirality
nonchir_mols = []
for mol in original_suppl:
    Chem.RemoveStereochemistry(mol)
    nonchir_mols.append(mol)

outname = outfile_namer(original_suppl_name, "nonchir")
outfilename = f"{outname}.sdf"
original_suppl_dir = "/".join(original_suppl_loc.split("/")[:-1])
os.chdir(original_suppl_dir)
sdf_writr(nonchir_mols, outfilename)

nonchir_path = os.path.abspath("{original_suppl_dir}/{outfilename}")
nonchir_suppl = get_mol_supplier(nonchir_path)
nonchir_biosynf = loop_over_supplier(nonchir_suppl)

os.chdir(runwd)
minlength = len(chir_biosynf)
if len(chir_biosynf) != len(nonchir_biosynf):
    minlength = min(len(chir_biosynf), len(nonchir_biosynf))
    with open(outfile_namer("log.txt"), "a") as f:
        f.write(
            "chir and nonchir biosynfoni have different lengths, taking {minlength} molecules"
        )


for i in range(minlength):
    if chir_biosynf[i] != nonchir_biosynf[i]:
        curr_smiles = Chem.MolToSmiles(nonchir_mols[i])
        with open(outfile_namer("chiraldifferences.txt"), "a") as g:
            g.write(f"{i}\t{chir_biosynf[i]}\t{nonchir_biosynf[i]}\t{curr_smiles}\n")

with open(outfile_namer("log.txt"), "a") as f:
    f.write("successfully finished checking molecules")
