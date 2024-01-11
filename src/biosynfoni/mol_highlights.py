from rdkit import Chem

from biosynfoni.concerto_fp import Biosynfoni
from biosynfoni.moldrawer import _get_highlight_loc_and_col
from biosynfoni.rdkfnx import BiosynfoniVersion


def mol_to_highlight_mapping(mol: Chem.Mol):
    version = BiosynfoniVersion("default")
    bsf = Biosynfoni(mol, substructure_set=version.substructures)
    matches = bsf.matches
