from rdkit import Chem

from biosynfoni.concerto_fp import Biosynfoni
from biosynfoni.moldrawing import _get_highlight_loc_and_col, pathway_colours, Palette
from biosynfoni.rdkfnx import BiosynfoniVersion


def mol_to_highlight_mapping(mol: Chem.Mol):
    """
    input: mol (Chem.Mol)

    output to be used like:
    drawing = rdMolDraw2D.MolDraw2DSVG()
    drawing.DrawMolecule(
        mol,
        highlightAtoms=atom_indexes,
        highlightBonds=bond_indexes,
        highlightAtomColors=atom_colors,
        highlightBondColors=bond_colors,
    )
    """
    version = BiosynfoniVersion("default")
    print(type(version))
    fp = Biosynfoni(mol, substructure_set=version.substructures)
    matches = fp.matches
    colors = version.get_subs_colors()
    atoms_bonds, atomcols_bondcols = _get_highlight_loc_and_col(mol, matches, colors)
    atom_indexes, bond_indexes = atoms_bonds
    atom_colors, bond_colors = atomcols_bondcols
    return atom_indexes, bond_indexes, atom_colors, bond_colors
