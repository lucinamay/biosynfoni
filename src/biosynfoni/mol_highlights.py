from rdkit import Chem

from biosynfoni.concerto_fp import Biosynfoni
from biosynfoni.moldrawing import _get_highlight_loc_and_col, pathway_colours, Palette
from biosynfoni.rdkfnx import BiosynfoniVersion


def mol_to_highlight_mapping(mol: Chem.Mol):
    """
    input: mol (Chem.Mol)

    output: tuple of [0] atom_indexes, [1] bond_indexes, [2] atom_colors, [3] bond_colors

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
    # get get version information
    version = BiosynfoniVersion("default")
    colors_per_sub = version.get_subs_colors()
    # get fingerprint's matched atoms
    fp = Biosynfoni(
        mol,
        substructure_set=version.substructures,
        intrasub_overlap=False,
        intersub_overlap=False,
    )
    matches_atomidx = fp.matches

    # get highlight locations and colors and unpack
    atoms_bonds, atomcols_bondcols = _get_highlight_loc_and_col(
        mol, matches_atomidx, colors_per_sub
    )
    atom_indexes, bond_indexes = atoms_bonds
    atom_colors, bond_colors = atomcols_bondcols
    return atom_indexes, bond_indexes, atom_colors, bond_colors
