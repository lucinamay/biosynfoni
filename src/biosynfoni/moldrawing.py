"""

||||||||||||  ㅇㅅㅇ  |||||||||||||||
____________________________________

title: figuremaking.py              ||
created: 2023-10-27                 ||
author: Lucina-May Nollen           ||
institute: WUR Bioinformatics       ||
____________________________________

||||||||||||  ()()()  |||||||||||||||

description: functions for figuremaking
"""

from enum import Enum
import typing as ty
import logging

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D, MolsToGridImage

from biosynfoni.subkeys import substructureSmarts


pathway_colours = {
    "shikimate": (0.65, 0.51, 0.71),
    "acetate": (1, 0.55, 0.38),
    "mevalonate": (0.61, 0.76, 0.73),
    "methylerythritol_phosphate": (0.61, 0.76, 0.73),
    "sugar": (1, 0.77, 0.81),
    "amino_acid": (1, 0.92, 0.63),
}


class Palette(Enum):
    """
    Palette of colors for drawing molecules as RGB.
    """

    Red = (230, 25, 75)
    Blue = (0, 130, 200)
    Green = (60, 180, 75)
    Maroon = (128, 0, 0)
    Brown = (170, 110, 40)
    Olive = (128, 128, 0)
    Teal = (0, 128, 128)
    Navy = (0, 0, 128)
    Orange = (245, 130, 48)
    Yellow = (255, 225, 25)
    Lime = (210, 245, 60)
    Cyan = (70, 240, 240)
    Purple = (145, 30, 180)
    Magenta = (240, 50, 230)
    Pink = (255, 190, 212)
    Apricot = (255, 215, 180)
    Beige = (255, 250, 200)
    Mint = (170, 255, 195)
    Lavender = (220, 190, 255)

    def as_hex(self, alpha: ty.Optional[float] = None) -> str:
        """
        Return the hex code of the color with the given alpha value.

        Args:
            alpha (float): Alpha value of the color, default is None.

        Returns:
            str: Hex code of the color.
        """
        r, g, b = self.value
        hex_base = "#{:02x}{:02x}{:02x}".fomat(r, g, b)
        if alpha is not None:
            if alpha < 0:
                alpha = 0.0
            elif alpha > 1:
                alpha = 1.0
            return "{}{:02x}".format(hex_base, int(alpha * 255))
        else:
            return hex_base

    def normalize(
        self, minv: float = 0.0, maxv: float = 255.0
    ) -> ty.Tuple[float, float, float]:
        """
        Normalize the color values to the range [0, 1].

        Args:
            minv (float, optional): Minimum value of the range. Defaults to 0.0.
            maxv (float, optional): Maximum value of the range. Defaults to 255.0.

        Returns:
            ty.Tuple[float, float, float]: Normalized color values.
        """
        r, g, b = self.value
        return (
            round((r - minv) / (maxv - minv), 3),
            round((g - minv) / (maxv - minv), 3),
            round((b - minv) / (maxv - minv), 3),
        )


# =============


def _get_highlight_loc_and_col(
    mol: Chem.Mol,
    subs_matches_for_highlighting: list = [],  # list of lists of lists of atom indices
    subs_colors: ty.Optional[ty.List[ty.Tuple[float, float, float]]] = None,
) -> tuple[tuple[list, dict]]:
    """gets highlighting information"""
    palette = [c.normalize() for c in Palette]
    atoms_to_highlight, bonds_to_highlight = [], []
    atom_highlight_colors, bond_highlight_colors = {}, {}

    if subs_colors is None:
        subs_colors = [
            palette[i % len(palette)] for i in range(len(subs_matches_for_highlighting))
        ]
    else:
        assert len(subs_colors) == len(
            subs_matches_for_highlighting
        ), "Number of colors must match number of substructures."

    for sub_index, sub_matches in enumerate(subs_matches_for_highlighting):
        atom_indices = []
        color = subs_colors[sub_index]
        for match in sub_matches:
            # individual substructure-match indices
            match_indices = []

            if len(match) == 2:
                logging.debug(match)

            for atom_index in match:
                atom_indices.append(atom_index)
                atoms_to_highlight.append(atom_index)
                atom_highlight_colors[atom_index] = color
                # get only the indices for current match
                match_indices.append(atom_index)

            for bond in mol.GetBonds():
                bond_index = bond.GetIdx()
                start_atom = bond.GetBeginAtom().GetIdx()
                end_atom = bond.GetEndAtom().GetIdx()

                # if start_atom in atom_indices and end_atom in atom_indices:
                if start_atom in match_indices and end_atom in match_indices:
                    bonds_to_highlight.append(bond_index)
                    bond_highlight_colors[bond_index] = color
                    logging.debug(match, bonds_to_highlight)
    locations = (atoms_to_highlight, bonds_to_highlight)
    colors = (atom_highlight_colors, bond_highlight_colors)
    return locations, colors


# ================================================


def draw(
    mol: Chem.Mol,
    window_size: ty.Tuple[int, int] = (800, 800),
    background_color: ty.Optional[str] = None,
    highlight_atoms_bonds_mappings=None,
) -> str:
    """
    Draw a molecule with its substructures highlighted.

    Args:
        mol (Chem.Mol): Molecule.
        subs (ty.List[tuple[tuple[int]]]], optional): List of substructure indices (from GetStructureMatches) to highlight. Defaults to [].
        window_size (ty.Tuple[int, int], optional): Window size. Defaults to (800, 800).
        background_color (ty.Optional[str], optional): Background color. Defaults to None.

    Returns:
        str: SVG string.
    """

    drawing = rdMolDraw2D.MolDraw2DSVG(*window_size)

    # unpack
    (
        atoms_to_highlight,
        bonds_to_highlight,
        atom_highlight_colors,
        bond_highlight_colors,
    ) = highlight_atoms_bonds_mappings

    options = drawing.drawOptions()
    if background_color is not None:
        options.setBackgroundColour(background_color)
    options.useBWAtomPalette()
    drawing.DrawMolecule(
        mol,
        highlightAtoms=atoms_to_highlight,
        highlightBonds=bonds_to_highlight,
        highlightAtomColors=atom_highlight_colors,
        highlightBondColors=bond_highlight_colors,
    )

    drawing.FinishDrawing()
    svg_str = drawing.GetDrawingText().replace("svg:", "")

    return svg_str


def drawfp(
    subs_set: list[Chem.Mol],
    subs_colors: list,
    subs_labels: list[str],
    window_size: ty.Tuple[int, int] = (1000, 1000),
) -> str:
    mols, indexes = [], []
    subs_atoms_to_highlight, subs_bonds_to_highlight = [], []
    subs_atom_highlight_colors, subs_bond_highlight_colors = [], []
    for sub_index, sub_mol in enumerate(subs_set):
        if sub_mol:
            emulate_match_for_highlighting = [[[]] for each_sub in subs_set]
            emulate_match_for_highlighting[sub_index] = [
                [atom.GetIdx() for atom in sub_mol.GetAtoms()]
            ]
            logging.debug(emulate_match_for_highlighting)
            loc, col = _get_highlight_loc_and_col(
                sub_mol, emulate_match_for_highlighting, subs_colors
            )
            atoms_to_highlight, bonds_to_highlight = loc
            atom_highlight_colors, bond_highlight_colors = col
            subs_atoms_to_highlight.append(atoms_to_highlight)
            subs_bonds_to_highlight.append(bonds_to_highlight)
            subs_atom_highlight_colors.append(atom_highlight_colors)
            subs_bond_highlight_colors.append(bond_highlight_colors)

            mols.append(sub_mol)
            indexes.append(sub_index)

    successful_subs = [subs_labels[i] for i in indexes]
    names = [sub_label.replace("_", " ") for sub_label in successful_subs]
    logging.debug(len(mols), len(names), len(subs_atoms_to_highlight))
    logging.debug(subs_atoms_to_highlight)
    subs_atoms_to_highlight = [[], [], [], [], []]
    grid_image = MolsToGridImage(
        mols,
        molsPerRow=6,
        subImgSize=window_size,
        legends=names,
        highlightAtomLists=[[i for i in range(len(mol.GetAtoms()))] for mol in mols],
        highlightAtomColors=subs_atom_highlight_colors,
        highlightBondColors=subs_bond_highlight_colors,
        useSVG=True,
    )
    svg_str = grid_image.replace("svg:", "")
    logging.debug(grid_image[:10], svg_str, sep="\n")
    return svg_str
