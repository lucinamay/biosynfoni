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

import plotly.express as px
import matplotlib.pyplot as plt
import pandas as pd
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

from concerto_fp import get_biosynfoni
from def_biosynfoni import DEFAULT_BIOSYNFONI_VERSION
from def_biosynfoni import get_subsset as gss

# Based off of David's code ====================================================


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


def draw(
    mol: Chem.Mol,
    subs: ty.List[ty.List[tuple[tuple[int]]]] = [],
    window_size: ty.Tuple[int, int] = (800, 800),
    background_color: ty.Optional[str] = None,
    get_match_highlighting: bool = False,
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
    if get_match_highlighting:
        _, subs = get_biosynfoni(
            mol, version=DEFAULT_BIOSYNFONI_VERSION, return_matches=True
        )

    drawing = rdMolDraw2D.MolDraw2DSVG(*window_size)
    palette = [c.normalize() for c in Palette]

    atoms_to_highlight = []
    bonds_to_highlight = []
    atom_highlight_colors = {}
    bond_highlight_colors = {}

    for i, sub in enumerate(subs):
        atom_indices = []
        color = palette[i % len(palette)]

        for match in sub:
            # individual substructure-match indices
            match_indices = []
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
                    # print(match, bonds_to_highlight)

    options = drawing.drawOptions()
    if background_color is not None:
        options.setBackgroundColour(background_color)
    options.useBWAtomPalette()
    # options.addAtomIndices = True
    # options.addBondIndices = True
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
    version: str = DEFAULT_BIOSYNFONI_VERSION,
    window_size: ty.Tuple[int, int] = (1000, 1000),
) -> str:
    drawing = rdMolDraw2D.MolDraw2DSVG(*window_size)
    mols, indexes = [], []
    for i in range(len(gss(version))):
        if gss(version)[i]:
            mols.append(gss(version)[i])
            indexes.append(i)
    drawing.DrawMolecules(mols)
    drawing.FinishDrawing()
    svg_str = drawing.GetDrawingText().replace("svg:", "")
    return svg_str


