# -*- coding: utf-8 -*-
"""
|||||||||||||  ㅇㅅㅇ  ||||||||||||||||
____________________________________

title: biosyn distance comparison   ||
creaetd: 2023.09.28 09:49           ||
author: lucina-may nollen           ||
institute: WUR Bioinformatics       ||
____________________________________

||||||||||||  ()()()  |||||||||||||||

description: uses molecule pairs for comparison of distance estimation
between various fingerprints,
supports:
    - morgan
    - rdkit
    - maccs
    - biosynfoni
    - binosynfoni

also supports various distance metrics:
    - countanimoto
    - tanimoto
    - euclidean


# ideally, do not have any metacyc-related things in here, but in metacyc_extract.py
"""
# standard packages
import sys
import typing as ty
from enum import Enum, auto
from sys import argv

# packages requiring installing
import pandas as pd
import numpy as np
import plotly.express as px
from rdkit import Chem

# from plotly.offline import download_plotlyjs, init_notebook_mode, iplot
# own imports #fix to work with different folders

sys.path.append("../src/")

from biosynfoni import fingerprints as fp
from biosynfoni import picklr, jaropener, outfile_namer, save_version
from biosynfoni import entryfile_dictify as ann
from metacyc_better_taxonomy import BETTER_TAX
from utils import figuremaking as fm

# =========================== GLOBALS =========================
FP_FUNCTIONS = {
    "biosynfoni": fp.biosynfoni_getter,
    "binosynfoni": fp.binosynfoni_getter,
    "maccs": fp.maccs_getter,
    "morgan": fp.morgan_getter,
    "rdkit": fp.rdk_fp_getter,
}

SIM_FUNCTIONS = {
    "c_tanimoto": fp.countanimoto,
    "manhattan": fp.manhattan,
    "tanimoto_dist": fp.bitvect_to_tanimoto,
    "euclidean_dist": fp.bitvect_to_euclidean,
}

# for later:


"""if not isinstance(fingerprint, FingerprintType):
    raise TypeError(
        f"wrong signature for fingerprint: {type(fingerprint)} != FingerprintType" )
"""
# fingerprint.values.lower()


# =============================================================================

# =============================== distance fnx ================================


def _fps_to_distance(fp1, fp2, metric="c_tanimoto"):
    """for given molecule pair fp1 and fp2, gives distance"""
    # input handling & checking
    assert isinstance(fp1, np.ndarray) and isinstance(
        fp2, np.ndarray
    ), "please provide two fingerprints"
    distance = None  # init
    distance = SIM_FUNCTIONS[metric]([fp1, fp2])
    return distance


def similarity(
    mol1: Chem.Mol,
    mol2: Chem.Mol,
    fingerprint: str,
    metric: str = "c_tanimoto",
    check_moliness: bool = True,
) -> float:
    fingerprint = fingerprint.lower()
    metric = metric.lower()

    assert fingerprint in FP_FUNCTIONS.keys(), "unsupported fingerprint: {}".format(
        fingerprint
    )
    assert metric in SIM_FUNCTIONS.keys(), "unsupported metric: {}".format(metric)
    if check_moliness:
        assert isinstance(mol1, Chem.Mol) and isinstance(
            mol2, Chem.Mol
        ), "please provide two Chem.Mol-type molecules"

    if not (isinstance(mol1, Chem.Mol) and isinstance(mol2, Chem.Mol)):
        return np.nan

    fp1, fp2 = np.array([]), np.array([])  # init.
    fp1 = FP_FUNCTIONS[fingerprint](mol1)
    fp2 = FP_FUNCTIONS[fingerprint](mol2)

    return _fps_to_distance(fp1, fp2, metric=metric)


# =============================================================================

# =============================== getting pairs ===============================


def get_step_pairs(chain_mols: list, degree_of_sep: int = 1) -> list:
    """for given list of molecules originating from one reaction,
    gets the available(non-null) molecule pairs separated by
    degree of separation"""

    pairs = []
    for i in range(len(chain_mols) - degree_of_sep):
        if chain_mols[i] and chain_mols[i + degree_of_sep]:
            pairs.append([chain_mols[i], chain_mols[i + degree_of_sep]])
    return pairs


def yield_row_from_file(struct_loc: str):
    """yields rows from csv, cleaning up the rows in the process"""
    struct = pd.read_csv(struct_loc, sep="\t", header=None)
    # get headers
    nums = ["{}".format(i - 1) for i in range(len(struct.columns))]
    nums.pop(0)
    headers = ["pathway"] + nums
    struct.columns = headers
    # clean up
    struct.replace("", np.nan, inplace=True)  # just in case one was missed
    filtered = struct.dropna(thresh=2)  # drop completely empty ones
    final = filtered.replace(np.nan, "")  # replace back for later fnxs
    # return pathways
    for i in range(len(final.index)):
        yield final.iloc[i].tolist()


def dictify_pw_pairs(struct_loc: str, degree_of_sep: int = 1) -> dict[list[list[str]]]:
    """returns {'pathway-id':[[product,precursor],[precursor,preprecursor]]}"""
    dictionary = {}
    for row in yield_row_from_file(struct_loc):
        pathway_id = row[0]
        pairs = get_step_pairs(row[1:], degree_of_sep=degree_of_sep)
        dictionary[pathway_id] = pairs
    return dictionary


def _label_pairs(pairs_per_pathway: dict[list[list[str]]]) -> tuple[str, list[str]]:
    """returns list of (pathway,distances)"""
    labeled_distances = []
    for pw_id in pairs_per_pathway.keys():
        for distance in pairs_per_pathway[pw_id]:
            labeled_distances.append((pw_id, distance))
    return labeled_distances


def str_to_mol(repres: str, clean: bool = True) -> Chem.Mol:
    mol = ""
    if repres.startswith("InChI="):
        mol = Chem.MolFromInchi(repres)
    elif clean:
        mol = Chem.MolFromSmiles(
            repres.split(
                "[a ",
            )[0]
        )
    else:
        mol = Chem.MolFromSmiles(repres)
    if mol:
        return mol
    else:
        return ""


def get_pairs_dist_df(
    labeled_pairs: list[tuple[str, list[str]]], metric="c_tanimoto"
) -> pd.DataFrame:
    """only supports inchi_pairs at the moment"""
    prec_strings = [x[1][0] for x in labeled_pairs]
    prod_strings = [x[1][1] for x in labeled_pairs]
    pathways = [x[0] for x in labeled_pairs]

    df = pd.DataFrame(
        {"pathway": pathways, "precursor": prec_strings, "product": prod_strings}
    )

    prec_mols = [str_to_mol(x) for x in prec_strings]
    prod_mols = [str_to_mol(x) for x in prod_strings]
    molpairs = [[prec_mols[i], prod_mols[i]] for i in range(len(prec_mols))]

    # find way to convert to mol only once instead of continuously
    for fp_type in FP_FUNCTIONS.keys():
        df[fp_type] = [
            similarity(
                x[0], x[1], fingerprint=fp_type, metric=metric, check_moliness=False
            )
            for x in molpairs
        ]
    return df


def _get_comparison_df(struct_loc, max_steps=4, metric="c_tanimoto"):
    current_df = pd.DataFrame()
    for i in range(1, (max_steps + 1), 1):
        previous_df = current_df.copy()
        degree_of_sep = i
        dictionary = dictify_pw_pairs(struct_loc, degree_of_sep=degree_of_sep)
        labeled_pairs = _label_pairs(dictionary)
        # labeled_pairs = _label_pairs(dictionary)[:100]
        current_df = get_pairs_dist_df(labeled_pairs, metric=metric)
        current_df["reaction_steps"] = i
        current_df = pd.concat([previous_df, current_df], axis=0)
    return current_df


def annotate_pathways(df, pw_tax_file: str, tax_text_file: str) -> pd.DataFrame:
    ndf = df.copy()
    start_val_sep = " - "
    entry_sep = "//"
    encoding = "iso8859-1"
    pw_tax = ann(
        pw_tax_file,
        keyvals=("UNIQUE-ID", "TAXONOMIC-RANGE"),
        start_val_sep=start_val_sep,
        encoding=encoding,
        entry_sep=entry_sep,
    )

    class_text = ann(
        tax_text_file,
        keyvals=("UNIQUE-ID", "COMMON-NAME"),
        start_val_sep=start_val_sep,
        encoding=encoding,
        entry_sep=entry_sep,
    )
    print(len(pw_tax), len(class_text))
    ndf["tax_id"] = ndf["pathway"].replace(to_replace=pw_tax)
    ndf["taxonomic_range"] = ndf["tax_id"].replace(to_replace=class_text)
    ndf["taxonomy"] = ndf["taxonomic_range"].replace(to_replace=BETTER_TAX)
    ndf["pathway_name"] = ndf["pathway"].replace(to_replace=class_text)
    # ndf['pathway_name']=ndf['pathway_name'].fillna(ndf['pathway'])
    return ndf


def get_fp_combinations() -> list[tuple[str]]:
    combinations = []
    for firstkey in FP_FUNCTIONS.keys():
        for secondkey in FP_FUNCTIONS.keys():
            if firstkey != secondkey:
                if (firstkey, secondkey) not in combinations:
                    if (secondkey, firstkey) not in combinations:
                        combinations.append((firstkey, secondkey))
    return combinations


def get_square(
    df: pd.DataFrame, col1: str, col2: str, range1: tuple[float], range2: tuple[float]
) -> pd.DataFrame:
    """returns the molecule pair's inchis for 'dots' in a given square of the scatter plot"""
    square = df[
        (df[col1] > range1[0])
        & (df[col1] < range1[1])
        & (df[col2] > range2[0])
        & (df[col2] < range2[1])
    ]
    return square


def draw_molpair(
    pair: list[Chem.Mol], outfilename: str, highlighting: bool = True
) -> None:
    for i in range(len(pair)):
        svg_text = fm.draw(pair[i], get_match_highlighting=highlighting)
        with open(f"{outfilename}_{i}.svg", "w") as f:
            f.write(svg_text)
    return None


def draw_squares(
    square_df: pd.DataFrame,
    pair_columns: str = ("precursor", "product"),
    squarename: str = "origin",
    highlighting: bool = True,
) -> None:
    """draws the molecules in the squares
    input: (pd.DataFrame) square_df -- the dataframe containing the squares
    (str) pair_columns -- the name of the column containing the molecule pairs
    (str) squarename -- the name of the square
    """
    for i in square_df.index:
        pair = [
            str_to_mol(square_df[pair_columns[0]][i]),
            str_to_mol(square_df[pair_columns[1]][i]),
        ]
        pathway = square_df["pathway"][i]
        outfilename = outfile_namer(f"{squarename}_{pathway}")
        draw_molpair(pair, outfilename, highlighting=highlighting)
    return None


def main():
    # struct_loc = argv[1]
    struct_loc = "../../../scripts/0927_metacyc_reactions.tsv"
    degree_of_sep = 1
    pw_tax_file = "../../../metacyc/pathways_taxid.txt"
    tax_text_file = "../../../metacyc/cleaner_classes.dat"
    biosynfoni_version = fp.DEFAULT_BIOSYNFONI_VERSION
    metric = "c_tanimoto"
    # get all the reaction pairs
    comparison_df = _get_comparison_df(struct_loc, max_steps=4)
    annotated_df = annotate_pathways(comparison_df, pw_tax_file, tax_text_file)

    annotated_df["stepnum"] = annotated_df["reaction_steps"].apply(str)
    fp_combs = get_fp_combinations()
    for com in fp_combs:
        fm.df_scatterplot(
            annotated_df,
            com[0],
            com[1],
            figtitle=f"{metric} for different reaction step numbers",
            filename=outfile_namer(f"{com[0]}_{com[1]}_{metric}"),
            color="stepnum",
            color_discrete_map=fm.COLOUR_DICT["stepnum"],
            hover_data=["pathway_name", "taxonomic_range", "taxonomy"],
            marginal_x="box",
            marginal_y="box",
        )
    for fp_type in ["rdkit", "maccs", "morgan"]:
        one_step = annotated_df[annotated_df["stepnum"] == "1"]
        left_bottom = get_square(one_step, "biosynfoni", fp_type, (0, 0.2), (0, 0.2))
        left_top = get_square(one_step, "biosynfoni", fp_type, (0, 0.2), (0.8, 1))
        right_bottom = get_square(one_step, "biosynfoni", fp_type, (0.8, 1), (0, 0.2))
        # right_top = get_square(one_step, "biosynfoni", fp_type, (0.8, 1), (0.8, 1))
        draw_squares(left_bottom, squarename=f"{fp_type}_origin")
        draw_squares(left_top, squarename=f"{fp_type}_left_top")

    # draw_squares(right_bottom, squarename='right_bottom')
    # draw_squares(right_top, squarename='right_top')
    # plot_two_cols(annotated_df,'biosynfoni', 'rdkit',symbol=None)

    # save the biosynfoni version for reference

    save_version(biosynfoni_version)

    #
    # annotate_pathways(df, pw_tax_file, tax_text_file)
    # plot_two_cols(df,col1,col2,colour_label='taxonomic_range',shape_label='',\
    #              hover_data=['pathway_name','taxonomic_range'])


if __name__ == "__main__":
    main()
