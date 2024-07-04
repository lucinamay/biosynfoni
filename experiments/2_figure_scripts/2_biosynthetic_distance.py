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
    - overlap_binosynfoni  

also supports various distance metrics:
    - countanimoto
    - tanimoto
    - euclidean


# ideally, do not have any metacyc-related things in here, but in metacyc_extract.py
"""
# standard packages
import os, argparse, logging
from enum import Enum
from functools import partial

# packages requiring installing
import numpy as np
import pandas as pd
from rdkit import Chem
from tqdm import tqdm

# from plotly.offline import download_plotlyjs, init_notebook_mode, iplot
# own imports #fix to work with different folders

# from metacyc_better_taxonomy import BETTER_TAX
from experiments.figure_producing.utils import figures as fm
from utils import set_style

from biosynfoni.subkeys import defaultVersion
from biosynfoni import fingerprints as fp
from biosynfoni import moldrawing
from biosynfoni.rdkfnx import save_version
from biosynfoni.inoutput import outfile_namer, output_direr
from biosynfoni.inoutput import entryfile_dictify as ann
from biosynfoni import get_highlight_mapping

# =========================== GLOBALS =========================
FP_FUNCTIONS = {
    "biosynfoni": fp.biosynfoni_getter,
    "interoverlap_biosynfoni": partial(fp.biosynfoni_getter, intersub_overlap=True),
    "overlap_biosynfoni": partial(
        fp.biosynfoni_getter, intersub_overlap=True, intrasub_overlap=True
    ),
    "overlap_binosynfoni": partial(
        fp.binosynfoni_getter, intersub_overlap=True, intrasub_overlap=True
    ),
    "maccs": fp.maccs_getter,
    "morgan": fp.morgan_getter,
    "rdkit": fp.rdk_fp_getter,
    "maccsynfoni": fp.maccsynfoni_getter,
    "overlap": partial(
        fp.bino_maccs_getter, intersub_overlap=True, intrasub_overlap=True
    ),
}

SIM_FUNCTIONS = {
    "c_tanimoto": fp.countanimoto,
    "cosine_sim": fp.cosine_sim,
    "manhattan_sim": fp.manhattan_sim,
    "tanimoto_sim": fp.tanimoto_sim,
    "euclidean_sim": fp.euclidean_sim,
}

# for later:


"""if not isinstance(fingerprint, FingerprintType):
    raise TypeError(
        f"wrong signature for fingerprint: {type(fingerprint)} != FingerprintType" )
"""
# fingerprint.values.lower()


BETTER_TAX = {
    "Viridiplantae": "Viridiplantae",
    "Bacteria <bacteria>": "Bacteria",
    "Pseudomonadota": "Bacteria",
    "Eukaryota": "Eukaryota",
    "Magnoliopsida": "Viridiplantae",
    "cellular organisms": "Cellular organisms",
    "Bryophyta": "Viridiplantae",
    "Pimpinella": "Viridiplantae",
    "Brassicaceae": "Viridiplantae",
    "Equisetum": "Viridiplantae",
    "Spermatophyta": "Viridiplantae",
    "Fabaceae": "Viridiplantae",
    "Tracheophyta": "Viridiplantae",
    "Fungi": "Fungi",
    "Mammalia": "Metazoa",
    "Gunneridae": "Viridiplantae",
    "Capsicum": "Viridiplantae",
    "asterids": "Viridiplantae",
    "Solanaceae": "Viridiplantae",
    "Rutaceae": "Viridiplantae",
    "Amaryllidaceae": "Viridiplantae",
    "Ephedra": "Viridiplantae",
    "Ranunculales": "Viridiplantae",
    "Papaveraceae": "Viridiplantae",
    "Berberidaceae": "Viridiplantae",
    "Caryophyllaceae": "Viridiplantae",
    "Homo": "Metazoa",
    "Archaea": "Archaea",
    "Methanococci": "Archaea",
    "Penicillium": "Fungi",
    "Lathyrus": "Viridiplantae",
    "Embryophyta": "Viridiplantae",
    "Clavicipitaceae": "Fungi",
    "Streptomyces": "Bacteria",
    "Solanum": "Viridiplantae",
    "Cyanobacteriota": "Bacteria",
    "Bacilli": "Bacteria",
    "Amygdaloideae": "Viridiplantae",
    "Poaceae": "Viridiplantae",
    "Lamiaceae": "Viridiplantae",
    "Vertebrata <vertebrates>": "Metazoa",
    "Allioideae": "Viridiplantae",
    "Caryophyllales": "Viridiplantae",
    "Drosophila <flies,genus>": "Metazoa",
    "Streptomycetales": "Bacteria",
    "Chlorophyta": "Viridiplantae",
    "Euphyllophyta": "Viridiplantae",
    "Insecta": "Metazoa",
    "Lepidoptera": "Metazoa",
    "Pinidae": "Viridiplantae",
    "Gossypium": "Viridiplantae",
    "Santalum": "Viridiplantae",
    "Shewanellaceae": "Bacteria",
    "Thermoprotei": "Archaea",
    "Chromatiaceae": "Bacteria",
    "Iridaceae": "Viridiplantae",
    "eudicotyledons": "Viridiplantae",
    "lamiids": "Viridiplantae",
    "Myxococcales": "Bacteria",
    "Actinomycetota": "Bacteria",
    "Shewanella": "Bacteria",
    "Cyperus": "Viridiplantae",
    "Cannabaceae": "Viridiplantae",
    "Cannabis": "Viridiplantae",
    "Hypericum": "Viridiplantae",
    "Metazoa": "Metazoa",
    "Sphingomonadales": "Bacteria",
    "Haptophyta": "Viridiplantae",
    "Actinomycetes": "Bacteria",
    "Apocynaceae": "Viridiplantae",
    "Protostomia": "Metazoa",
    "Enterobacterales": "Bacteria",
    "Erythroxylaceae": "Viridiplantae",
    "Halobacteria": "Archaea",
    "Escherichia coli": "Bacteria",
    "Streptomycetaceae": "Bacteria",
    "Bacillota": "Bacteria",
    "Coffea": "Viridiplantae",
    "Opisthokonta": "Opisthokonta",  # animal/fung
    "Boraginaceae": "Viridiplantae",
    "Sorghum": "Viridiplantae",
    "Mesangiospermae": "Viridiplantae",
    "Rosaceae": "Viridiplantae",
    "Marchantiophyta": "Viridiplantae",
    "Rhodophyta": "Viridiplantae",
    "Brassicales": "Viridiplantae",
    "Bacteroidales": "Bacteria",
    "Helicobacter pylori": "Bacteria",
    "Corynebacteriales": "Bacteria",
    "Ascomycota": "Fungi",
    "Haemodoraceae": "Viridiplantae",
    "Lecanorineae": "Fungi",
    "Lactobacillales": "Bacteria",
    "Theaceae": "Viridiplantae",
    "Rubiaceae": "Viridiplantae",
    "Rhizobiaceae": "Bacteria",
    "Methanobacteria": "Archaea",
    "Populus": "Viridiplantae",
    "Porphyromonas gingivalis": "Bacteria",
    "Lampyridae": "Metazoa",
    "Gammaproteobacteria": "Bacteria",
    "Mycobacteriaceae": "Bacteria",
    "Enterobacteriaceae": "Bacteria",
    "Arthropoda": "Metazoa",
    "Thermococci": "Archaea",
    "Chlamydia": "Bacteria",
    "Mimoseae": "Viridiplantae",
    "Bacillariophyta": "Viridiplantae",
    "Apiales": "Viridiplantae",
    "Avena": "Viridiplantae",
}


# =============================================================================


def cli():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "structures_path",
        type=str,
        help="path to the metacyc reaction chains with structure strings",
    )
    parser.add_argument(
        "-m",
        "--metric",
        choices=SIM_FUNCTIONS.keys(),
        default="c_tanimoto",
        help="distance metric to use",
    )
    parser.add_argument(
        "-l",
        "--lessblock",
        action="store_true",
        default=False,
        help="if passed, removes only intra-substructure overlapping matches (for biosynfoni)",
    )
    parser.add_argument(
        "-n",
        "--noblock",
        default=False,
        action="store_true",
        help="if passed, allows all overlapping matches (for biosynfoni)",
    )
    parser.add_argument(
        "-b",
        "--binary",
        default=False,
        action="store_true",
        help="if passed, uses binary fingerprints (a.k.a. overlap_binosynfoni)",
    )
    parser.add_argument(
        "-a",
        "--annotate",
        type=str,
        nargs=2,
        default=None,
        help="paths to the pathawy_to_taxonomy and taxonomy_to_humanlang files. if passed, annotates the pathways with taxonomic information",
    )
    args = parser.parse_args()

    # args.add_text = ""
    # args.intersub_overlap, args.intrasub_overlap = False, False
    # if args.lessblock:
    #     args.intersub_overlap = True
    #     args.intrasub_overlap = False
    #     args.add_text = "lessblock"
    # if args.noblock:
    #     args.intersub_overlap, args.intrasub_overlap = True, True
    #     args.add_text = "noblock"

    return args


# =============================== distance fnx ================================


def _fps_to_distance(fp1, fp2, metric="c_tanimoto") -> float:
    """for given molecule pair fp1 and fp2, gives distance"""
    # input handling & checking
    assert isinstance(fp1, np.ndarray) and isinstance(
        fp2, np.ndarray
    ), "please provide two fingerprints"
    distance = None  # init
    distance = SIM_FUNCTIONS[metric]([fp1, fp2])
    return distance


def _get_fp(
    mol1: Chem.Mol,
    mol2: Chem.Mol,
    fingerprint: str,
    *args,
    **kwargs,
) -> float:
    """args and kwargs are passed to get_biosynfoni"""
    fingerprint = fingerprint.lower()

    assert fingerprint in FP_FUNCTIONS.keys(), f"unsupported fingerprint: {fingerprint}"

    if not (isinstance(mol1, Chem.Mol) and isinstance(mol2, Chem.Mol)):
        return np.nan

    fp1, fp2 = np.array([]), np.array([])  # init.

    if fingerprint in [
        "biosynfoni",
        "overlap_binosynfoni",
        "maccsynfoni",
        "overlap",
    ]:
        fp1 = FP_FUNCTIONS[fingerprint](mol1, *args, **kwargs)
        fp2 = FP_FUNCTIONS[fingerprint](mol2, *args, **kwargs)
    else:
        fp1 = FP_FUNCTIONS[fingerprint](mol1)
        fp2 = FP_FUNCTIONS[fingerprint](mol2)

    return fp1, fp2


def similarity(
    fp1: np.array,
    fp2: np.array,
    metric: str = "c_tanimoto",
) -> float:
    """args and kwargs are passed to get_biosynfoni"""
    assert metric in SIM_FUNCTIONS.keys(), f"unsupported metric: {metric}"

    if not (isinstance(fp1, np.ndarray) and isinstance(fp1, np.ndarray)):
        return np.nan

    return _fps_to_distance(fp1, fp2, metric=metric)


# =============================================================================

# =============================== getting pairs ===============================


# def __get_step_pairs(chain_mols: list, degree_of_sep: int = 1) -> list:
#     """for given list of molecules originating from one reaction,
#     gets the available(non-null) molecule pairs separated by
#     degree of separation"""

#     pairs = []
#     for i in range(len(chain_mols) - degree_of_sep):
#         if chain_mols[i] and chain_mols[i + degree_of_sep]:
#             pairs.append([chain_mols[i], chain_mols[i + degree_of_sep]])
#     return pairs


# def __clean_structure_file(struct_loc: str) -> pd.DataFrame:
#     """yields rows from csv, cleaning up the rows in the process"""
#     struct = pd.read_csv(struct_loc, sep="\t", header=None)
#     # get headers
#     nums = ["{}".format(i - 1) for i in range(len(struct.columns))]
#     nums.pop(0)
#     headers = ["pathway"] + nums
#     struct.columns = headers
#     # clean up
#     struct.replace("", np.nan, inplace=True)  # just in case one was missed
#     filtered = struct.dropna(thresh=2)  # drop completely empty ones
#     final = filtered.replace(np.nan, "")  # replace back for later fnxs
#     return final


# def _dictify_pw_pairs(struct_loc: str, degree_of_sep: int = 1) -> dict[list[list[str]]]:
#     """returns {'pathway-id':[[product,precursor],[precursor,preprecursor]]}"""
#     dictionary = {}
#     struct_df = __clean_structure_file(struct_loc)
#     for _, row in tqdm(struct_df.iterrows(), desc="getting pairs"):
#         pathway_id = row[0]
#         pairs = __get_step_pairs(row[1:], degree_of_sep=degree_of_sep)
#         dictionary[pathway_id] = pairs
#     return dictionary


# def _add_pair_distances(
#     pairs_per_pathway: dict[list[list[str]]],
# ) -> tuple[str, list[str]]:
#     """returns list of (pathway,distances)"""
#     labeled_distances = []
#     for pw_id in pairs_per_pathway.keys():
#         for distance in pairs_per_pathway[pw_id]:
#             labeled_distances.append((pw_id, distance))
#     return labeled_distances


def _str_to_mol(repres: str, clean: bool = True) -> Chem.Mol:
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
        return np.nan


# def _get_pairs_dist_df(
#     labeled_pairs: list[tuple[str, list[str]]],
#     fingerprints: list[str] = FP_FUNCTIONS.keys(),
#     metric="c_tanimoto",
#     *args,
#     **kwargs,
# ) -> pd.DataFrame:
#     """gets the distances between molecules, for the given fingerprints
#     passes any other args and kwargs to biosynfoni-derived fp functions in similarity
#     only supports inchi_pairs at the moment"""
#     prec_strings = [x[1][0] for x in labeled_pairs]
#     prod_strings = [x[1][1] for x in labeled_pairs]
#     pathways = [x[0] for x in labeled_pairs]

#     df = pd.DataFrame(
#         {"pathway": pathways, "precursor": prec_strings, "product": prod_strings}
#     )

#     prec_mols = [_str_to_mol(x) for x in prec_strings]
#     prod_mols = [_str_to_mol(x) for x in prod_strings]
#     molpairs = [[prec_mols[i], prod_mols[i]] for i in range(len(prec_mols))]

#     # find way to convert to mol only once instead of continuously
#     for fp_type in fingerprints:
#         df[fp_type] = [
#             similarity(
#                 x[0],
#                 x[1],
#                 fingerprint=fp_type,
#                 metric=metric,
#                 check_moliness=False,
#                 *args,
#                 **kwargs,
#             )
#             for x in molpairs
#         ]
#     return df


# def get_comparison_df(struct_loc, max_steps=4, metric="c_tanimoto", *args, **kwargs):
#     df = pd.DataFrame()
#     for i in tqdm(range(max_steps), desc="getting comparison df"):
#         degree_of_sep = i + 1
#         previous_df = df.copy()
#         pw_pairs = _dictify_pw_pairs(struct_loc, degree_of_sep=degree_of_sep)
#         labeled_pairs = _add_pair_distances(pw_pairs)
#         # labeled_pairs = _label_pairs(dictionary)[:100]
#         df = _get_pairs_dist_df(labeled_pairs, metric=metric, *args, **kwargs)
#         df["reaction_steps"] = degree_of_sep
#         df = pd.concat([previous_df, df], axis=0)
#     return df


def _pairs_in_chain(chain_indexes: list[int], separation: int):
    # get all unique pairs separated by separation in chain_indexes
    pairs = []
    for i in range(len(chain_indexes) - separation):
        pairs.append([chain_indexes[i], chain_indexes[i + separation]])
    return pairs


def _explode_into_sep_i_pairs(df: pd.DataFrame, separation: int) -> pd.DataFrame:
    pair_indexes = _pairs_in_chain(df.columns, separation=separation)
    new_df = pd.DataFrame()
    # make df with pair_indexes and pws as rows, and the two mols as 2 columns
    for pair in pair_indexes:
        mol1 = df[pair[0]]
        mol2 = df[pair[1]]
        new_df = pd.concat([new_df, pd.DataFrame({"mol1": mol1, "mol2": mol2})], axis=0)
    # drop rows where either mol is null
    new_df = new_df.dropna()
    # sort by index
    new_df = new_df.sort_index()

    return new_df


def _get_random_pairs(df: pd.DataFrame, fraction: float = 0.25) -> pd.DataFrame:
    random_pairs = pd.DataFrame()
    for col in df.columns:
        random_pairs[col] = df[col].sample(frac=fraction).values
    # random_pairs = pd.DataFrame({"mol1":df["mol1"].sample(frac=fraction).values, "mol2":df["mol"].sample(frac=fraction).values})
    return random_pairs


def _turn_to_mol(df: pd.DataFrame, col: str) -> pd.DataFrame:
    df[col] = df[col].apply(_str_to_mol)
    return df


def pairs_per_separation(df: pd.DataFrame, max_separation: int = 4) -> pd.DataFrame:
    # get all pairs separated by separation in chain_indexes
    pair_dfs = pd.DataFrame()
    separations = []
    for i in range(max_separation):
        this_separation_df = _explode_into_sep_i_pairs(df, separation=i + 1)
        this_separation_df["separation"] = i + 1
        pair_dfs = pd.concat([pair_dfs, this_separation_df], axis=0)

    # add in random pairs, 1/4 of current pairs
    random_pairs = _get_random_pairs(pair_dfs, fraction=0.25)
    random_pairs["separation"] = -1
    pair_dfs = pd.concat([pair_dfs, random_pairs], axis=0).astype("category")
    # turn to mol
    pair_dfs = _turn_to_mol(pair_dfs, "mol1")
    pair_dfs = _turn_to_mol(pair_dfs, "mol2")
    # dropna
    pair_dfs = pair_dfs.dropna()
    # add unique indexes for each row
    # pair_dfs = pair_dfs.reset_index()
    return pair_dfs


def add_fp_to_df(
    df: pd.DataFrame, fp_types: list = FP_FUNCTIONS.keys(), *args, **kwargs
) -> pd.DataFrame:
    for fp_type in tqdm(fp_types, desc="getting fps"):
        df[fp_type] = df.apply(
            lambda x: _get_fp(
                x["mol1"],
                x["mol2"],
                fingerprint=fp_type,
                *args,
                **kwargs,
            ),
            axis=1,
        )
    return df


def add_similarities_per_fp(
    molpair_df: pd.DataFrame,
    fp_names: list[str] = FP_FUNCTIONS.keys(),
    metric: str = "c_tanimoto",
    *args,
    **kwargs,
) -> pd.DataFrame:
    """gets the distances between molecules, for the given fingerprints
    passes any other args and kwargs to biosynfoni-derived fp functions in similarity
    only supports inchi_pairs at the moment"""

    df = molpair_df.copy()

    for fp_type in tqdm(fp_names, desc="getting similarity"):
        df[fp_type] = df.apply(
            lambda x: similarity(
                x[fp_type][0],
                x[fp_type][1],
                metric=metric,
            ),
            axis=1,
        )
        df[fp_type] = df[fp_type].replace(-1, np.nan)
    return df


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
    logging.info(len(pw_tax), len(class_text))
    ndf["tax_id"] = ndf["pathway"].replace(to_replace=pw_tax)
    ndf["taxonomic_range"] = ndf["tax_id"].replace(to_replace=class_text)
    ndf["taxonomy"] = ndf["taxonomic_range"].replace(to_replace=BETTER_TAX)
    ndf["pathway_name"] = ndf["pathway"].replace(to_replace=class_text)
    # ndf['pathway_name']=ndf['pathway_name'].fillna(ndf['pathway'])
    return ndf


def get_fp_combinations() -> list[tuple[str, str]]:
    combinations = []
    for firstkey in FP_FUNCTIONS.keys():
        for secondkey in FP_FUNCTIONS.keys():
            if firstkey != secondkey:
                if (firstkey, secondkey) not in combinations:
                    if (secondkey, firstkey) not in combinations:
                        combinations.append((firstkey, secondkey))
    return combinations


def get_square(
    df: pd.DataFrame,
    col1: str,
    col2: str,
    range1: tuple[float, float],
    range2: tuple[float, float],
) -> pd.DataFrame:
    """returns the molecule pair's inchis for 'dots' in a given square of the scatter plot"""
    square = df.loc[
        (df[col1] >= range1[0])
        & (df[col1] <= range1[1])
        & (df[col2] >= range2[0])
        & (df[col2] <= range2[1])
    ]
    return square


def draw_molpair(
    pair: list[Chem.Mol], annotation: str = "", highlighting: bool = True
) -> None:
    for i in range(len(pair)):
        highlighting_info = None
        if highlighting:
            highlighting_info = get_highlight_mapping(mol=pair[i])
        svg_text = moldrawing.draw(
            pair[i], highlight_atoms_bonds_mappings=highlighting_info
        )
        if annotation:
            svg_text = svg_text.replace(
                "</svg>",
                f'<text x="30" y="30" font-size="20" font-family="montserrat">{annotation}</text></svg>',
            )
        with open(f"{annotation}_{i}.svg", "w") as f:
            f.write(svg_text)
    return None


def draw_squares(
    square_df: pd.DataFrame,
    pair_columns: tuple[str, str] = ("mol1", "mol2"),
    squarename: str = "origin",
    highlighting: bool = True,
) -> None:
    """draws the molecules in the squares
    input: (pd.DataFrame) square_df -- the dataframe containing the squares
    (str) pair_columns -- the name of the column containing the molecule pairs
    (str) squarename -- the name of the square
    """
    if square_df.empty:
        return None
    for _, row in tqdm(
        square_df.iterrows(),
        desc=f"drawing {squarename} squares",
        total=square_df.shape[0],
        position=1,
    ):
        pair = [row[pair_columns[0]], row[pair_columns[1]]]
        pathway = row["pathway"]
        outfilename = outfile_namer(f"{squarename}_{pathway}")
        draw_molpair(pair, annotation=outfilename, highlighting=highlighting)
    return None


def loopsquares(
    df: pd.DataFrame,
    x_fp: str = "biosynfoni",
    y_fps: list[str] = ["rdkit", "maccs", "morgan"],
    size: int = 0.2,
) -> None:
    for i, y_fp in tqdm(
        enumerate(y_fps), desc="looping squares", leave=False, position=0
    ):
        _, iwd = output_direr(f"./{x_fp}_{y_fp}_squares")
        min_val, max_val = 0.0, 1.0
        min_border = 0.0 + size
        max_border = 1.0 - size
        left_bottom = get_square(
            df,
            x_fp,
            y_fp,
            (min_val, min_border),
            (min_val, min_border),
        )
        left_top = get_square(
            df,
            x_fp,
            y_fp,
            (min_val, min_border),
            (max_border, max_val),
        )
        right_bottom = get_square(
            df,
            x_fp,
            y_fp,
            (max_border, max_val),
            (min_val, min_border),
        )
        right_top = get_square(
            df,
            x_fp,
            y_fp,
            (max_border, max_val),
            (max_border, max_val),
        )
        exactly_middle = get_square(
            df,
            x_fp,
            y_fp,
            (0.5, 0.5),
            (min_val, max_val),
        )
        draw_squares(left_bottom, squarename=f"{y_fp}_origin")
        draw_squares(left_top, squarename=f"{y_fp}_left_top")
        draw_squares(right_bottom, squarename=f"{y_fp}_right_bottom")
        draw_squares(right_top, squarename=f"{y_fp}_right_top")
        if i == 0:
            draw_squares(exactly_middle, squarename=f"{x_fp}_middle")
        os.chdir(iwd)
    return None


def main():
    # @TODO: incorporate the distances along the longest chain
    set_style()
    logging.getLogger(__name__).setLevel(logging.INFO)
    logging.info("==========\nbiosyn distance\n==========")
    # struct_loc = "../../../scripts/0927_metacyc_reactions.tsv"

    args = cli()

    df = pd.read_csv(args.structures_path, sep="\t", header=None, index_col=0)
    df.index.name = "pathway"
    logging.info(f"read {df.shape[0]} pathways from {args.structures_path}")

    old_n_rows = df.shape
    df.replace("", np.nan, inplace=True)
    df = df.dropna(thresh=2)  # drop where not at least 1 pair of mol
    logging.info(
        f"{old_n_rows[0]-df.shape[0]} pathways dropped due to lack of mol pairs\n\n"
    )

    _, iwd = output_direr("./biosynthetic_distance")  # move to outputdir

    pairs = pairs_per_separation(df)

    dist_df = f"{outfile_namer(args.metric)}.tsv"
    mols_pkl = f"{outfile_namer('mols')}.pkl"
    if not os.path.exists(dist_df) or not os.path.exists(mols_pkl):
        df = add_fp_to_df(
            pairs,
            fp_types=FP_FUNCTIONS.keys(),
        )
        df = add_similarities_per_fp(df, metric=args.metric)
        mols = df.copy()
        # save mols as pickle
        mols.to_pickle(f"{outfile_namer('mols')}.pkl")
        # remove mols from df
        df = df.drop(columns=["mol1", "mol2"])
        df.to_csv(dist_df, sep="\t", index=True)

    mols = pd.read_pickle(f"{outfile_namer('mols')}.pkl")
    df = mols
    df["pathway"] = df.index

    logging.debug(df.shape, mols.shape, df.columns)
    logging.debug(df, pairs, pairs.columns)

    # if args.annotate:
    #     # pw_tax_file = "../../../metacyc/pathways_taxid.txt"
    #     # tax_text_file = "../../../metacyc/cleaner_classes.dat"
    #     pw_tax_file, tax_text_file = args.annotate
    #     annotated_df = annotate_pathways(comparison_df, pw_tax_file, tax_text_file)

    df["stepnum"] = df["separation"].apply(str)

    logging.info("getting scatterplots...")
    fp_combs = get_fp_combinations()
    for combination in tqdm(fp_combs, desc="getting scatterplots"):
        scatter = fm.scatter_boxplots(
            df,
            col_x=combination[0],
            col_y=combination[1],
            figtitle=f"{args.metric} for different reaction step numbers",
            color_by="stepnum",
        )
        filename = outfile_namer(f"{combination[0]}_{combination[1]}_{args.metric}.png")
        fm.savefig(scatter, filename)

    onestep = df[df["stepnum"] == "1"]
    onestep.to_csv(f'{outfile_namer("onestep")}.tsv', sep="\t", index=False)
    logging.info("getting squares...")
    for fp in [
        "biosynfoni",
        "overlap_binosynfoni",
        "overlap_biosynfoni",
        "interoverlap_biosynfoni",
    ]:
        loopsquares(
            onestep,
            fp,
            ["rdkit", "maccs", "morgan", "maccsynfoni", "overlap"],
            size=0.2,
        )
    # loopsquares(
    #     onestep,
    #     "biosynfoni",
    #     ["rdkit", "maccs", "morgan", "maccsynfoni", "overlap"],
    #     size=0.2,
    # )

    # save the biosynfoni version for reference
    logging.info("saving current biosynfoni version...")
    save_version(defaultVersion)

    os.chdir(iwd)
    logging.info("done\nbyebye")
    exit(0)
    return None


if __name__ == "__main__":
    main()
