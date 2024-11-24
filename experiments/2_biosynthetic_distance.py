#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
from pathlib import Path
from typing import Generator


import networkx as nx  # tried with 3.3, now 2.2 --> @TODO: check if it works with 2.2
import numpy as np
import pandas as pd
from rdkit import Chem
from tqdm import tqdm

# from biosynfoni.inoutput import *
from biosynfoni import fingerprints as fp  # for similarity functions
from helper import ChangeDirectory


def compound_graph(reaction_info: pd.DataFrame) -> nx.DiGraph:
    """Returns graph of reagents per MetaCyc pathway ID"""
    compound_graph = nx.DiGraph()
    for _, row in reaction_info.iterrows():
        compound_graph.add_edge(
            row["left"], row["right"], reaction_id=row["reaction_id"]
        )
    return compound_graph


def longest_path(pathway_graph) -> list:
    if not nx.is_directed_acyclic_graph(pathway_graph):
        longest_path = []
        for node in pathway_graph.nodes:
            for path in nx.all_simple_paths(
                pathway_graph, source=node, target=pathway_graph.nodes
            ):
                if len(path) > len(longest_path):
                    longest_path = path
        return longest_path
    else:
        return nx.dag_longest_path(pathway_graph)  # type: ignore


def generate_reaction_chains(pathways: pd.DataFrame) -> dict:
    longest_paths = {
        pathway_id: longest_path(compound_graph(reaction_info))
        for pathway_id, reaction_info in pathways.groupby("pathway_id")
    }
    return longest_paths
    # return pd.DataFrame.from_dict(longest_paths, orient="index")


# =========================== GLOBALS =========================


SIM_FUNCTIONS = {
    "c_tanimoto": fp.countanimoto,
    "cosine_sim": fp.cosine_sim,
    "manhattan_sim": fp.manhattan_sim,
    "tanimoto_sim": fp.tanimoto_sim,
    "euclidean_sim": fp.euclidean_sim,
}

# =============================================================================

def id_mol_dict(struct_loc: str) -> dict:
    return {mol.GetProp("compound_id"): mol for mol in Chem.SDMolSupplier(struct_loc)}  # type: ignore


def id_fp_dict(struct_loc: str, fp_loc: str, fp_names: list) -> dict:
    fp_files = [
        np.loadtxt(fp_loc / f"metacyc_{name}.csv", dtype=int, delimiter=",")  # type: ignore
        for name in fp_names
    ]
    ids = id_mol_dict(struct_loc).keys()
    return {
        id_: {fp_name: fp for fp_name, fp in zip(fp_names, fps)}
        for id_, *fps in zip(ids, *fp_files)
    }


def _generate_pairs(chain, separation) -> Generator:
    for i in range(len(chain) - separation):
        yield (chain[i], chain[i + separation])


def _get_random_pairs(df: pd.DataFrame, fraction: float = 0.25) -> pd.DataFrame:
    random_pairs = pd.DataFrame()
    for col in df.columns:
        random_pairs[col] = df[col].sample(frac=fraction).values
    return random_pairs


def pairs_per_separation(chains: dict, max_separation: int = 4) -> pd.DataFrame:
    """returns a dataframe with all pairs separated by separation (i.e. amount
    of biosynthetic steps between compounds) in chain_indexes"""
    pairs_records = []
    for sep in range(1, max_separation + 1):
        pairs = {pid: _generate_pairs(chain, sep) for pid, chain in chains.items()}
        pairs_records.extend(
            [
                {"pathway_id": pid, "mol1": pair[0], "mol2": pair[1], "separation": sep}
                for pid, pair_gen in pairs.items()
                for pair in pair_gen
            ]
        )
    pairs_records = pd.DataFrame.from_records(pairs_records)

    # add in random pairs from molecules, in equal amounts of
    # rows for of each separation, i.e. if 4 separations, fraction is 4
    random_pairs = _get_random_pairs(pairs_records, fraction=1 / max_separation)
    random_pairs["separation"] = -1
    # as category to reduce the time to map fingerprints
    return pd.concat([pairs_records, random_pairs], axis=0).astype("category")


def all_similarity_scores(
    pairs_df: pd.DataFrame,
    id_to_fps: dict,
    metric: str = "c_tanimoto",
) -> pd.DataFrame:
    """gets the distances between molecules, for the given fingerprints
    passes any other args and kwargs to biosynfoni-derived fp functions in similarity
    only supports inchi_pairs at the moment"""

    fp_names = list(id_to_fps.values())[0].keys()
    for fp_name in tqdm(fp_names, desc="getting similarity"):
        fp_dic = {k: v[fp_name] for k, v in id_to_fps.items()}
        df = pairs_df.copy().astype(str)
        df["fp1"] = df["mol1"].apply(lambda x: fp_dic[x] if x in fp_dic else None)
        df["fp2"] = df["mol2"].apply(lambda x: fp_dic[x] if x in fp_dic else None)
        # df = df.dropna()
        # pairs_df = pairs_df.loc[df.index]
        df_none = df.isna().sum(axis=1) > 0
        df = df[~df_none]
        pairs_df = pairs_df[~df_none]
        pairs_df[fp_name] = [
            SIM_FUNCTIONS[metric]([x[0], x[1]]) for x in df[["fp1", "fp2"]].values
        ]
        pairs_df[fp_name].replace(-1, np.nan, inplace=True)
    return pairs_df


def _chain_from_sim_matrix(sim_matrix, start=None):
    if start is None:
        start = np.argmin(sim_matrix.sum(axis=1))
    return np.argsort(sim_matrix[start])[::-1]


def reconstruct_pathway(pathway_chain, id_to_fp):
    right_order = pathway_chain
    random_order = np.random.permutation(pathway_chain)
    reconstructions = {"true": right_order}
    invalid_pathway = False
    if not all([id_ in id_to_fp for id_ in pathway_chain]):
        invalid_pathway = True
    for fp_name in list(id_to_fp.values())[0].keys():
        if invalid_pathway:
            path = None
            tipped_path_start = None
            tipped_path_end = None
        else:
            fps = [id_to_fp[id_][fp_name] for id_ in random_order]

            sim_matrix = np.zeros((len(fps), len(fps)))
            # fill the matrix in one vectorized operation
            for i in range(len(fps)):
                for j in range(len(fps)):
                    sim_matrix[i, j] = sim_matrix[j, i] = SIM_FUNCTIONS["c_tanimoto"](
                        [fps[i], fps[j]]
                    )
            # sim_matrix = sim_matrix + sim_matrix.T - np.diag(sim_matrix.diagonal())
            # np.fill_diagonal(sim_matrix, 1)

            path = [
                random_order[i] for i in _chain_from_sim_matrix(sim_matrix)
            ]  # as it includes itself
            # get index for the item in random_order that is first in right_order
            start_index = np.where(random_order == right_order[0])[0][0]
            end_index = np.where(random_order == right_order[-1])[0][0]
            tipped_path_start = [
                random_order[i]
                for i in _chain_from_sim_matrix(sim_matrix, start=start_index)
            ]
            tipped_path_end = [
                random_order[i]
                for i in _chain_from_sim_matrix(sim_matrix, start=end_index)
            ][::-1]

        reconstructions.update(
            {
                f"{fp_name}_independent": path,
                f"{fp_name}_f_start": tipped_path_start,
                f"{fp_name}_f_end": tipped_path_end,
            }
        )
    return reconstructions


def main():
    input_dir = Path(sys.argv[1]).resolve(strict=True)  # root/data/input
    fp_dir = input_dir.parent.parent / "fps"  # root/fps
    chain_length = 8

    pathways = pd.read_csv(
        input_dir / "metacyc_pathways.tsv",
        sep="\t",
        index_col="pathway_id",
        header=0,
    )
    compound_data = pd.read_csv(
        input_dir / "metacyc_compounds.tsv", sep="\t", index_col="compound_id", header=0
    )

    pw_chains = generate_reaction_chains(pathways)
    pw_chains = {k: v for k, v in pw_chains.items() if len(v) > chain_length}
    pairs_df = pairs_per_separation(pw_chains, max_separation=chain_length)
    # id_to_mol = id_mol_dict(input_dir / "metacyc.sdf")
    id_to_fps = id_fp_dict(
        input_dir / "metacyc.sdf", fp_dir, ["bsf", "rdk", "maccs", "morgan"]
    )
    similarities = all_similarity_scores(pairs_df, id_to_fps, metric="c_tanimoto")
    with ChangeDirectory(input_dir.parent.parent / "output"):
        pd.DataFrame.from_dict(pw_chains, orient="index").to_csv(
            "biosynthetic_chains.tsv", sep="\t"
        )
        similarities.to_csv("biosynthetic_distances.tsv", sep="\t", index=False)

        # reconstructions:
        reconstructions = pd.DataFrame.from_records(
            [reconstruct_pathway(chain, id_to_fps) for chain in pw_chains.values()]
        )
        reconstructions["pathway"] = [pwid for pwid in pw_chains.keys()]
        reconstructions.set_index("pathway", inplace=True)
        reconstructions.to_csv("reconstructed_pathways.tsv", sep="\t", header=True)
    exit(0)


if __name__ == "__main__":
    main()


# # later for visualisation:

# def get_square(
#     df: pd.DataFrame,
#     col1: str,
#     col2: str,
#     range1: tuple[float, float],
#     range2: tuple[float, float],
# ) -> pd.DataFrame:
#     """returns the molecule pair's inchis for 'dots' in a given square of the scatter plot"""
#     square = df.loc[
#         (df[col1] >= range1[0])
#         & (df[col1] <= range1[1])
#         & (df[col2] >= range2[0])
#         & (df[col2] <= range2[1])
#     ]
#     return square


# def draw_molpair(
#     pair: list[Chem.Mol], annotation: str = "", highlighting: bool = True
# ) -> None:
#     for i in range(len(pair)):
#         highlighting_info = None
#         if highlighting:
#             highlighting_info = get_highlight_mapping(mol=pair[i])
#         svg_text = moldrawing.draw(
#             pair[i], highlight_atoms_bonds_mappings=highlighting_info
#         )
#         if annotation:
#             svg_text = svg_text.replace(
#                 "</svg>",
#                 f'<text x="30" y="30" font-size="20" font-family="montserrat">{annotation}</text></svg>',
#             )
#         with open(f"{annotation}_{i}.svg", "w") as f:
#             f.write(svg_text)
#     return None


# def draw_squares(
#     square_df: pd.DataFrame,
#     pair_columns: tuple[str, str] = ("mol1", "mol2"),
#     squarename: str = "origin",
#     highlighting: bool = True,
# ) -> None:
#     """draws the molecules in the squares
#     input: (pd.DataFrame) square_df -- the dataframe containing the squares
#     (str) pair_columns -- the name of the column containing the molecule pairs
#     (str) squarename -- the name of the square
#     """
#     if square_df.empty:
#         return None
#     for _, row in tqdm(
#         square_df.iterrows(),
#         desc=f"drawing {squarename} squares",
#         total=square_df.shape[0],
#         position=1,
#     ):
#         pair = [row[pair_columns[0]], row[pair_columns[1]]]
#         pathway = row["pathway"]
#         outfilename = outfile_namer(f"{squarename}_{pathway}")
#         draw_molpair(pair, annotation=outfilename, highlighting=highlighting)
#     return None


# def loopsquares(
#     df: pd.DataFrame,
#     x_fp: str = "biosynfoni",
#     y_fps: list[str] = ["rdkit", "maccs", "morgan"],
#     size: int = 0.2,
# ) -> None:
#     for i, y_fp in tqdm(
#         enumerate(y_fps), desc="looping squares", leave=False, position=0
#     ):
#         _, iwd = output_direr(f"./{x_fp}_{y_fp}_squares")
#         min_val, max_val = 0.0, 1.0
#         min_border = 0.0 + size
#         max_border = 1.0 - size
#         left_bottom = get_square(
#             df,
#             x_fp,
#             y_fp,
#             (min_val, min_border),
#             (min_val, min_border),
#         )
#         left_top = get_square(
#             df,
#             x_fp,
#             y_fp,
#             (min_val, min_border),
#             (max_border, max_val),
#         )
#         right_bottom = get_square(
#             df,
#             x_fp,
#             y_fp,
#             (max_border, max_val),
#             (min_val, min_border),
#         )
#         right_top = get_square(
#             df,
#             x_fp,
#             y_fp,
#             (max_border, max_val),
#             (max_border, max_val),
#         )
#         exactly_middle = get_square(
#             df,
#             x_fp,
#             y_fp,
#             (0.5, 0.5),
#             (min_val, max_val),
#         )
#         draw_squares(left_bottom, squarename=f"{y_fp}_origin")
#         draw_squares(left_top, squarename=f"{y_fp}_left_top")
#         draw_squares(right_bottom, squarename=f"{y_fp}_right_bottom")
#         draw_squares(right_top, squarename=f"{y_fp}_right_top")
#         if i == 0:
#             draw_squares(exactly_middle, squarename=f"{x_fp}_middle")
#         os.chdir(iwd)
#     return None


# def biosynthetic_distance_analysis(pairs_df, metric):
#     df = pd.read_csv(structures_path, sep="\t", header=None, index_col=0)
#     logging.info(
#         f"{old_n_rows[0]-df.shape[0]} pathways dropped due to lack of mol pairs\n\n"
#     )

#     _, iwd = output_direr("./biosynthetic_distance")  # move to outputdir

#     pairs = pairs_per_separation(df)

#     dist_df = f"{outfile_namer(metric)}.tsv"
#     mols_pkl = f"{outfile_namer('mols')}.pkl"
#     if not os.path.exists(dist_df) or not os.path.exists(mols_pkl):
#         df = add_fp_to_df(
#             pairs,
#             fp_types=FP_FUNCTIONS.keys(),
#         )
#         df = get_all_similarity_scores(df, metric=metric)
#         mols = df.copy()
#         # save mols as pickle
#         mols.to_pickle(f"{outfile_namer('mols')}.pkl")
#         # remove mols from df
#         df = df.drop(columns=["mol1", "mol2"])
#         df.to_csv(dist_df, sep="\t", index=True)

#     mols = pd.read_pickle(f"{outfile_namer('mols')}.pkl")
#     df = mols
#     df["pathway"] = df.index

#     logging.debug(df.shape, mols.shape, df.columns)
#     logging.debug(df, pairs, pairs.columns)

#     # if args.annotate:
#     #     # pw_tax_file = "../../../metacyc/pathways_taxid.txt"
#     #     # tax_text_file = "../../../metacyc/cleaner_classes.dat"
#     #     pw_tax_file, tax_text_file = args.annotate
#     #     annotated_df = annotate_pathways(comparison_df, pw_tax_file, tax_text_file)

#     df["stepnum"] = df["separation"].apply(str)

#     logging.info("getting scatterplots...")
#     fp_combs = list(itertools.combinations(fp_names, 2))
#     for combination in tqdm(fp_combs, desc="getting scatterplots"):
#         scatter = fm.scatter_boxplots(
#             df,
#             col_x=combination[0],
#             col_y=combination[1],
#             figtitle=f"{args.metric} for different reaction step numbers",
#             color_by="stepnum",
#         )
#         filename = outfile_namer(f"{combination[0]}_{combination[1]}_{args.metric}.png")
#         fm.savefig(scatter, filename)

#     onestep = df[df["stepnum"] == "1"]
#     onestep.to_csv(f'{outfile_namer("onestep")}.tsv', sep="\t", index=False)
#     logging.info("getting squares...")
#     for fp in [
#         "biosynfoni",
#         "overlap_binosynfoni",
#         "overlap_biosynfoni",
#         "interoverlap_biosynfoni",
#     ]:
#         loopsquares(
#             onestep,
#             fp,
#             ["rdkit", "maccs", "morgan", "maccsynfoni", "overlap"],
#             size=0.2,
#         )
#     # loopsquares(
#     #     onestep,
#     #     "biosynfoni",
#     #     ["rdkit", "maccs", "morgan", "maccsynfoni", "overlap"],
#     #     size=0.2,
#     # )

#     # save the biosynfoni version for reference
#     logging.info("saving current biosynfoni version...")
#     save_version(defaultVersion)

#     os.chdir(iwd)
#     logging.info("done\nbyebye")
#     exit(0)
#     return None
