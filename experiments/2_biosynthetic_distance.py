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

from helper import ChangeDirectory
from biosynfoni import fingerprints as fp


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


def pairs_per_separation(chains: dict, max_separation: int) -> pd.DataFrame:
    """returns a dataframe with all pairs separated by separation (i.e. amount
    of biosynthetic steps between compounds) in chain_indexes"""
    pairs_records = []
    for sep in range(1, max_separation + 1):
        pairs = {
            pid: list(_generate_pairs(chain, sep)) for pid, chain in chains.items()
        }
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
        pairs_df[fp_name] = np.nan
        fp_dic = {k: v[fp_name] for k, v in id_to_fps.items()}
        df = pairs_df.copy().astype(str)
        df["fp1"] = df["mol1"].apply(lambda x: fp_dic[x] if x in fp_dic else np.nan)
        df["fp2"] = df["mol2"].apply(lambda x: fp_dic[x] if x in fp_dic else np.nan)
        similarities = [
            (
                SIM_FUNCTIONS[metric]([x[0], x[1]])
                if not isinstance(x[0], float) and not isinstance(x[1], float)
                else np.nan
            )
            for x in df[["fp1", "fp2"]].values
        ]
        pairs_df[fp_name] = similarities
        # pairs_df[fp_name] = np.where(similarities == -1, np.nan, similarities)
        # pairs_df[fp_name] = pairs_df[fp_name].replace(
        #     -1, np.nan
        # )
    # drop rows for which all fp_names are nan
    return pairs_df.dropna(subset=fp_names, how="all")
    # return pairs_df


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
            tsp_path = None
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

            # get shortest travelling salesman path
            tsp_path = nx.approximation.traveling_salesman_problem(
                nx.from_numpy_array(1 - sim_matrix), cycle=False
            )
            tsp_path = [random_order[i] for i in tsp_path]

        reconstructions.update(
            {
                f"{fp_name}_independent": path,
                f"{fp_name}_f_start": tipped_path_start,
                f"{fp_name}_f_end": tipped_path_end,
                f"{fp_name}_tsp": tsp_path,
            }
        )
    return reconstructions


def num_chains_per_length(pw_chains):
    lengths = [i for i in range(11)]
    num_chains = []
    for i in lengths:
        chain_length = i
        pw_chains = {k: v for k, v in pw_chains.items() if len(v) > chain_length}
        num_chains.append(len(pw_chains))
    num_chains_per_length = np.array(list(zip(lengths, num_chains)))
    np.savetxt(
        "num_chains_per_length.csv", num_chains_per_length, delimiter=",", fmt="%d"
    )


def main():
    input_dir = Path(sys.argv[1]).resolve(strict=True)  # root/data/input
    fp_dir = input_dir.parent.parent / "fps"  # root/fps
    chain_length = 6

    pathways = pd.read_csv(
        input_dir / "metacyc_pathways.tsv",
        sep="\t",
        index_col="pathway_id",
        header=0,
    )
    # compound_data = pd.read_csv(
    #     input_dir / "metacyc_compounds.tsv", sep="\t", index_col="compound_id", header=0
    # )

    pw_chains = generate_reaction_chains(pathways)

    with ChangeDirectory(input_dir.parent.parent / "output"):
        num_chains_per_length(pw_chains)
    pw_chains = {k: v for k, v in pw_chains.items() if len(v) > chain_length}
    pairs_df = pairs_per_separation(pw_chains, max_separation=chain_length)
    id_to_fps = id_fp_dict(
        input_dir / "metacyc.sdf", fp_dir, ["bsf", "rdk", "maccs", "morgan"]
    )
    similarities = all_similarity_scores(pairs_df, id_to_fps, metric="c_tanimoto")
    with ChangeDirectory(input_dir.parent.parent / "output"):
        pd.DataFrame.from_dict(pw_chains, orient="index").to_csv(
            "biosynthetic_chains.tsv", sep="\t"
        )
        similarities.to_csv("biosynthetic_distances.tsv", sep="\t", index=False)

        reconstructions = pd.DataFrame.from_records(
            [reconstruct_pathway(chain, id_to_fps) for chain in pw_chains.values()]
        )
        reconstructions["pathway"] = [pwid for pwid in pw_chains.keys()]
        reconstructions.set_index("pathway", inplace=True)
        reconstructions.to_csv("reconstructed_pathways.tsv", sep="\t", header=True)
    exit(0)


if __name__ == "__main__":
    main()
