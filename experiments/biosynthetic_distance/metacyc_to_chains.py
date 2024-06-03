import sys
import subprocess
from sys import argv

import numpy as np
import pandas as pd

# sys.path.append("../src/")
from biosynfoni.inoutput import readr, entry_parser, per_entry


def extract_linestartswith(
    lines: list, starts: list = [], remove_starts: bool = True, strip_extra: str = " - "
) -> dict[list[str]]:
    """within a collection of lines, only extracts the lines starting with
    terms in starts, returning them as a list to accomodate multiple values
    per entry"""
    extraction = {}
    for start in starts:
        extraction[start] = []
    for line in lines:
        for start in starts:
            if line.startswith(start):
                extraction[start].append(
                    line.strip().replace(start, "").replace(strip_extra, "").strip()
                )
    return extraction


def get_all_pathways(
    info_loc: str,
    info2extract: list[str],
    ignore_start: str = "#",
    entry_sep="//",
    encoding="iso8859-1",
    remove_starts=True,
    strip_extra: str = " - ",
) -> pd.DataFrame:
    entries = entry_parser(
        readr(info_loc, ignore_start=ignore_start, encoding=encoding), sep=entry_sep
    )

    all_vals = per_entry(
        entries,
        extract_linestartswith,
        starts=info2extract,
        remove_starts=remove_starts,
        strip_extra=strip_extra,
    )

    # get data frame, will have doubles if rxnslists (because of 'expanding')
    # df = get_normalised_db(all_vals)
    df = pd.DataFrame.from_records(all_vals)

    return df


def get_reaction_df(df: pd.DataFrame) -> pd.DataFrame:
    """from df per pathway, gets df per reaction-layout. splits the reaction-list and reaction-layout columns into separate columns, removes waters and empty cells"""
    r_df = df.copy()
    # make all column headers lowercase
    r_df.columns = r_df.columns.str.lower()
    r_df.columns = r_df.columns.str.replace("-", "_")
    # separate row for each reaction
    r_df = r_df.explode("reaction_layout")

    r_df["reaction_id"] = (
        r_df["reaction_layout"].str.split(" ", expand=True)[0].str.strip("(")
    )
    r_df["left"] = (
        r_df["reaction_layout"]
        .str.split(":", expand=True)[1]
        .str.replace("LEFT-PRIMARIES", "")
        .str.strip(") (")
    )
    r_df["direction"] = (
        r_df["reaction_layout"].str.split(":", expand=True)[3].str.strip(") (")
    )
    r_df["right"] = (
        r_df["reaction_layout"]
        .str.split("RIGHT-PRIMARIES", expand=True)[1]
        .str.strip("))")
        .str.strip()
    )

    # turn empty cells into None
    r_df = r_df.replace(r"^\s*$", np.nan, regex=True)
    r_df.dropna(subset=["left", "right"], inplace=True)

    # remove water from left and right
    r_df["left"] = r_df["left"].str.replace("WATER", "")
    r_df["right"] = r_df["right"].str.replace("WATER", "")

    # turn empty cells into None
    r_df = r_df.replace(r"^\s*$", np.nan, regex=True)
    r_df.dropna(subset=["left", "right"], inplace=True)

    # get only pathways with 4 or more reactions in reaction-list
    r_df = r_df[r_df["reaction_list"].str.len() > 3]

    # convert pwd-id from list to string
    assert r_df[
        r_df["unique_id"].str.len() > 1
    ].empty, "there are pathways with more than one unique-id"

    r_df["unique_id"] = r_df["unique_id"].str[0]

    # rename column to pathway_id
    r_df.rename(columns={"unique_id": "pathway_id"}, inplace=True)
    return r_df[["reaction_id", "pathway_id", "left", "direction", "right"]]


def main():
    pathway_path = "/Users/lucina-may/thesis/metacyc/pathways.dat"
    info2extract = [
        "UNIQUE-ID",
        "REACTION-LIST",
        "SPECIES",
        "TAXONOMIC-RANGE",
        "REACTION-LAYOUT",
    ]
    pathways = get_all_pathways(pathway_path, info2extract)
    reactions = get_reaction_df(pathways)

    return None


if __name__ == "__main__":
    main()
