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


def df_scatterplot(
    df: pd.DataFrame,
    col_x: str,
    col_y: str,
    figtitle: str,
    filename: str = "scatterplot",
    auto_open: bool = True,
    *args,
    **kwargs,
) -> None:
    fig = px.scatter(df, x=col_x, y=col_y, *args, **kwargs)
    fig.update_layout(title=figtitle)
    fig.write_html(
        f"{filename}.html",
        auto_open=auto_open,
    )
    # fig.close()
    return None


COLOUR_DICT = {
    "taxonomy": {
        "Viridiplantae": px.colors.qualitative.Plotly[7],  # green
        "Bacteria": px.colors.qualitative.D3[9],  # light red
        "Fungi": px.colors.qualitative.Plotly[3],  # purple
        "Metazoa": px.colors.qualitative.Plotly[4],  # orange
        "Archaea": px.colors.qualitative.T10[7],  # pink
        "Eukaryota": px.colors.qualitative.Set3[8],
        "Cellular organisms": px.colors.qualitative.Plotly[9],
        "Opisthokonta": px.colors.qualitative.Plotly[8],
    },
    "stepnum": {
        "1": px.colors.qualitative.Pastel[0],  # blue,
        "2": px.colors.qualitative.Pastel[4],  # green,
        "3": px.colors.qualitative.Pastel[1],  # yellow,
        "4": px.colors.qualitative.Pastel[2],  # orange
    },
    "pathways": {
        "shikimate": px.colors.qualitative.Plotly[3],  # purple
        "acetate": px.colors.qualitative.Pastel[2],  # orange,
        "mevalonate": px.colors.qualitative.Pastel[4],  # green,
        "methylerythritol": px.colors.qualitative.Pastel[0],  # blue
        "sugar": px.colors.qualitative.T10[7],  # pink
    },
    "class": {
        "Terpenoids": px.colors.qualitative.Plotly[7],  # green
        "Alkaloids": "lightblue",  # purple
        "Shikimates and Phenylpropanoids": px.colors.qualitative.Plotly[3],
        "Fatty acids": px.colors.qualitative.Plotly[4],  # orange
        "Carbohydrates": "lightpink",  # pink
        "Polyketides": px.colors.qualitative.Prism[7],  # light red
        "Amino acids and Peptides": "bisque",
        "No NP-Classifier prediction": "grey",
        "None": "grey",
        "Synthetic": "black",
    },
}


def scatter_3d(df, col1, col2, col3, m="o"):
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")

    ax.scatter(df[col1], df[col2], df[col3], marker=m)

    ax.set_xlabel("Dim 1")
    ax.set_ylabel("Dim 2")
    ax.set_zlabel("Dim 3")

    plt.show()
    return None


def violins(df):
    df = px.data.tips()
    fig = px.violin(
        df,
        y="count",
        x="fingerprint",
        color="class",
        box=True,
        points="all",
        hover_data=df.columns,
    )
    fig.show()

    return None
