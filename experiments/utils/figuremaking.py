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
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
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


def heatmap(
    data: np.array,
    row_labels: list,
    col_labels: list,
    ax=None,
    cbar_kw=None,
    cbarlabel="",
    **kwargs,
):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (M, N).
    row_labels
        A list or array of length M with the labels for the rows.
    col_labels
        A list or array of length N with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """

    if ax is None:
        ax = plt.gca()

    if cbar_kw is None:
        cbar_kw = {}

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # Show all ticks and label them with the respective list entries.
    ax.set_xticks(np.arange(data.shape[1]), labels=col_labels)
    ax.set_yticks(np.arange(data.shape[0]), labels=row_labels)

    # Let the horizontal axes labeling appear on bottom (default)
    ax.tick_params(top=False, bottom=True, labeltop=False, labelbottom=True)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=-30, ha="right", rotation_mode="anchor")

    # Turn spines off and create white grid.
    ax.spines[:].set_visible(False)

    ax.set_xticks(np.arange(data.shape[1] + 1) - 0.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0] + 1) - 0.5, minor=True)
    ax.grid(which="minor", color="w", linestyle="-", linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar


def annotate_heatmap(
    im,
    data=None,
    valfmt="{x:.2f}",
    textcolors=("black", "white"),
    threshold=None,
    **textkw,
):
    """
    A function to annotate a heatmap.

    Parameters
    ----------
    im
        The AxesImage to be labeled.
    data
        Data used to annotate.  If None, the image's data is used.  Optional.
    valfmt
        The format of the annotations inside the heatmap.  This should either
        use the string format method, e.g. "$ {x:.2f}", or be a
        `matplotlib.ticker.Formatter`.  Optional.
    textcolors
        A pair of colors.  The first is used for values below a threshold,
        the second for those above.  Optional.
    threshold
        Value in data units according to which the colors from textcolors are
        applied.  If None (the default) uses the middle of the colormap as
        separation.  Optional.
    **kwargs
        All other arguments are forwarded to each call to `text` used to create
        the text labels.
    """

    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max()) / 2.0

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center", verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = mpl.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
            text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
            texts.append(text)

    return texts
