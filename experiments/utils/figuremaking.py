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

# import plotly.express as px
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# def df_scatterplot(
#     df: pd.DataFrame,
#     col_x: str,
#     col_y: str,
#     figtitle: str,
#     filename: str = "scatterplot",
#     auto_open: bool = True,
#     *args,
#     **kwargs,
# ) -> None:
#     fig = px.scatter(df, x=col_x, y=col_y, *args, **kwargs)
#     fig.update_layout(title=figtitle)
#     fig.write_html(
#         f"{filename}.html",
#         auto_open=auto_open,
#     )
#     # fig.close()
#     return None


def _set_ax_boxplot_i_colour(ax_boxplot, i, colour, inner_alpha=0.6):
    translucent = mpl.colors.to_rgba(colour, inner_alpha)

    ax_boxplot["boxes"][i].set_facecolor(translucent)
    ax_boxplot["boxes"][i].set_edgecolor(colour)
    ax_boxplot["medians"][i].set_color(colour)
    ax_boxplot["whiskers"][i * 2].set_color(colour)
    ax_boxplot["whiskers"][i * 2 + 1].set_color(colour)
    ax_boxplot["caps"][i * 2].set_color(colour)
    ax_boxplot["caps"][i * 2 + 1].set_color(colour)
    ax_boxplot["fliers"][i].set_markeredgecolor(translucent)
    return ax_boxplot


def scatter_boxplots(
    df: pd.DataFrame,
    col_x: str,
    col_y: str,
    figtitle: str,
    color_by: str = "stepnum",
) -> plt.Figure:
    fig = plt.figure()
    # fig, ax = plt.subplots()
    # add gridspec for subplots

    gs = fig.add_gridspec(
        2,
        2,
        width_ratios=(4, 1),
        height_ratios=(1, 4),
        # left=0.1, right=0.9, bottom=0.1, top=0.9,
        wspace=-1,
        hspace=-5,
    )

    # Create the Axes.
    ax = fig.add_subplot(gs[1, 0])
    legax = fig.add_subplot(gs[0, 1])
    # ax.set_xlim(-0.05, 1.05)
    # ax.set_ylim(-0.05, 1.05)

    all_data_x = [
        np.array(df[df[color_by] == category][col_x].to_numpy(dtype=float))
        for category in df[color_by].unique()
    ]
    all_data_x = [x[~np.isnan(x)] for x in all_data_x]
    # all_data_x = np.random.normal(100, 10, 3426)
    all_data_y = [
        np.array(df[df[color_by] == category][col_y].tolist())
        for category in df[color_by].unique()
    ]
    all_data_y = [y[~np.isnan(y)] for y in all_data_y]
    ax_xobs, ax_yobs = [], []
    print(all_data_x)
    print(len(all_data_x))

    xax = fig.add_subplot(gs[0, 0], sharex=ax)
    yax = fig.add_subplot(gs[1, 1], sharey=ax)

    xax.tick_params(
        length=0, labelbottom=False, labelsize=5
    )  # , labelleft=False, labelbottom=False, labelsize=0)
    yax.tick_params(
        length=0, labelrotation=-30, labelleft=False, labelsize=5
    )  # , labelbottom=False, labelsize=0)
    legax.tick_params(length=0, labelleft=False, labelbottom=False, labelsize=0)
    ax.tick_params(length=0)  # , labelleft=False, labelbottom=False)

    labels = [
        f"{category}" if category != "-1" else "nonce"
        for category in df[color_by].unique()
    ]
    xplot = xax.boxplot(
        all_data_x, vert=False, patch_artist=True, labels=labels
    )  # labels=df[color_by].tolist())#, patch_artist=True
    yplot = yax.boxplot(
        all_data_y, vert=True, patch_artist=True, labels=labels
    )  # labels=df[color_by].tolist())#, patch_artist=True)
    print(xplot["boxes"])
    i = 0
    for category in df[color_by].unique():
        colour = COLOUR_DICT[color_by][category]

        label = f"{category}" if category != "-1" else "None"

        # ax.scatter(x, y, c=color, s=scale, label=color,alpha=0.3, edgecolors='none')
        scatterplot = ax.scatter(
            x=col_x,
            y=col_y,
            data=df[df[color_by] == category],
            c=colour,
            label=label,
            alpha=0.5,
            edgecolors="none",
            zorder=3,
        )

        alpha = 0.6

        _set_ax_boxplot_i_colour(xplot, i, colour, inner_alpha=alpha)
        _set_ax_boxplot_i_colour(yplot, i, colour, inner_alpha=alpha)
        leg = legax.scatter(
            x=col_x,
            y=col_y,
            data=df[df[color_by] == category][0:0],
            c=colour,
            label=label,
            alpha=0.5,
            edgecolors="none",
            s=10,
        )
        leg.set_facecolor(mpl.colors.to_rgba(colour, alpha=alpha))

        i += 1

    # ==================================================

    legax.legend(loc="lower left", prop={"size": 6}, frameon=False)
    squareside = 0.2
    s_color = "#7A7979AA"
    s_color = mpl.colors.to_rgba(COLOR, alpha=0.3)
    linewidth = 1

    # ax.set_xticklabels([0,0.2,0.4,0.6,0.8,1.0])
    ax.set_xlabel(col_x, labelpad=10)
    ax.set_ylabel(col_y, labelpad=10)
    # ax_xobs[0].set_title(figtitle, loc="center", pad=20)
    xax.set_title(figtitle, loc="center", pad=20)

    ax.grid(True, alpha=0.3, linewidth=0.5, mouseover=True)
    gs.tight_layout(fig)

    plt.show()
    # plt.savefig(f"{filename}.png", dpi=500)
    # fig.close()
    return fig


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
        "1": "#081d58",  #'#57BAC0',  # navy,
        "2": "#225ea8",  #'#77BC4D',  # royal blue,
        "3": "#41b6c4",  #'#F3C55F',  # teal,
        "4": "#7fcdbb",  #'#F48861',  # turquoise,
        "-1": "#c7e9b4",  #'#797979',  # lemon green,
        "random pairs": "#c7e9b4",
        "control": "#c7e9b4",
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


# def violins(df):
#     df = px.data.tips()
#     fig = px.violin(
#         df,
#         y="count",
#         x="fingerprint",
#         color="class",
#         box=True,
#         points="all",
#         hover_data=df.columns,
#     )
#     fig.show()

#     return None


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
