import logging

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
import numpy as np
import pandas as pd


colourDict = {
    "taxonomy": {
        "Viridiplantae": "#B9C311",  # green
        "Bacteria": "#9BC2BA",  # teal
        "Fungi": "#CFC0D6",  # purple
        "Metazoa": "#FFAD61",  # light orange
        "Archaea": "#EB6737",  # soft red
        "Eukaryota": "#D2CEC4",  # pale grey
        "Cellular organisms": "#FFEAA0",  # yellow
        "Opisthokonta": "#FFC4CE",  # pink
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
        "shikimate": "#A783B6",  # purple
        "acetate": "#FF8B61",  # orange,
        "mevalonate": "#B9C311",  # green,
        "methylerythritol": "#6FB5C6",  # blue
        "sugar": "#FFC4CE",  # pink
        "amino": "#FFEAA0",  # yellow
        "amino_acid": "#FFEAA0",  # yellow
    },
    "class": {
        # "Terpenoids": "#9BC2BA",  # soft bluegreen
        "Terpenoids": "#B9C311",  # green
        "Alkaloids": "#B4CAD8",  # purple
        "Shikimates and Phenylpropanoids": "#A783B6",  # purple
        "Fatty acids": "#FF8B61",  # orange
        "Carbohydrates": "#FFC4CE",  # pink
        "Polyketides": "#C21100",  # soft red
        "Amino acids and Peptides": "#FFEAA0",  # yellow
        "No NP-Classifier prediction": "#797979",
        "None": "#595959",
        "Synthetic": "#393939",
        "Multiple": "#BBBBBB",
        # lowercase
        # "terpenoids": "#9BC2BA",  # soft bluegreen
        "terpenoids": "#B9C311",  # green
        "alkaloids": "#B4CAD8",  # purple
        "shikimates and phenylpropanoids": "#A783B6",  # purple
        "fatty acids": "#FF8B61",  # orange
        "carbohydrates": "#FFC4CE",  # pink
        "polyketides": "#C21100",  # soft red
        "amino acids and peptides": "#FFEAA0",  # yellow
        # chebi:
        "phenylpropanoid": "#A783B6",  # purple
        "fatty_acid": "#FF8B61",  # orange
        "polyketide": "#C21100",  # soft red
        "alkaloid": "#B4CAD8",  # purple
        # "isoprenoid": "#9BC2BA",  # soft bluegreen
        "isoprenoid": "#B9C311",  # green
        "carbohydrate": "#FFC4CE",  # pink
        "amino_acid": "#FFEAA0",  # yellow
        "synthetic": "#393939",  # grey
    },
    "NPClassifier prediction": {
        "Terpenoids": "#B9C311",  # green
        "Alkaloids": "#B4CAD8",  # purple
        "Shikimates and Phenylpropanoids": "#A783B6",  # purple
        "Fatty acids": "#FF8B61",  # orange
        "Carbohydrates": "#FFC4CE",  # pink
        "Polyketides": "#C21100",  # soft red
        "Amino acids and Peptides": "#FFEAA0",  # yellow
        "No NP-Classifier prediction": "#797979",
        "None": "#595959",
        "Synthetic": "#393939",
        "Multiple": "#BBBBBB",
    },
    "chebi class": {
        "phenylpropanoid": "#A783B6",  # purple
        "fatty_acid": "#FF8B61",  # orange
        "polyketide": "#C21100",  # soft red
        "alkaloid": "#B4CAD8",  # purple
        # "isoprenoid": "#9BC2BA",  # soft bluegreen
        "isoprenoid": "#B9C311",  # green
        "carbohydrate": "#FFC4CE",  # pink
        "amino_acid": "#FFEAA0",  # yellow
        "synthetic": "#393939",  # grey
    },
}


def savefig(fig, filename):
    """
    Save a figure to a file

        Args:
            fig (matplotlib.figure.Figure): figure to save
            filename (str): filename to save to

        Returns:
            None
    """
    if not "." in filename:
        filename = f"{filename}.png"
    fig.savefig(filename, dpi=500)
    fig.savefig(f"{filename.split('.')[0]}.png", dpi=500)
    plt.close(fig)
    return None


def custom_cmap(default_cmap, last_color=None, first_color=None):
    """
    Edit a colormap to have a specific first and last color

        Args:
            default_cmap (str): name of the colormap to edit
            last_color (bool): if True, the last color is set to white
            first_color (bool): if True, the first color is set to white

        Returns:
            matplotlib.colors.Colormap: the edited colormap

    """
    if not isinstance(default_cmap, mpl.colors.Colormap):
        cmap = mpl.colormaps[default_cmap]
    else:
        cmap = default_cmap
    # extract all colors from the .jet map
    cmaplist = [cmap(i) for i in range(cmap.N)]
    # force the first color entry to be grey
    if last_color:
        cmaplist[-1] = (1.0, 1.0, 0.95, 0)
    if first_color:
        cmaplist[0] = (0.95, 0.95, 0.95, 0)

    # create the new map
    cmap = mpl.colors.LinearSegmentedColormap.from_list("Custom cmap", cmaplist, cmap.N)
    return cmap


def cleanfmt(text):
    """
    Clean a string or list of strings to be used as labels in a plot

        Args:
            text (str or list): text to clean

        Returns:
            str or list: cleaned text

    Remarks:
        - replaces underscores with spaces
        - makes all text lowercase
    """
    if isinstance(text, str):
        return text.replace("_", " ").lower()
    elif isinstance(text, list):
        newtext = []
        for t in text:
            if isinstance(t, str):
                newtext.append(t.replace("_", " ").lower())
            else:
                newtext.append(t)
        return newtext
    else:
        return text


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

        Args:
            data (np.array): 2D array of data for heatmap
            row_labels (list): list of labels for the rows
            col_labels (list): list of labels for the columns
            ax (matplotlib.axes.Axes): a `matplotlib.axes.Axes` object, optional
            cbar_kw (dict): a dictionary with arguments to pass to `cbar` when creating the colorbar, optional
            cbarlabel (str): label for the colorbar, optional
            **kwargs: other arguments are forwarded to `imshow`

        Returns:
            matplotlib.image.AxesImage: the heatmap
            matplotlib.colorbar.Colorbar: the colorbar

    Remarks:
        - adapted from https://matplotlib.org/stable/gallery/images_contours_and_fields/image_annotated_heatmap.html

    reate a heatmap from a numpy array and two lists of labels.
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
    # cbar.ax.set_yticks([]) #remove tick labels
    cbar.outline.set_visible(False)

    # Show all ticks and label them with the respective list entries.
    ax.set_xticks(np.arange(data.shape[1]), labels=cleanfmt(col_labels), size=6)
    ax.set_yticks(np.arange(data.shape[0]), labels=cleanfmt(row_labels), size=6)
    # ax.set_yticks([])

    # Let the horizontal axes labeling appear on bottom (default)
    ax.tick_params(top=False, bottom=True, labeltop=False, labelbottom=True)
    ax.tick_params(pad=0.5)

    # Rotate the tick labels and set their alignment.
    plt.setp(
        ax.get_xticklabels(),
        rotation=90,
        ha="right",
        va="center",
        rotation_mode="anchor",
    )

    # Turn spines off and create white grid.
    ax.spines[:].set_visible(False)

    ax.set_xticks(np.arange(data.shape[1] + 1) - 0.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0] + 1) - 0.5, minor=True)
    ax.grid(which="minor", color="w", linestyle="-", linewidth=1, alpha=1.0)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar


def annotate_heatmap(
    im,
    data=None,
    valfmt="{x:.2f}",
    textcolors=("black", "white"),
    threshold=None,
    shift=0,
    **textkw,
):
    """
    Annotate a heatmap.

        Args:
            im (matplotlib.image.AxesImage): the heatmap to be labeled
            data (np.array): data used to annotate, optional
            valfmt (str): the format of the annotations inside the heatmap, optional
            textcolors (tuple): a pair of colors. The first is used for values below a threshold, the second for those above, optional
            threshold (float): value in data units according to which the colors from textcolors are applied. If None (the default) uses the middle of the colormap as separation, optional
            shift (int): shift the annotations by a certain amount, optional
            **kwargs: all other arguments are forwarded to each call to `text` used to create the text labels

        Returns:
            list: the text labels

    Remarks:
        - adapted from https://matplotlib.org/stable/gallery/images_contours_and_fields/image_annotated_heatmap.html

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


def _set_ax_boxplot_i_colour(
    ax_boxplot: mpl.container.BarContainer,
    i: int,
    colour: str,
    inner_alpha: float = 0.6,
):
    """
    Set the colour of a boxplot element

        Args:
            ax_boxplot (matplotlib.container.BarContainer): the boxplot to change
            i (int): the index of the element to change
            colour (str): the colour to change to
            inner_alpha (float): the alpha of the inner colour, optional. Default is 0.6

        Returns:
            matplotlib.container.BarContainer: the changed boxplot

    """
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
    *args,
    **kwargs,
) -> plt.Figure:
    """
    Make a scatterplot with boxplots on the axes

        Args:
            df (pd.DataFrame): dataframe to plot
            col_x (str): column to plot on the x-axis
            col_y (str): column to plot on the y-axis
            figtitle (str): title of the figure
            color_by (str): column to colour by, optional. Default is "stepnum"
            *args: other arguments to pass to scatterplot
            **kwargs: other keyword arguments to pass to scatterplot
        Returns:
            plt.Figure: the figure

    """
    fig = plt.figure()

    gs = fig.add_gridspec(
        2,
        2,
        width_ratios=(4, 1),
        height_ratios=(1, 4),
        # left=0.1, right=0.9, bottom=0.1, top=0.9,
        wspace=-1,
        # hspace=-5,
    )

    # Create the Axes.
    # for scatterplot
    sc_ax = fig.add_subplot(gs[1, 0])
    # for legend
    legax = fig.add_subplot(gs[0, 1])

    # Get Data
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
    logging.debug(all_data_x)
    logging.debug(len(all_data_x))

    top_bp_ax = fig.add_subplot(gs[0, 0], sharex=sc_ax)
    right_bp_ax = fig.add_subplot(gs[1, 1], sharey=sc_ax)

    top_bp_ax.tick_params(
        length=0, labelbottom=False, labelsize=5
    )  # , labelleft=False, labelbottom=False, labelsize=0)
    right_bp_ax.tick_params(
        length=0, labelrotation=-30, labelleft=False, labelsize=5
    )  # , labelbottom=False, labelsize=0)
    legax.tick_params(length=0, labelleft=False, labelbottom=False, labelsize=0)
    sc_ax.tick_params(length=0)  # , labelleft=False, labelbottom=False)

    labels = [
        f"{category}" if category != "-1" else "control"
        for category in df[color_by].unique()
    ]
    xplot = top_bp_ax.boxplot(
        all_data_x, vert=False, patch_artist=True, labels=labels
    )  # labels=df[color_by].tolist())#, patch_artist=True
    yplot = right_bp_ax.boxplot(
        all_data_y, vert=True, patch_artist=True, labels=labels
    )  # labels=df[color_by].tolist())#, patch_artist=True)
    logging.info(xplot["boxes"])
    i = 0
    for category in df[color_by].unique():
        colour = colourDict[color_by][category]

        label = f"{category}" if category != "-1" else "control"

        # ax.scatter(x, y, c=color, s=scale, label=color,alpha=0.3, edgecolors='none')
        scatterplot = sc_ax.scatter(
            x=col_x,
            y=col_y,
            data=df[df[color_by] == category],
            c=colour,
            label=label,
            alpha=0.5,
            edgecolors="none",
            zorder=3,
            *args,
            **kwargs,
        )

        # change boxplot colours to match the category colour
        alpha = 0.6
        _set_ax_boxplot_i_colour(xplot, i, colour, inner_alpha=alpha)
        _set_ax_boxplot_i_colour(yplot, i, colour, inner_alpha=alpha)

        # scatter empty df, to get legend in right format in right position
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

    sc_ax.set_xlabel(cleanfmt(col_x), labelpad=10)
    sc_ax.set_ylabel(cleanfmt(col_y), labelpad=10)
    top_bp_ax.set_title(cleanfmt(figtitle), loc="center", pad=20)

    sc_ax.grid(True, alpha=0.3, linewidth=0.5, mouseover=True)
    gs.tight_layout(fig)

    return fig


def scatter(
    df: pd.DataFrame,
    col_x: str,
    col_y: str,
    figtitle: str,
    color_by: str = "stepnum",
    *args,
    **kwargs,
) -> plt.Figure:
    """
    Make a scatterplot

        Args:
            df (pd.DataFrame): dataframe to plot
            col_x (str): column to plot on the x-axis
            col_y (str): column to plot on the y-axis
            figtitle (str): title of the figure
            color_by (str): column to colour by, optional. Default is "stepnum"
            *args: other arguments to pass to scatterplot
            **kwargs: other keyword arguments to pass to scatterplot
        Returns:
            plt.Figure: the figure
    """
    # check if figsize is given in kwargs, if not, set default figsize
    if "figsize" not in kwargs:
        kwargs["figsize"] = (2, 2)
    fig = plt.figure(figsize=kwargs["figsize"])
    # remove figsize from kwargs, so it doesn't get passed to scatterplot
    kwargs.pop("figsize", None)
    fig, ax = plt.subplots()

    # Get Data
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

    # ax_xobs, ax_yobs = [], []
    logging.debug(all_data_x)
    logging.debug(len(all_data_x))

    ax.tick_params(length=0)  # , labelleft=False, labelbottom=False)

    labels = [
        f"{category}" if category != "-1" else "control"
        for category in df[color_by].unique()
    ]

    for i, category in enumerate(df[color_by].unique()):
        colour = colourDict[color_by][category]
        label = f"{category}" if category != "-1" else "control"

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
            *args,
            **kwargs,
        )

    # ==================================================

    ax.legend(loc="lower left", prop={"size": 6}, frameon=True, edgecolor="#FFFFFFAA")

    ax.set_xlabel(cleanfmt(col_x), labelpad=10)
    ax.set_ylabel(cleanfmt(col_y), labelpad=10)
    ax.set_title(cleanfmt(figtitle), loc="center", pad=20)

    ax.grid(True, alpha=0.3, linewidth=0.5, mouseover=True)
    fig.tight_layout()

    return fig


def scatter_3d(df, col1, col2, col3, m="o"):
    """
    Make a 3D scatterplot

        Args:
            df (pd.DataFrame): dataframe to plot
            col1 (str): column to plot on the x-axis
            col2 (str): column to plot on the y-axis
            col3 (str): column to plot on the z-axis
            m (str): marker to use, optional. Default is "o"
        Returns:
            None
    Remarks:
        - under construction
    """
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")

    ax.scatter(df[col1], df[col2], df[col3], marker=m)

    ax.set_xlabel("Dim 1")
    ax.set_ylabel("Dim 2")
    ax.set_zlabel("Dim 3")

    plt.show()
    return None


def set_label_colors(ticklabels: list, colors: list) -> str:
    """
    Set the colours of labels

        Args:
            ticklabels (list): list of labels to change
            colors (list): list of colours to change to

        Returns:
            None
    """
    for label, color in zip(ticklabels, colors):
        plt.setp(
            label,
            backgroundcolor=color,
            bbox=dict(
                facecolor=color,
                alpha=0.5,
                # boxstyle="round, rounding_size=0.8",
                boxstyle="round, rounding_size=0.7",
                edgecolor="none",
            ),
        )
    return None


def cat_to_colour(categories: list, col_dict: dict) -> list[str]:
    """
    Convert a list of categories to a list of colours

        Args:
            categories (list): list of categories to convert
            col_dict (dict): dictionary with categories as keys and colours as values

        Returns:
            list: list of colours
    """
    colors = []
    for category in categories:
        if isinstance(category, list):
            if len(category) > 0:
                category = category[0]  # get colour for first one
            else:
                category = ""
        if category in col_dict:
            color = col_dict[category]
        else:
            color = "grey"
        colors.append(color)
    return colors


def set_label_colors_from_categories(
    ticklabels: list, categories: list, col_dict: list
) -> str:
    """
    Set the colours of labels from a list of categories

        Args:
            ticklabels (list): list of labels to change
            categories (list): list of categories to convert
            col_dict (dict): dictionary with categories as keys and colours as values

        Returns:
            None
    """
    colors = cat_to_colour(categories, col_dict)
    set_label_colors(ticklabels, colors)
    return None


def two_gradient_cmap():
    top = mpl.colormaps["Greys"].resampled(128)
    middle = mpl.colormaps["Greys"].resampled(128)
    bottom = mpl.colormaps["Greens"].resampled(128)

    newcolors = np.vstack(
        (top(np.linspace(0, 0.4, 128)), bottom(np.linspace(0.3, 1, 128)))
    )
    return mpl.ListedColormap(newcolors, name="GreenGrey")


def _triangulation_for_triheatmap(M, N):

    xv, yv = np.meshgrid(
        np.arange(-0.5, M), np.arange(-0.5, N)
    )  # vertices of the little squares
    xc, yc = np.meshgrid(
        np.arange(0, M), np.arange(0, N)
    )  # centers of the little squares
    x = np.concatenate([xv.ravel(), xc.ravel()])
    y = np.concatenate([yv.ravel(), yc.ravel()])
    cstart = (M + 1) * (N + 1)  # indices of the centers

    north = [
        (i + j * (M + 1), i + 1 + j * (M + 1), cstart + i + j * M)
        for j in range(N)
        for i in range(M)
    ]
    east = [
        (i + 1 + j * (M + 1), i + 1 + (j + 1) * (M + 1), cstart + i + j * M)
        for j in range(N)
        for i in range(M)
    ]
    south = [
        (i + 1 + (j + 1) * (M + 1), i + (j + 1) * (M + 1), cstart + i + j * M)
        for j in range(N)
        for i in range(M)
    ]
    west = [
        (i + (j + 1) * (M + 1), i + j * (M + 1), cstart + i + j * M)
        for j in range(N)
        for i in range(M)
    ]
    return [Triangulation(x, y, triangles) for triangles in [north, east, south, west]]


def triheatmap(
    data1: np.array,
    data2: np.array,
    row_labels: list = [],
    col_labels: list = [],
    ax=None,
    cbar_kw=None,
    cbarlabel1="% of natural products",
    cbarlabel2="% of compounds",
    title="Substructure counts for natural products vs synthetic compounds",
    colours=("Greens", "Purples"),
    normalise=False,
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

    values = [data1, data2, data2, data1]
    triangle_ = _triangulation_for_triheatmap(data1.shape[1], data1.shape[0])
    cmaps = [colours[0], colours[1], colours[1], colours[0]]
    # Plot the heatmap
    # im = ax.imshow(data, **kwargs)
    if normalise:
        norms = [plt.Normalize(0, 100) for _ in range(4)]
        imgs = [
            ax.tripcolor(t, np.ravel(val), cmap=cmap, norm=norm)
            for t, val, cmap, norm in zip(triangle_, values, cmaps, norms)
        ]
    else:
        imgs = [
            ax.tripcolor(t, np.ravel(val), cmap=cmap)
            for t, val, cmap in zip(triangle_, values, cmaps)
        ]

    # Create colorbar
    cbar = ax.figure.colorbar(imgs[1], ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel2, rotation=-90, va="bottom")
    cbar.ax.set_yticks([])
    cbar.outline.set_visible(False)

    # Show all ticks and label them with the respective list entries.
    ax.set_xticks(np.arange(data1.shape[1]), labels=col_labels, size=6)
    # ax.set_yticks(np.arange(data.shape[0]), labels=row_labels, size = 6)
    ax.set_yticks([])

    # Let the horizontal axes labeling appear on bottom (default)
    ax.tick_params(top=False, bottom=True, labeltop=False, labelbottom=True)
    ax.tick_params(pad=0.5)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=30, ha="right", rotation_mode="anchor")

    # Turn spines off and create white grid.
    ax.spines[:].set_visible(False)

    ax.invert_yaxis()
    ax.margins(x=0, y=0)
    ax.set_xticks(np.arange(data1.shape[1] + 1) - 0.5, minor=True)
    ax.set_yticks(np.arange(data1.shape[0] + 1) - 0.5, minor=True)
    ax.set_aspect("equal", "box")  # square cells

    ax.grid(which="minor", color="w", linestyle="-", linewidth=1, alpha=1.0)
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.set_xlabel("substructure")
    ax.set_ylabel("counts", labelpad=10)
    ax.set_title(title, loc="center", pad=20)
    plt.tight_layout()

    return imgs, cbar
