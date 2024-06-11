import argparse, logging
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


def cli():
    """
    Parse command line arguments
    """
    parser = argparse.ArgumentParser(
        description="Create a heatmap from a confusion matrix"
    )
    parser.add_argument(
        "matrix_path",
        type=str,
        help="path to confusion matrix",
    )
    args = parser.parse_args()
    return args


def set_style() -> None:
    """
    Set the style of the plot to biostylefoni
    """
    # get path of this script
    script_path = os.path.dirname(os.path.realpath(__file__))
    parent_path = os.path.dirname(script_path)
    utils_path = os.path.join(parent_path, "utils")
    print(utils_path)
    style_path = os.path.join(utils_path, "biostylefoni.mplstyle")
    # set style
    plt.style.use(style_path)
    return None


def parse_cms_files(m_path: str) -> tuple[np.array, np.array]:
    """
    Parse the confusion matrices and names from a file

    Args:
        m_path (str): path to file with confusion matrices

    Returns:
        tuple[np.array, np.array]: array of confusion matrices and array of names (i.e. )
    """
    cms = np.loadtxt(
        m_path,
        delimiter="\t",
        dtype=int,
        skiprows=1,
        usecols=(1, 2, 3, 4),
    )
    # for each row, split the array into a 2x2 array
    cms = cms.reshape(cms.shape[0], 2, 2)

    # now get only the 'index names'
    cm_names = np.loadtxt(
        m_path,
        delimiter="\t",
        dtype=str,
        skiprows=1,
        usecols=(0),
    )
    return cms, cm_names


def get_matrices(cms: np.array) -> tuple[list[np.array], list[np.array]]:
    """
    Take an array of confusion matrices and return a list of matrices and a list of normalised matrices

    Args:
        cms (np.array): array of confusion matrices

    Returns:
        tuple[list[np.array], list[np.array]]: a list of matrices and a list of normalised matrices
    """
    matrices = []
    norm_matrices = []
    perc_matrices = []
    for i in range(cms.shape[0]):
        # # make random matrix with values between 0 and 100000
        # matrix = np.random.randint(0, 100000, size=(2, 2))
        matrix = cms[i]
        # normalise matrix
        norm_matrix = matrix / matrix.sum(axis=1, keepdims=True)
        # turn normalised matrix into percentages
        perc_matrix = norm_matrix * 100
        # append matrices to list
        matrices.append(matrix)
        norm_matrices.append(norm_matrix)
        perc_matrices.append(perc_matrix)
    assert len(matrices) == len(perc_matrices), "#matrices don't match"
    return matrices, perc_matrices


def main():
    # get arguments
    args = cli()

    # set style
    set_style()

    # get absolute path from relative path
    m_path = os.path.abspath(args.matrix_path)

    # get folder name
    folder = m_path.split("/")[-3:-1]
    ml_input = "_".join(folder)

    # read in the confusion matrix and names
    cms, cm_names = parse_cms_files(m_path)
    # print(cms, cm_names)

    # get the matrices and the percentage versions for each category
    matrices, perc_matrices = get_matrices(cms)
    assert len(matrices) == len(cm_names), "#matrices and #categories don't match"

    # make subplots
    fig, axs = plt.subplots(
        1, len(matrices), figsize=(len(matrices), 2), dpi=500
    )  # , sharey=True) #sharing y makes the y axis ticks appear in each subplot

    # make a heatmap in each subplot
    for i, ax in enumerate(axs):
        cmap_name = "Greys"
        cmap = mpl.colormaps[cmap_name]
        # plot heatmap
        im = ax.imshow(perc_matrices[i], cmap=cmap_name, vmin=0, vmax=100)
        # remove ticks
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        # set title
        title = cm_names[i].replace("_", "\n")
        axtitle = ax.set_title(title, fontsize=7, fontweight=600, wrap=True)
        # force the wrap line width to be shorter
        axtitle._get_wrap_line_width = lambda: 600.0  #  wrap to 600 screen pixels
        # annotate values in each box, with dark text for light background and light text for dark background
        fontweight = 500
        a_size = 5
        for j in range(2):
            for k in range(2):
                if perc_matrices[i][j][k] < 50:
                    text = ax.text(
                        k,
                        j,
                        f"{round(matrices[i][j][k], 2)}\n({round(perc_matrices[i][j][k], 2)}%)",
                        ha="center",
                        va="center",
                        color=cmap(1.0),
                        fontsize=a_size,
                        fontweight=fontweight,
                    )
                else:
                    text = ax.text(
                        k,
                        j,
                        # round(matrices[i][j][k], 2),
                        f"{round(matrices[i][j][k], 2)}\n({round(perc_matrices[i][j][k], 2)}%)",
                        ha="center",
                        va="center",
                        color=cmap(0.0),
                        fontsize=a_size,
                        fontweight=fontweight,
                    )

    # set ticks
    axs[0].set_yticks([0, 1])
    axs[0].set_yticklabels(["P", "N"], fontsize=6, fontweight=500)
    axs[0].set_xticks([0, 1])
    axs[0].set_xticklabels(["P", "N"], fontsize=6, fontweight=500)

    # set y label
    axs[0].set_ylabel("truth", fontsize=7, fontweight=600)
    axs[0].set_xlabel("prediction", fontsize=7, fontweight=600)

    # set common title
    fig.suptitle(
        f"confusion matrices for multilabel RF on {ml_input}",
        fontsize=9,
        fontweight=600,
    )

    # set suptitle on y axis (for later when looping across all folders)
    # fig.text(0.02, 0.5, "confusion matrices for multilable RF on", fontsize=8, fontweight=600, rotation=90, va="center")

    # set tight layout
    plt.tight_layout()

    # set colorbar
    cbar = fig.colorbar(im, ax=axs, shrink=0.8, orientation="horizontal")
    cbar.ax.tick_params(labelsize=6)
    # add label to colorbar
    cbar.ax.set_xlabel("% of compounds", fontsize=7)

    # get path to folder where confusion matrices are, using os.path.dirname
    save_path = "/".join(m_path.split("/")[:-1])
    logging.info(m_path, save_path)
    # save figure
    plt.savefig(
        m_path.replace("_matrix.txt", "_heatmap.png"), dpi=500, bbox_inches="tight"
    )
    plt.clf()
    plt.close()
    exit(0)
    return None


if __name__ == "__main__":
    main()
