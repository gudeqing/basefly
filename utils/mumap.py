import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

"""
希望基于python模仿maf的瀑布图
未开发好
"""

def get_test_data():
    genes = ["A", "B", "C"]
    samples = ["Sample1", "Sample2", "Sample3", "Sample4"]

    harvest = np.array([
        ['SNV', 'Indel', 'Insertion','SNV'],
        ['SNV', 'Indel', 'Insertion','deletion'],
        ['SNV', 'Indel', 'Insertion','SNV'],
    ])
    data = pd.DataFrame(harvest, columns=samples, index=genes)
    data.index.name = 'gene'
    data.columns.name = 'sample'

    # 直方图
    count_df = data.apply(lambda x: x.value_counts())
    count_df.index.name = 'mutation_type'
    print(count_df)
    print(count_df.index)
    print(count_df.columns)
    return data, count_df


def heatmap(df:pd.DataFrame, ax=None, cbar_kw={}, show_cbar=False, cbarlabel="", **kwargs):
    if not ax:
        ax = plt.gca()
    row_labels = list(df.index)
    col_labels = list(df.columns)

    # Plot the heatmap
    im = ax.imshow(df.values, **kwargs)

    # Create colorbar
    if show_cbar:
        cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
        cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")
    else:
        cbar = None

    # Show all ticks and label them with the respective list entries.
    ax.xaxis.set_ticks(np.arange(df.shape[1]))
    ax.xaxis.set_ticklabels(col_labels)
    ax.yaxis.set_ticks(np.arange(df.shape[0]))
    ax.yaxis.set_ticklabels(row_labels)

    # Let the horizontal axes labeling appear on top or bottom.
    ax.tick_params(top=True, bottom=True, labeltop=False, labelbottom=True)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=30, ha="right", rotation_mode="anchor")

    # Turn spines off and create white grid.
    ax.spines[:].set_visible(False)

    ax.set_xticks(np.arange(df.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(df.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar


def annotate_heatmap(im, data=None, valfmt="{x:.2f}",
                     textcolors=("black", "white"),
                     threshold=None, **textkw):
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
        threshold = im.norm(data.max())/2.

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center", verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
            text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
            texts.append(text)

    return texts


def layout():
    fig = plt.figure(constrained_layout=False)
    ax0 = fig.add_axes(rect=[0.01, 0.01, 0.98, 0.7])
    ax1 = fig.add_axes(rect=[0.01, 0.72, 0.98, 0.27])
    return ax0, ax1


def pipeline():
    ax0, ax1 = layout()
    data, count = get_test_data()
    value_types = sorted(set(data.values.flatten()))
    value_mapper = dict(zip(value_types, range(len(value_types))))
    data = data.applymap(lambda x:value_mapper[x])
    data.index = [x + ' 100%' for x in data.index]
    print(data)
    norm = matplotlib.colors.BoundaryNorm(range(len(value_types)+1), len(value_types))
    im, cbar = heatmap(data, ax=ax0, show_cbar=True,
                       cmap=plt.get_cmap("PiYG", len(value_types)),
                       norm=norm, cbarlabel="mutation_type")

    if cbar:
        cbar.ax.yaxis.set_ticks(np.arange(len(value_types))+0.5)
        cbar.ax.yaxis.set_ticklabels(value_types)

    ax0_position = ax0.get_position()
    ax1_position = ax1.get_position()
    ax1_position.x0 = ax0_position.x0
    ax1_position.x1 = ax0_position.x1
    ax1.set_position(ax1_position)
    print(ax1.get_position())

    count.T.plot.bar(ax=ax1, stacked=True)
    ax1.xaxis.set_ticklabels([])
    ax1.set_xlabel('')
    ax1.legend(loc=(1.01, 0.1))
    plt.savefig('tmp.pdf', bbox_inches='tight')

if __name__ == '__main__':
    pipeline()
