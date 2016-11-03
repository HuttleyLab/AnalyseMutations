#!/usr/bin/env python
from configparser import SafeConfigParser
import os
from matplotlib import pyplot
from matplotlib.ticker import FuncFormatter, MultipleLocator
from numpy import array
import click
from scitrack import CachingLogger
from mutation_motif import util, draw, logo

LOGGER = CachingLogger(create_dir=True)


def read_fig_config(path):
    LOGGER.input_file(path.name)
    parser = SafeConfigParser()
    parser.optionxform = str  # stops automatic conversion to lower case
    parser.readfp(path)

    nrows = int(parser.get('fig setup', 'num_rows'))
    ncols = int(parser.get('fig setup', 'num_cols'))
    figsize = eval(parser.get('fig setup', 'figsize'))
    col_labels = [l.strip()
                  for l in parser.get('fig setup', 'col_labels').split(',')]
    row_labels = [l.strip()
                  for l in parser.get('fig setup', 'row_labels').split(',')]

    # get individual plot config
    axis_cfg = {}
    for ax in ['x', 'y']:
        for attr in ['label', 'tick']:
            label = "%s%s_fontsize" % (ax, attr)
            axis_cfg[label] = int(parser.get('1-way plot', label))
        label = label = "%slabel_pad" % ax
        axis_cfg[label] = float(parser.get('1-way plot', label))

    # now get path for each section, converting section head into python
    # indices
    json_paths = {}
    for col in range(ncols):
        for row in range(nrows):
            sect = "%d,%d" % (col + 1, row + 1)

            path = parser.get(sect, "path")
            json_paths[(row, col)] = path  # because of how numpy arrays work
    return ncols, nrows, figsize, col_labels, row_labels, json_paths, axis_cfg


@click.command()
@click.argument("fig_config", type=click.File())
@click.argument("output_path", required=False, type=click.Path())
def main(fig_config, output_path):
    # we read in the config file and determine number of rows and columns
    # paths, headings, etc ..
    # then create the figure and axes and call the mutation_motif drawing code
    if not output_path:
        output_path = "delme.pdf"
    else:
        util.makedirs(os.path.dirname(output_path))

    ncols, nrows, figsize, col_labels, row_labels, paths, axis_cfg = \
        read_fig_config(fig_config)

    fig, axes = pyplot.subplots(nrows=nrows, ncols=ncols, figsize=figsize,
                                sharex=True, sharey=True)
    figwidth = fig.get_figwidth()
    try:
        axes[0]
    except TypeError:
        axes = array([[axes]])

    adaptive_y = 0
    plottable = {}
    for coord in paths:
        data = util.load_loglin_stats(paths[coord])
        positions = list(data)
        positions.sort()
        heights, characters, indices = draw.get_plot_data(data, positions)
        adaptive_y = max(adaptive_y, logo.est_ylim(heights))
        plottable[coord] = dict(char_heights=heights,
                                characters=characters,
                                position_indices=indices,
                                figwidth=figwidth,
                                verbose=False)

    ylim = adaptive_y
    for coord in plottable:
        kwargs = plottable[coord]
        kwargs["ax"] = axes[coord]
        kwargs["ylim"] = ylim
        fig = logo.draw_multi_position(**kwargs)

    xformat = FuncFormatter(draw.format_float(1e-3, float_places=2))

    for col in range(ncols):
        top_ax = axes[0, col]
        top_ax.set_title(col_labels[col], fontsize=axis_cfg["xlabel_fontsize"],
                         weight="bold", y=1.1)
        btm_ax = axes[-1, col]
        for xticklabel in btm_ax.get_xticklabels():
            xticklabel.set_fontsize(axis_cfg["xtick_fontsize"])
            xticklabel.set_rotation(0)
        btm_ax.set_xlabel("Position", fontsize=axis_cfg["xlabel_fontsize"],
                          weight="bold")
        btm_ax.xaxis.labelpad = axis_cfg['xlabel_pad']

    for row in range(nrows):
        lft_ax = axes[row, 0]
        for yticklabel in lft_ax.get_yticklabels():
            yticklabel.set_fontsize(axis_cfg["ytick_fontsize"])
            yticklabel.set_rotation(0)

        lft_ax.yaxis.set_major_formatter(FuncFormatter(xformat))
        lft_ax.yaxis.labelpad = axis_cfg['ylabel_pad']
        lft_ax.set_ylabel(row_labels[row], rotation=0,
                          fontsize=axis_cfg['ylabel_fontsize'],
                          weight="bold")

    fig.tight_layout()
    fig.savefig(output_path)
    click.secho("Wrote %s" % output_path)


if __name__ == "__main__":
    main()
