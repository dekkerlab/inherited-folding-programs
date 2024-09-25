import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd


from cooltools.api.dotfinder import bp_to_bins, generate_tiles_diag_band, get_adjusted_expected_tile_some_nans, clust_2D_pixels
from cooltools.lib.numutils import LazyToeplitz
from functools import partial, reduce
import logging
import multiprocess as mp
import matplotlib as mpl
import matplotlib.lines as lines
from matplotlib.lines import Line2D
from matplotlib.patches import ConnectionPatch, Rectangle
from mpl_toolkits.axes_grid1 import Divider, Size
from mpl_toolkits.axes_grid1.inset_locator import BboxConnector

from scipy.ndimage import convolve
logging.basicConfig(level=logging.INFO)
import warnings
from cytoolz import merge
import cooltools
import bbi


from bioframe.core import checks, construction
from bioframe.core.specs import _get_default_colnames, _verify_columns
from bioframe.core.stringops import parse_region


def merge_nested_intervals(starts, ends, overlap_frac=0.5):
    """
    Merge overlapping intervals.

    Parameters
    ----------
    starts, ends : numpy.ndarray
        Interval coordinates. Warning: if provided as pandas.Series, indices
        will be ignored.

    Returns
    -------
    cluster_ids : numpy.ndarray
        The indices of interval clusters that each interval belongs to.
    cluster_starts : numpy.ndarray
    cluster_ends : numpy.ndarray
        The spans of the merged intervals.

    Notes
    -----
    From
    https://stackoverflow.com/questions/43600878/merging-overlapping-intervals/58976449#58976449
    """

    for vec in [starts, ends]:
        if isinstance(vec, pd.Series):
            warnings.warn(
                "One of the inputs is provided as pandas.Series and its index "
                "will be ignored.",
                SyntaxWarning,
                stacklevel=2,
            )

    starts = np.asarray(starts)
    ends = np.asarray(ends)

    order = np.lexsort([ends, starts])
    starts, ends = starts[order], ends[order]

    ends = np.maximum.accumulate(ends)
    cluster_borders = np.zeros(len(starts) + 1, dtype=bool)
    cluster_borders[0] = True
    cluster_borders[-1] = True

    # ... MODIFICATION ...
    # cluster borders are where start of next interval is beyond the end of previous one
    # or the end of next interval "sticks out" relative to the end of the previous one ...
    _overlaps = ( ends[:-1] - starts[1:] )
    _sizes = (ends - starts)
    # _stickout_ends = ( ends[1:] - ends[:-1] )
    # _stickout_starts = ( starts[1:] - starts[:-1] )
    # assert (_stickout_starts >= 0).all()
    # _interval_sizes = ( ends - starts )
    # break clusters when there is no overlap at all, or the next interval mostly sticks out ...
    # cluster_borders[1:-1] = (_overlaps < 0) | ( _stickout_sizes > stickout_frac*_overlaps)
    # ..... or turn it around like so - 
    # when ovelap is very small compare to both stickout situations ...
    # basically we want to keep merging if an overlap is considerably larger either of the sizes ! (not stickouts)
    # so write that down as a statement and negate that - that should be it !
    cluster_borders[1:-1] = (
        (_overlaps < 0) | \
        (
            ( _overlaps < overlap_frac*_sizes[1:] ) & \
            ( _overlaps < overlap_frac*_sizes[:-1] )
        )
    )
    # cluster_borders[1:-1] = (_stickout_sizes > stickout_frac)
    # cluster_borders[1:-1] = (starts[1:] >= ends[:-1]) | (ends[1:] > ends[:-1])
    # ... END MODIFICATION ...

    cluster_ids_sorted = np.cumsum(cluster_borders)[:-1] - 1
    cluster_ids = np.full(starts.shape[0], -1)
    cluster_ids[order] = cluster_ids_sorted

    cluster_starts = starts[:][cluster_borders[:-1]]
    cluster_ends = ends[:][cluster_borders[1:]]

    return cluster_ids, cluster_starts, cluster_ends


def merge_nested(df, overlap_frac=0.5, cols=None, on=None):
    """
    Merge overlapping intervals.

    This returns a new dataframe of genomic intervals, which have the genomic
    coordinates of the interval cluster groups from the input dataframe. Also
    :func:`cluster()`, which returns the assignment of intervals to clusters
    prior to merging.

    Parameters
    ----------
    df : pandas.DataFrame

    overlap_frac : float
        How big an overlap between intervals should be in order for them to be
        merged. As a fraction of interval sizes.

    cols : (str, str, str) or None
        The names of columns containing the chromosome, start and end of the
        genomic intervals. The default values are 'chrom', 'start', 'end'.

    on : None or list
        List of column names to perform clustering on independently, passed as
        an argument to df.groupby before clustering. Default is None.
        An example useage would be to pass ``on=['strand']``.

    Returns
    -------
    df_merged : pandas.DataFrame
        A pandas dataframe with coordinates of merged clusters.

    Notes
    -------
    Resets index.

    """

    # Allow users to specify the names of columns containing the interval coordinates.
    ck, sk, ek = _get_default_colnames() if cols is None else cols
    checks.is_bedframe(df, raise_errors=True, cols=[ck, sk, ek])

    df = df.copy()
    df.reset_index(inplace=True, drop=True)

    # Find overlapping intervals for groups specified by on=[] (default on=None)
    group_list = [ck]
    if on is not None:
        if not isinstance(on, list):
            raise ValueError("on=[] must be None or list")
        if ck in on:
            raise ValueError("on=[] should not contain chromosome colnames")
        _verify_columns(df, on)
        group_list += on
    df_groups = df.groupby(group_list, observed=True).groups

    clusters = []

    for group_keys, df_group_idxs in df_groups.items():
        if pd.isna(pd.Series(group_keys)).any():
            continue
        if df_group_idxs.empty:
            continue

        df_group = df.loc[df_group_idxs]
        (
            cluster_ids_group,
            cluster_starts_group,
            cluster_ends_group,
        ) = merge_nested_intervals(
            df_group[sk].values.astype(np.int64),
            df_group[ek].values.astype(np.int64),
            overlap_frac=overlap_frac,
        )
        interval_counts = np.bincount(cluster_ids_group)
        n_clusters = cluster_starts_group.shape[0]

        ## Storing chromosome names causes a 2x slowdown. :(
        if isinstance(group_keys, str):
            group_keys = (group_keys,)
        clusters_group = {}
        for col in group_list:
            clusters_group[col] = pd.Series(
                data=np.full(n_clusters, group_keys[group_list.index(col)]),
                dtype=df[col].dtype,
            )
        clusters_group[sk] = cluster_starts_group
        clusters_group[ek] = cluster_ends_group
        clusters_group["n_intervals"] = interval_counts
        clusters_group = pd.DataFrame(clusters_group)

        clusters.append(clusters_group)

    df_nans = pd.isnull(df[[sk, ek, *group_list]]).any(axis=1)
    df_has_nans = df_nans.sum()
    if df_has_nans:
        nan_intervals = pd.DataFrame(
            [pd.NA] * df_has_nans,
            columns=["n_intervals"],
            index=df.loc[df_nans].index,
        )
        clusters.append(
            pd.concat(
                [df.loc[df_nans], nan_intervals],
                axis=1,
            )
        )

    clusters = pd.concat(clusters).reset_index(drop=True)
    if df_has_nans:
        clusters = clusters.astype(
            {sk: pd.Int64Dtype(), ek: pd.Int64Dtype(), "n_intervals": pd.Int64Dtype()}
        )

    # reorder cluster columns to have chrom,start,end first
    clusters_names = list(clusters.keys())
    clusters = clusters[
        [ck, sk, ek] + [col for col in clusters_names if col not in [ck, sk, ek]]
    ]

    return clusters



def get_stack(bigwig_path, df, kind="start", flank=5000, nbins=200, chrom_col="chrom", start_col="start", end_col="end"):
    """
    extract stackups ...
    """
    if kind=="start":
        _start = df.eval(f"({start_col}) - {flank}")
        _end = df.eval(f"({start_col}) + {flank}")
    elif kind=="end":
        _start = df.eval(f"({end_col}) - {flank}")
        _end = df.eval(f"({end_col}) + {flank}")
    elif kind=="mid":
        _start = df.eval(f"({end_col}+{start_col})//2 - {flank}")
        _end = df.eval(f"({end_col}+{start_col})//2 + {flank}")
    else:
        raise("kind can only be start,end or mid")

    return bbi.stackup(
        bigwig_path,
        df[chrom_col],
        _start,
        _end,
        bins=nbins,
    )

# update this t5o the Anne-Laure's version later on ...
def show_stacks(stackups, lims=None, order_idx=None, flank=5000, nbins=200, len_per_unit_depth = 0.005, _aspect=0.15):
    """
    show 4 stackups in the dict ...
    """

    if order_idx is None:
        # first key
        _k, *_ = stackups.keys()
        # get stack depth
        _stack_depth = stackups[_k].shape[0]
        # order based on depth ...
        order_idx = np.arange(_stack_depth)

    # len_per_unit_depth = 0.005
    width_per_unit_stack = 3.33

    fig, axs = plt.subplots(
        nrows=2,
        ncols=len(stackups),
        figsize=[width_per_unit_stack*len(stackups), len_per_unit_depth*len(order_idx)],
        # sharey=True,
        sharex=True,
        height_ratios=[0.2,1]
    )
    for ax, (name, _stack) in zip(axs[0], stackups.items()):
        ax.plot(_stack[order_idx].mean(axis=0))
        vmin = lims[name]["vmin"]
        vmax = lims[name]["vmax"]
        ax.set_ylim(vmin,vmax)
        ax.set_title(name)
        ax.set_xticks([0,nbins/2, nbins])
        # ax.set_xticklabels([f"-{flank//1000}kb", "", f"{flank//1000}kb"])

    for ax, (name, _stack) in zip(axs[1], stackups.items()):
        ax.imshow(
            _stack[order_idx],
            norm=mpl.colors.LogNorm() if lims is None else mpl.colors.LogNorm(**lims[name]),
            # interpolation="none",
            aspect=_aspect,
            cmap="Blues"
        )
        ax.set_xticks([0,nbins/2, nbins])
        ax.set_xticklabels([f"-{flank//1000}kb", "", f"{flank//1000}kb"])

    #fig.suptitle(bw_fname)



def saddleplot(
    track,
    saddledata,
    n_bins,
    vrange=None,
    qrange=(0.0, 1.0),
    cmap="coolwarm",
    scale="log",
    vmin=0.5,
    vmax=2,
    color=None,
    title=None,
    xlabel=None,
    ylabel=None,
    clabel=None,
    fig=None,
    fig_kws=None,
    heatmap_kws=None,
    margin_kws=None,
    cbar_kws=None,
    subplot_spec=None,
):
    """
    Generate a saddle plot.
    Parameters
    ----------
    track : pd.DataFrame
        See cooltools.digitize() for details.
    saddledata : 2D array-like
        Saddle matrix produced by `make_saddle`. It will include 2 flanking
        rows/columns for outlier signal values, thus the shape should be
        `(n+2, n+2)`.
    cmap : str or matplotlib colormap
        Colormap to use for plotting the saddle heatmap
    scale : str
        Color scaling to use for plotting the saddle heatmap: log or linear
    vmin, vmax : float
        Value limits for coloring the saddle heatmap
    color : matplotlib color value
        Face color for margin bar plots
    fig : matplotlib Figure, optional
        Specified figure to plot on. A new figure is created if none is
        provided.
    fig_kws : dict, optional
        Passed on to `plt.Figure()`
    heatmap_kws : dict, optional
        Passed on to `ax.imshow()`
    margin_kws : dict, optional
        Passed on to `ax.bar()` and `ax.barh()`
    cbar_kws : dict, optional
        Passed on to `plt.colorbar()`
    subplot_spec : GridSpec object
        Specify a subregion of a figure to using a GridSpec.
    Returns
    -------
    Dictionary of axes objects.
    """

#     warnings.warn(
#         "Generating a saddleplot will be deprecated in future versions, "
#         + "please see https://github.com/open2c_examples for examples on how to plot saddles.",
#         DeprecationWarning,
#     )

    from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
    from matplotlib.colors import Normalize, LogNorm
    from matplotlib import ticker
    import matplotlib.pyplot as plt

    class MinOneMaxFormatter(ticker.LogFormatter):
        def set_locs(self, locs=None):
            self._sublabels = set([vmin % 10 * 10, vmax % 10, 1])

        def __call__(self, x, pos=None):
            if x not in [vmin, 1, vmax]:
                return ""
            else:
                return "{x:g}".format(x=x)

    track_value_col = track.columns[3]
    track_values = track[track_value_col].values

    digitized_track, binedges = cooltools.digitize(
        track, n_bins, vrange=vrange, qrange=qrange
    )
    x = digitized_track[digitized_track.columns[3]].values.astype(int).copy()
    x = x[(x > -1) & (x < len(binedges) + 1)]

    # Old version
    # hist = np.bincount(x, minlength=len(binedges) + 1)

    groupmean = track[track.columns[3]].groupby(digitized_track[digitized_track.columns[3]]).mean()

    if qrange is not None:
        lo, hi = qrange
        binedges = np.linspace(lo, hi, n_bins + 1)

    # Barplot of mean values and saddledata are flanked by outlier bins
    n = saddledata.shape[0]
    X, Y = np.meshgrid(binedges, binedges)
    C = saddledata
    if (n - n_bins) == 2:
        C = C[1:-1, 1:-1]
        groupmean = groupmean[1:-1]

    # Layout
    if subplot_spec is not None:
        GridSpec = partial(GridSpecFromSubplotSpec, subplot_spec=subplot_spec)
    grid = {}
    gs = GridSpec(
        nrows=3,
        ncols=3,
        width_ratios=[0.2, 1, 0.1],
        height_ratios=[0.2, 1, 0.1],
        wspace=0.05,
        hspace=0.05,
    )

    # Figure
    if fig is None:
        fig_kws_default = dict(figsize=(5, 5))
        fig_kws = merge(fig_kws_default, fig_kws if fig_kws is not None else {})
        fig = plt.figure(**fig_kws)

    # Heatmap
    if scale == "log":
        norm = LogNorm(vmin=vmin, vmax=vmax)
    elif scale == "linear":
        norm = Normalize(vmin=vmin, vmax=vmax)
    else:
        raise ValueError("Only linear and log color scaling is supported")

    grid["ax_heatmap"] = ax = plt.subplot(gs[4])
    heatmap_kws_default = dict(cmap="coolwarm", rasterized=True)
    heatmap_kws = merge(
        heatmap_kws_default, heatmap_kws if heatmap_kws is not None else {}
    )
    img = ax.pcolormesh(X, Y, C, norm=norm, **heatmap_kws)
    plt.gca().yaxis.set_visible(False)

    # Margins
    margin_kws_default = dict(edgecolor="k", facecolor=color, linewidth=1)
    margin_kws = merge(margin_kws_default, margin_kws if margin_kws is not None else {})
    # left margin hist
    grid["ax_margin_y"] = plt.subplot(gs[3], sharey=grid["ax_heatmap"])

    plt.barh(
        binedges, height=1/len(binedges), width=groupmean, align="edge", **margin_kws
    )

    plt.xlim(plt.xlim()[1], plt.xlim()[0])  # fliplr
    plt.ylim(hi, lo)
    plt.gca().spines["top"].set_visible(False)
    plt.gca().spines["bottom"].set_visible(False)
    plt.gca().spines["left"].set_visible(False)
    plt.gca().xaxis.set_visible(False)
    # top margin hist
    grid["ax_margin_x"] = plt.subplot(gs[1], sharex=grid["ax_heatmap"])

    plt.bar(
        binedges, width=1/len(binedges), height=groupmean, align="edge", **margin_kws
    )

    plt.xlim(lo, hi)
    # plt.ylim(plt.ylim())  # correct
    plt.gca().spines["top"].set_visible(False)
    plt.gca().spines["right"].set_visible(False)
    plt.gca().spines["left"].set_visible(False)
    plt.gca().xaxis.set_visible(False)
    plt.gca().yaxis.set_visible(False)

#     # Colorbar
    grid["ax_cbar"] = plt.subplot(gs[5])
    cbar_kws_default = dict(fraction=0.8, label=clabel or "")
    cbar_kws = merge(cbar_kws_default, cbar_kws if cbar_kws is not None else {})
    if scale == "linear" and vmin is not None and vmax is not None:
        grid["ax_cbar"] = cb = plt.colorbar(img, **cbar_kws)
        # cb.set_ticks(np.arange(vmin, vmax + 0.001, 0.5))
        # # do linspace between vmin and vmax of 5 segments and trunc to 1 decimal:
        decimal = 10
        nsegments = 5
        cd_ticks = np.trunc(np.linspace(vmin, vmax, nsegments) * decimal) / decimal
        cb.set_ticks(cd_ticks)
    else:
        print('cbar')

        cb = plt.colorbar(img, format=MinOneMaxFormatter(), cax=grid["ax_cbar"], **cbar_kws)
        cb.ax.yaxis.set_minor_formatter(MinOneMaxFormatter())

    # extra settings
    grid["ax_heatmap"].set_xlim(lo, hi)
    grid["ax_heatmap"].set_ylim(hi, lo)
    grid['ax_heatmap'].grid(False)
    if title is not None:
        grid["ax_margin_x"].set_title(title)
    if xlabel is not None:
        grid["ax_heatmap"].set_xlabel(xlabel)
    if ylabel is not None:
        grid["ax_margin_y"].set_ylabel(ylabel)

    return grid



# function to draw kernels:
def draw_kernel(kernel, axis=None, kernel_name="default", cmap='viridis'):
    if axis is None:
        f, axis = plt.subplots()
    # kernel:
    imk = axis.imshow(
                    kernel[::-1,::-1],  # flip it, as in convolution
                    alpha=0.85,
                    cmap=cmap,
                    interpolation='nearest')
    # draw a square around the target pixel:
    x0 = kernel.shape[0] // 2 - 0.5
    y0 = kernel.shape[1] // 2 - 0.5
    rect = patches.Rectangle((x0, y0), 1, 1, lw=1, ec='r', fc='r')
    axis.add_patch(rect)

    # clean axis:
    axis.set_xticks([])
    axis.set_yticks([])
    axis.set_xticklabels('',visible=False)
    axis.set_yticklabels('',visible=False)
    axis.set_title(f"{kernel_name} kernel", fontsize=16)
    # add a checkerboard to highlight pixels:
    checkerboard = np.add.outer(range(kernel.shape[0]),
                                range(kernel.shape[1])) % 2
    # show it:
    axis.imshow(checkerboard,
            cmap='gray',
            interpolation='nearest',
            alpha=0.3)

    return imk


# create a functions that would return a series of rectangles around called dots
# in a specific region, and exposing importnat plotting parameters
def rectangles_around_dots(dots_bins_df, binstart_ij, loc="upper", lw=1, ec="cyan", fc="none"):
    rectangle_kwargs = dict(lw=lw, ec=ec, fc=fc)
    s1, s2 = binstart_ij
    for b1, b2 in dots_bins_df[["bin1_id", "bin2_id"]].itertuples(index=False):
        width1 = 1
        width2 = 1
        b1 -= s1
        b2 -= s2
        if loc == "upper":
            yield patches.Rectangle((b2, b1), width2, width1, **rectangle_kwargs)
        elif loc == "lower":
            yield patches.Rectangle((b1, b2), width1, width2, **rectangle_kwargs)
        else:
            raise ValueError("loc has to be uppper or lower")


def get_adjusted_expected_tile_OE( origin_ij, observed, expected, bal_weights, kernels ):
    # extract origin_ij coordinate of this tile:
    io, jo = origin_ij
    # let's extract full matrices and ice_vector:
    O_raw = observed  # raw observed, no need to copy, no modifications.
    E_bal = np.copy(expected)
    # 'bal_weights': ndarray or a couple of those ...
    if isinstance(bal_weights, np.ndarray):
        v_bal_i = bal_weights
        v_bal_j = bal_weights
    elif isinstance(bal_weights, (tuple, list)):
        v_bal_i, v_bal_j = bal_weights
    else:
        raise ValueError(
            "'bal_weights' must be an numpy.ndarray"
            "for slices of a matrix with diagonal-origin or"
            "a tuple/list of a couple of numpy.ndarray-s"
            "for a slice of matrix with an arbitrary origin."
        )
    # prepare matrix of balancing weights Ci*Cj
    bal_weights_ij = np.outer(v_bal_i, v_bal_j)
    # balanced observed, from raw-observed
    # by element-wise multiply:
    O_bal = np.multiply(O_raw, bal_weights_ij)
    # O_bal is separate from O_raw memory-wise.

    # OE:
    with np.errstate(divide="ignore", invalid="ignore"):
        OE_bal = np.divide(O_bal, E_bal)

    # fill lower triangle of O_bal and E_bal with NaNs
    # in order to prevent peak calling from the lower triangle
    OE_bal[np.tril_indices_from(OE_bal, k=(io - jo) - 1)] = np.nan

    # raw E_bal: element-wise division of E_bal[i,j] and
    # v_bal[i]*v_bal[j]:
    E_raw = np.divide(E_bal, bal_weights_ij)

    # a matrix of NaNs/Inf-s
    N_bal = ~np.isfinite(OE_bal)
    # fill them in with zeroes, preventing
    # NaNs during convolution:
    OE_bal[N_bal] = 0.0
    # we are going to accumulate all the results
    # into a DataFrame, keeping NaNs, and other
    # unfiltered results (even the lower triangle for now):
    i, j = np.indices(O_raw.shape)
    # pack it into DataFrame to accumulate results:
    peaks_df = pd.DataFrame(
        {
            "bin1_id": i.ravel() + io,
            "bin2_id": j.ravel() + jo,
            "count": O_raw.ravel(),
            "expected_raw": E_raw.ravel(),  # new column for better exploration
            "oe": OE_bal.ravel(),  # regular bserved/expected per pixel
        }
    )

    # common kwargs for convolution:
    conv_kwargs = dict(mode="constant", origin=0)

    with np.errstate(divide="ignore", invalid="ignore"):
        for kernel_name, kernel in kernels.items():
            # a matrix filled with the kernel-weighted sums
            # based on a balanced observed/expected matrix:
            KOE = convolve(OE_bal, kernel, cval=0.0, **conv_kwargs)
            # get number of NaNs near each pixel (kernel's nonzero footprint)
            kernel_footprint = (kernel != 0).astype(int)
            NN = convolve( N_bal.astype(int), kernel_footprint, cval=1, **conv_kwargs)
            # get number of zeros near each pixel (zeros+nans that would be)
            NZN = convolve(np.isclose(OE_bal, 0.0).astype(int), kernel_footprint, cval=1, **conv_kwargs)
            # now finally, E_raw*(KOE), as the
            # locally-adjusted expected with raw counts as values:
            local_adjustment_factor = KOE / (kernel.sum() - NN)  # average OE in the kernel footprint
            Ek_raw = np.multiply(E_raw, local_adjustment_factor)
            #
            logging.debug(
                f"Convolution with kernel {kernel_name} is done for tile @ {io} {jo}."
            )
            # accumulation into single DataFrame:
            # store locally adjusted expected for each kernel
            # and number of NaNs in the footprint of each kernel
            peaks_df[f"la_exp.{kernel_name}.value"] = Ek_raw.ravel()
            peaks_df[f"la_exp.{kernel_name}.nnans"] = NN.ravel()
            peaks_df[f"la_exp.{kernel_name}.zeros"] = (NZN-NN).ravel()
            # NEW STUFF FOR DEEPER DIVE: store convolution ratio itself
            peaks_df[f"convolution_ratio.{kernel_name}"] = local_adjustment_factor.ravel()
            # division by KE=0 has to be treated separately:
            peaks_df[f"safe_division.{kernel_name}"] = np.isfinite(local_adjustment_factor.ravel())
            # do all the filter/logic/masking etc on the complete DataFrame ...

    return peaks_df


# derived from cooltools.api.dotfinder.score_tile to return more columns in the output:
def score_tile_custom_cols_OE(
    tile_cij,
    clr,
    expected_indexed,
    expected_value_col,
    clr_weight_name,
    kernels,
    max_nans_tolerated,
    band_to_cover,
    cols_to_return = None,
):
    # unpack tile's coordinates
    region_name, tile_span_i, tile_span_j = tile_cij
    tile_start_ij = (tile_span_i[0], tile_span_j[0])
    # we have to do it for every tile, because
    # region_name is not known apriori (maybe move outside)
    # use .loc[region, region] for symmetric cis regions to conform with expected v1.0
    lazy_exp = LazyToeplitz(
        expected_indexed.loc[region_name, region_name][expected_value_col].to_numpy()
    )
    # RAW observed matrix slice:
    observed = clr.matrix(balance=False)[slice(*tile_span_i), slice(*tile_span_j)]
    # expected as a rectangular tile :
    expected = lazy_exp[slice(*tile_span_i), slice(*tile_span_j)]
    # slice of balance_weight for row-span and column-span :
    bal_weight_i = clr.bins()[slice(*tile_span_i)][clr_weight_name].to_numpy()
    bal_weight_j = clr.bins()[slice(*tile_span_j)][clr_weight_name].to_numpy()
    # do the convolutions
    result = get_adjusted_expected_tile_OE(
        origin_ij=tile_start_ij,
        observed=observed,
        expected=expected,
        bal_weights=(bal_weight_i, bal_weight_j),
        kernels=kernels,
    )
    # Post-processing filters
    # (0) keep only upper-triangle pixels:
    upper_band = result["bin1_id"] < result["bin2_id"]
    # (1) exclude pixels that connect loci further than 'band_to_cover' apart:
    is_inside_band = result["bin1_id"] > (result["bin2_id"] - band_to_cover)
    # (2) identify pixels that pass number of NaNs compliance test for ALL kernels:
    does_comply_nans = np.all(
        result[[f"la_exp.{k}.nnans" for k in kernels]] < max_nans_tolerated, axis=1
    )
    # (3) keep pixels without nan/infinite local adjustment factors for all kernel
    finite_values_only = result[[f"safe_division.{k}" for k in kernels]].all(axis="columns")
    #
    # so, selecting inside band and nNaNs compliant results:
    res_df = result[upper_band & is_inside_band & does_comply_nans & finite_values_only].reset_index(
        drop=True
    )
    # return larger subset of columns ...
    if cols_to_return is None:
        cols_to_return = ["bin1_id", "bin2_id", "count"]
        cols_to_return += [f"la_exp.{k}.value" for k in kernels]

    return res_df[ cols_to_return ].astype(dtype={f"la_exp.{k}.value": "float64" for k in kernels})

def score_pixels_only_OE(
    clr,
    expected_indexed,
    expected_value_col,
    clr_weight_name,
    tiles,
    kernels,
    max_nans_tolerated,
    loci_separation_bins,
    nproc,
    cols_to_return = None,
):
    logging.info(f"convolving {len(tiles)} tiles to build histograms for lambda-bins")
    # to score per tile - a function of a single argument :
    _job = partial(
        score_tile_custom_cols_OE,
        clr=clr,
        expected_indexed=expected_indexed,
        expected_value_col=expected_value_col,
        clr_weight_name=clr_weight_name,
        kernels=kernels,
        max_nans_tolerated=max_nans_tolerated,
        band_to_cover=loci_separation_bins,
        cols_to_return=cols_to_return,
    )
    # standard multiprocessing implementation
    if nproc > 1:
        logging.info(f"creating a Pool of {nproc} workers to tackle {len(tiles)} tiles")
        pool = mp.Pool(nproc)
        map_ = pool.imap
        map_kwargs = dict(chunksize=int(np.ceil(len(tiles) / nproc)))
    else:
        logging.info("fallback to serial implementation.")
        map_ = map
        map_kwargs = {}
    try:
        # consider using
        # https://github.com/mirnylab/cooler/blob/9e72ee202b0ac6f9d93fd2444d6f94c524962769/cooler/tools.py#L59
        scored_df_chunks = map_(_job, tiles, **map_kwargs)
    finally:
        if nproc > 1:
            pool.close()

    return scored_df_chunks


def get_adjusted_expected_tile_more_columns( origin_ij, observed, expected, bal_weights, kernels ):
    # extract origin_ij coordinate of this tile:
    io, jo = origin_ij
    # let's extract full matrices and ice_vector:
    O_raw = observed  # raw observed, no need to copy, no modifications.
    E_bal = np.copy(expected)
    # 'bal_weights': ndarray or a couple of those ...
    if isinstance(bal_weights, np.ndarray):
        v_bal_i = bal_weights
        v_bal_j = bal_weights
    elif isinstance(bal_weights, (tuple, list)):
        v_bal_i, v_bal_j = bal_weights
    else:
        raise ValueError(
            "'bal_weights' must be an numpy.ndarray"
            "for slices of a matrix with diagonal-origin or"
            "a tuple/list of a couple of numpy.ndarray-s"
            "for a slice of matrix with an arbitrary origin."
        )
    # prepare matrix of balancing weights Ci*Cj
    bal_weights_ij = np.outer(v_bal_i, v_bal_j)
    # balanced observed, from raw-observed
    # by element-wise multiply:
    O_bal = np.multiply(O_raw, bal_weights_ij)
    # O_bal is separate from O_raw memory-wise.

    # fill lower triangle of O_bal and E_bal with NaNs
    # in order to prevent peak calling from the lower triangle
    # and also to provide fair locally adjusted expected
    # estimation for pixels very close to diagonal, whose
    # "donuts"(kernels) would be crossing the main diagonal.
    # The trickiest thing here would be dealing with the origin_ij: io,jo.
    O_bal[np.tril_indices_from(O_bal, k=(io - jo) - 1)] = np.nan
    E_bal[np.tril_indices_from(E_bal, k=(io - jo) - 1)] = np.nan

    # raw E_bal: element-wise division of E_bal[i,j] and
    # v_bal[i]*v_bal[j]:
    E_raw = np.divide(E_bal, bal_weights_ij)

    # let's calculate a matrix of common NaNs
    # shared between observed and expected:
    # check if it's redundant ? (is NaNs from O_bal sufficient? )
    N_bal = np.logical_or(np.isnan(O_bal), np.isnan(E_bal))
    # fill in common nan-s with zeroes, preventing
    # NaNs during convolution:
    O_bal[N_bal] = 0.0
    E_bal[N_bal] = 0.0
    E_raw[N_bal] = 0.0
    # think about usinf copyto and where functions later:
    # https://stackoverflow.com/questions/6431973/how-to-copy-data-from-a-numpy-array-to-another
    # #
    # we are going to accumulate all the results
    # into a DataFrame, keeping NaNs, and other
    # unfiltered results (even the lower triangle for now):
    i, j = np.indices(O_raw.shape)
    # pack it into DataFrame to accumulate results:
    with np.errstate(divide="ignore", invalid="ignore"):
        peaks_df = pd.DataFrame(
            {
                "bin1_id": i.ravel() + io,
                "bin2_id": j.ravel() + jo,
                "count": O_raw.ravel(),
                "expected_raw": E_raw.ravel(),  # new column for better exploration
                "oe": (O_bal/E_bal).ravel(),  # regular bserved/expected per pixel
            }
        )

    # common kwargs for convolution:
    conv_kwargs = dict(mode="constant", origin=0)

    with np.errstate(divide="ignore", invalid="ignore"):
        for kernel_name, kernel in kernels.items():
            # a matrix filled with the kernel-weighted sums
            # based on a balanced observed matrix:
            KO = convolve(O_bal, kernel, cval=0.0, **conv_kwargs)
            # a matrix filled with the kernel-weighted sums
            # based on a balanced expected matrix:
            KE = convolve(E_bal, kernel, cval=0.0, **conv_kwargs)
            # get number of NaNs in a vicinity of every
            # pixel (kernel's nonzero footprint)
            # based on the NaN-matrix N_bal.
            # N_bal is shared NaNs between O_bal E_bal,
            kernel_footprint = (kernel != 0).astype(int)
            NN = convolve( N_bal.astype(int), kernel_footprint, cval=1, **conv_kwargs)
            # now finally, E_raw*(KO/KE), as the
            # locally-adjusted expected with raw counts as values:
            local_adjustment_factor = np.divide(KO, KE)
            Ek_raw = np.multiply(E_raw, local_adjustment_factor)
            #
            logging.debug(
                f"Convolution with kernel {kernel_name} is done for tile @ {io} {jo}."
            )
            # accumulation into single DataFrame:
            # store locally adjusted expected for each kernel
            # and number of NaNs in the footprint of each kernel
            peaks_df[f"la_exp.{kernel_name}.value"] = Ek_raw.ravel()
            peaks_df[f"la_exp.{kernel_name}.nnans"] = NN.ravel()
            # NEW STUFF FOR DEEPER DIVE: store convolution ratio itself
            peaks_df[f"convolution_ratio.{kernel_name}"] = local_adjustment_factor.ravel()
            # division by KE=0 has to be treated separately:
            peaks_df[f"safe_division.{kernel_name}"] = np.isfinite(local_adjustment_factor.ravel())
            # do all the filter/logic/masking etc on the complete DataFrame ...

    return peaks_df


# derived from cooltools.api.dotfinder.score_tile to return more columns in the output:
def score_tile_custom_cols(
    tile_cij,
    clr,
    expected_indexed,
    expected_value_col,
    clr_weight_name,
    kernels,
    max_nans_tolerated,
    band_to_cover,
    cols_to_return = None,
):
    # unpack tile's coordinates
    region_name, tile_span_i, tile_span_j = tile_cij
    tile_start_ij = (tile_span_i[0], tile_span_j[0])
    # we have to do it for every tile, because
    # region_name is not known apriori (maybe move outside)
    # use .loc[region, region] for symmetric cis regions to conform with expected v1.0
    lazy_exp = LazyToeplitz(
        expected_indexed.loc[region_name, region_name][expected_value_col].to_numpy()
    )
    # RAW observed matrix slice:
    observed = clr.matrix(balance=False)[slice(*tile_span_i), slice(*tile_span_j)]
    # expected as a rectangular tile :
    expected = lazy_exp[slice(*tile_span_i), slice(*tile_span_j)]
    # slice of balance_weight for row-span and column-span :
    bal_weight_i = clr.bins()[slice(*tile_span_i)][clr_weight_name].to_numpy()
    bal_weight_j = clr.bins()[slice(*tile_span_j)][clr_weight_name].to_numpy()
    # do the convolutions
    result = get_adjusted_expected_tile_more_columns(
        origin_ij=tile_start_ij,
        observed=observed,
        expected=expected,
        bal_weights=(bal_weight_i, bal_weight_j),
        kernels=kernels,
    )
    # Post-processing filters
    # (0) keep only upper-triangle pixels:
    upper_band = result["bin1_id"] < result["bin2_id"]
    # (1) exclude pixels that connect loci further than 'band_to_cover' apart:
    is_inside_band = result["bin1_id"] > (result["bin2_id"] - band_to_cover)
    # (2) identify pixels that pass number of NaNs compliance test for ALL kernels:
    does_comply_nans = np.all(
        result[[f"la_exp.{k}.nnans" for k in kernels]] < max_nans_tolerated, axis=1
    )
    # (3) keep pixels without nan/infinite local adjustment factors for all kernel
    finite_values_only = result[[f"safe_division.{k}" for k in kernels]].all(axis="columns")
    #
    # so, selecting inside band and nNaNs compliant results:
    res_df = result[upper_band & is_inside_band & does_comply_nans & finite_values_only].reset_index(
        drop=True
    )
    # return larger subset of columns ...
    if cols_to_return is None:
        cols_to_return = ["bin1_id", "bin2_id", "count"]
        cols_to_return += [f"la_exp.{k}.value" for k in kernels]

    return res_df[ cols_to_return ].astype(dtype={f"la_exp.{k}.value": "float64" for k in kernels})


def score_pixels_only(
    clr,
    expected_indexed,
    expected_value_col,
    clr_weight_name,
    tiles,
    kernels,
    max_nans_tolerated,
    loci_separation_bins,
    nproc,
    cols_to_return = None,
):
    logging.info(f"convolving {len(tiles)} tiles to build histograms for lambda-bins")
    # to score per tile - a function of a single argument :
    _job = partial(
        score_tile_custom_cols,
        clr=clr,
        expected_indexed=expected_indexed,
        expected_value_col=expected_value_col,
        clr_weight_name=clr_weight_name,
        kernels=kernels,
        max_nans_tolerated=max_nans_tolerated,
        band_to_cover=loci_separation_bins,
        cols_to_return=cols_to_return,
    )
    # standard multiprocessing implementation
    if nproc > 1:
        logging.info(f"creating a Pool of {nproc} workers to tackle {len(tiles)} tiles")
        pool = mp.Pool(nproc)
        map_ = pool.imap
        map_kwargs = dict(chunksize=int(np.ceil(len(tiles) / nproc)))
    else:
        logging.info("fallback to serial implementation.")
        map_ = map
        map_kwargs = {}
    try:
        # consider using
        # https://github.com/mirnylab/cooler/blob/9e72ee202b0ac6f9d93fd2444d6f94c524962769/cooler/tools.py#L59
        scored_df_chunks = map_(_job, tiles, **map_kwargs)
    finally:
        if nproc > 1:
            pool.close()

    return scored_df_chunks



# modified clustering step that return all pixels associated with the cluster ...
def clustering_step_all(
    scored_df,
    dots_clustering_radius,
    assigned_regions_name="region",
):
    """
    Group together adjacent significant pixels into clusters after
    the lambda-binning multiple hypothesis testing by iterating over
    assigned regions and calling `clust_2D_pixels`.

    Parameters
    ----------
    scored_df : pandas.DataFrame
        DataFrame with enriched pixels that are ready to be
        clustered and are annotated with their genomic  coordinates.
    dots_clustering_radius : int
        Birch-clustering threshold.
    assigned_regions_name : str | None
        Name of the column in scored_df to use for grouping pixels
        before clustering. When None, full chromosome clustering is done.
    Returns
    -------
    centroids : pandas.DataFrame
        Pixels from 'scored_df' annotated with clustering information.

    Notes
    -----
    'dots_clustering_radius' in Birch clustering algorithm corresponds to a
    double the clustering radius in the "greedy"-clustering used in HiCCUPS

    """
    # make sure provided pixels are annotated with genomic corrdinates and raw counts column is present:
    if not {"chrom1", "chrom2", "start1", "start2"}.issubset(scored_df):
        raise ValueError("Scored pixels provided for clustering are not annotated")

    scored_df = scored_df.copy()
    if (
        not assigned_regions_name in scored_df.columns
    ):  # If input scores are not annotated by regions:
        logging.warning(
            f"No regions assigned to the scored pixels before clustering, using chromosomes"
        )
        scored_df[assigned_regions_name] = np.where(
            scored_df["chrom1"] == scored_df["chrom2"], scored_df["chrom1"], np.nan
        )

    # cluster within each regions separately and accumulate the result:
    pixel_clust_list = []
    scored_pixels_by_region = scored_df.groupby(assigned_regions_name, observed=True)
    for region, _df in scored_pixels_by_region:
        logging.info(f"clustering enriched pixels in region: {region}")
        # Using genomic corrdinated for clustering, not bin_id
        pixel_clust = clust_2D_pixels(
            _df,
            threshold_cluster=dots_clustering_radius,
            bin1_id_name="start1",
            bin2_id_name="start2",
        )
        pixel_clust_list.append(pixel_clust)
    logging.info("Clustering is complete")

    # concatenate clustering results ...
    # indexing information persists here ...
    if not pixel_clust_list:
        logging.warning("No clusters found for any regions! Output will be empty")
        empty_output = pd.DataFrame(
            [],
            columns=list(scored_df.columns)
            + [
                assigned_regions_name + "1",
                assigned_regions_name + "2",
                "c_label",
                "c_size",
                "cstart1",
                "cstart2",
            ],
        )
        return empty_output  # Empty dataframe with the same columns as anticipated
    else:
        pixel_clust_df = pd.concat(
            pixel_clust_list, ignore_index=False
        )  # Concatenate the clustering results for different regions

    # now merge pixel_clust_df and scored_df DataFrame ...
    # TODO make a more robust merge here
    df = pd.merge(
        scored_df, pixel_clust_df, how="left", left_index=True, right_index=True
    )
    # TODO check if next str-cast is neccessary
    df[assigned_regions_name + "1"] = df[assigned_regions_name].astype(str)
    df[assigned_regions_name + "2"] = df[assigned_regions_name].astype(str)
    # report ALL pixels per cluster instead of centroids ...
    # # report only centroids with highest Observed:
    # chrom_clust_group = df.groupby(
    #     [assigned_regions_name + "1", assigned_regions_name + "2", "c_label"],
    #     observed=True,
    # )
    # centroids = df.loc[
    #     chrom_clust_group[obs_raw_name].idxmax()
    # ]  # Select the brightest pixel in the cluster
    return df


import tempfile
import subprocess

def to_bigbed3(
    df,
    outpath,
    chromsizes,
    chrom_col = "chrom",
    start_col = "start",
    end_col = "end",
    cmd = "bedToBigBed",

):
    """
    bioframe's to_bigbed can only hanlde bed6
    actual ucsc-bedToBigBed can do bed3 easily !
    hence our little reimplemntation...
    """

    col_names = {
        chrom_col: "chrom",
        start_col: "start",
        end_col: "end",
    }
    columns = list(col_names.values())

    bed = df[[chrom_col,start_col,end_col]].rename(columns=col_names).copy()
    bed["chrom"] = bed["chrom"].astype(str)  # just in case ...
    bed = bed.sort_values(["chrom", "start", "end"])

    with tempfile.NamedTemporaryFile(suffix=".bed") as f, tempfile.NamedTemporaryFile(
        "wt", suffix=".chrom.sizes"
    ) as cs:

        chromsizes.to_csv(
            cs,
            sep="\t",
            header=False,
        )
        cs.flush()

        bed.to_csv(
            f.name,
            sep="\t",
            columns=columns,
            index=False,
            header=False,
            na_rep="nan",
        )

        p = subprocess.run(
            [cmd, f"-type=bed3", f.name, cs.name, outpath],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    return




#########################
#
#
#  stackup plotting functions ...
#
#
##########################

def _same_val(iterable_dict, func, return_val=True):
    """
    helper func checking if all values in 
    [func(elem) for elem in iterable_dict] are
    identical ...
    """
    _iter = iter(iterable_dict)
    _k = next(_iter)
    the_value = func(_k)
    # now check all remaining values ...
    if return_val:
        return the_value, all( func(_k) == the_value for _k in _iter )
    else:
        return all( func(_k) == the_value for _k in _iter )


def _get_hms_shape(hms_dict):
    """
    make sure heatmap-dicts are good for plotting ...
    i.e. they have the same length, width ? etc

    return number of stacks and the shape ...
    """
    the_shape, _shape_is_same = _same_val(hms_dict, lambda k: hms_dict[k].shape)
    if not _shape_is_same:
         raise ValueError('not all stacks have same shape !')
    #
    hm_height, hm_width = the_shape
    return len(hms_dict), hm_height, hm_width


def _get_hms_nested_shape(hms_dict_dict):
    """
    inspect a nested dictionary
    make sure heatmap-dicts are good for plotting ...
    i.e. they have the same length inside a group,
    they have the same width (*all*)
    number of stack per group is the same!
    ...

    return a bunch of info on the nested dict of stacks:
     - number of groups
     - stack samples per group (should be the same across groups)
     - heights of stacks per group (should be the same within each group)
     - stack width - should be the same across all groups/stacks
    """
    #
    num_groups = len(hms_dict_dict)
    group_hm_heights = {}
    group_hm_widths = {}
    group_amount = {}
    group_sample_keys = {}
    #
    for group_k, hms_dict in hms_dict_dict.items():
        # make sure shape of stackup in a group is identical:
        try:
            num_stacks, height, width = _get_hms_shape(hms_dict)
        except ValueError:
            raise ValueError(f'not all stacks in {group_k} have the same shape !')
        # accumulate heights, widths and amount per group ...
        group_hm_heights[group_k] = height
        group_hm_widths[group_k] = width
        group_amount[group_k] = num_stacks
        group_sample_keys[group_k] = list(hms_dict)
    # make sure widths are the same across groups:
    the_width, _width_is_same = _same_val(group_hm_widths, lambda k: group_hm_widths[k])
    if not _width_is_same:
         raise ValueError(f'not all groups of stacks have the same width !')
    # make sure amount of stacks is the same across groups:
    num_stacks_per_group, _amount_is_same = _same_val(group_amount, lambda k: group_amount[k])
    if not _amount_is_same:
         raise ValueError(f'not all groups have the same amount of stacks in them !')
    # make sure stacks (their sample keys) are the same across groups:
    stack_samples, _samples_are_same = _same_val(group_sample_keys, lambda k: group_sample_keys[k])
    if not _samples_are_same:
         raise ValueError(f'not all groups have the same stack samples !')
    # returning
    return num_groups, stack_samples, group_hm_heights, the_width


def _fillmissing_hms(hms_dict, how="col mean"):
    """
    modify heatmaps in the hms_dict by filling
    out the missing values ...

    there should be option on the fill values
    e.g. average from rows, average from columns
    a fixed value
    smth else ?
    how:
        "row mean"
        "col mean"
        float
    """
    hms_dict_filled = {}
    for k, hm in hms_dict.items():
            missing = ~np.isfinite(hm)
            if how == "col mean":
                mu = np.nanmean(hm, axis=0, keepdims=True)
            elif how == "row mean":
                mu = np.nanmean(hm, axis=1, keepdims=True)
            else:
                mu = float(how)  # should be inside try
            hms_dict_filled[k] = np.where(missing, mu, hm)
    return hms_dict_filled


def _get_profiles(hms_dict, scales):
    """
    for each stackup in the hms_dict generate an
    average profile ...
    """
    profile_hm = {}
    for k, hm in hms_dict.items():
        _scale = scales[k]
        if _scale == "log":
            profile_hm[k] = np.nanmean(hm, axis=0)
        else:
            profile_hm[k] = np.nanmean(hm, axis=0)
    # returning
    return profile_hm


def _get_norms(scales, vlims):
    """
    given scales and vlims - return norms !
    """
    norms = {}
    for k, scale in scales.items():
        # try to extract vmin , vmax ...
        try:
            vmin, vmax = vlims[k]
        except Exception:
            vmin, vmax = None, None
        # depending on scale ...
        if scale == "log":
            norms[k] = mpl.colors.LogNorm(vmin, vmax)
        else:
            norms[k] = mpl.colors.Normalize(vmin, vmax)
    # returning ...
    return norms


def _get_norm(scale, vlims):
    """
    given a scale and vlims - return norm !
    """
    # try to extract vmin , vmax ...
    try:
        vmin, vmax = vlims
    except Exception:
        vmin, vmax = None, None
    # depending on scale ...
    if scale == "log":
        return mpl.colors.LogNorm(vmin, vmax)
    else:
        return mpl.colors.Normalize(vmin, vmax)


def plot_stackups_lite(
    extra_plots,
    hms_dict,  # heatmaps dict, that controls the order, the amount etc etc ...
    scales,
    vlims,
    titles,
    cmaps,
    binsizes,
    fillmissing=False,
    len_per_thousand=1.2,
    width_per_stack=3.5,
    extra_height=3.0,  # height that goes toward the profile and colorbar
    # **plot_kwargs,
    fig_fontsize=25,
    **imshow_kwargs,
):
    """
    plot a buch of stackups ...
    """
    if extra_plots is None:
        extra_plots = []
    # how many stackups are there and what are the sizes ...
    num_stackups, stackup_items, hm_width = _get_hms_shape(hms_dict)
    num_cols = num_stackups + len(extra_plots)
    # figure out - how tall is this stackup
    stackup_height = stackup_items*len_per_thousand/1_000
    fig_height = stackup_height + extra_height
    fig_width = width_per_stack*num_cols
    fig = plt.figure(
        figsize=(fig_width, fig_height + 2*.8 + 2*.8),
        facecolor="white",
        tight_layout=False,
        constrained_layout=False,
    )
    # width_ratios = equal ...
    gs = fig.add_gridspec(
        nrows = 3,  # average plot, stackup and colorbar ...
        ncols = num_cols,
        height_ratios = [
            0.95*extra_height/fig_height,  # average profile
            (fig_height-extra_height)/fig_height,  # stackup itself
            0.05*extra_height/fig_height,  # colorbar ...
        ],
        left=0,
        right=1,
        top=1 - .8/(fig_height + 2*.8 + 2*.8),
        bottom=.8/(fig_height + 2*.8 + 2*.8),
        hspace=(3/fig_height)*.8,
        wspace=0.07,
    )

    ax_profile = {}
    ax_stackup = {}
    ax_xtra = {}  # default dict in pythoh - asre ordered already ...
    ax_cbar = {}

    # associate hms_dict elements with the gridspec items ...
    for idx, k in enumerate(hms_dict):
        ax_profile[k] = fig.add_subplot(gs[0,idx])
        ax_stackup[k] = fig.add_subplot(gs[1,idx])
        ax_cbar[k] = fig.add_subplot(gs[2,idx])

    hms_dict = _fillmissing_hms(hms_dict, how="col mean") if fillmissing else hms_dict
    profile_hm = _get_profiles(hms_dict, scales)
    norms = _get_norms(scales, vlims)

    # start plotting ...
    for k, hm in hms_dict.items():
        ax_profile[k].plot(profile_hm[k], linewidth=4)
        stack_hm = ax_stackup[k].imshow(
                          hm,
                          norm=norms[k],
                          aspect="auto",
                          cmap=cmaps[k],
                          **imshow_kwargs,
        )
        # beautify ...
        first_bin = -.5
        center_bin = hm_width/2 - .5
        last_bin = hm_width - .5
        flank_in_kb = int( (center_bin+.5)*binsizes[k]/1000 )
        flank_ticks = [first_bin, center_bin, last_bin]
        flank_ticklabels = [-flank_in_kb, 0, flank_in_kb]

        # average profile axes ...
        ax_profile[k].set_title(titles[k],fontsize=fig_fontsize)
        ax_profile[k].set_yscale(scales[k])
        ax_profile[k].set_ylim(vlims[k])
        ax_profile[k].tick_params(axis="y", length=0, direction="in", pad=-5)
        ax_profile[k].minorticks_off()
        ax_profile[k].set_yticks(vlims[k])
        ax_profile[k].set_yticklabels(vlims[k],fontsize=fig_fontsize)
        for _tidx, tick in enumerate(ax_profile[k].yaxis.get_majorticklabels()):
            tick.set_horizontalalignment("left")
            if _tidx == 0:
                tick.set_verticalalignment("bottom")
            elif _tidx == 1:
                tick.set_verticalalignment("top")
        ax_profile[k].tick_params(axis="x", length=6)
        ax_profile[k].set_xlim(first_bin, last_bin)
        ax_profile[k].set_xticks(flank_ticks)
        ax_profile[k].set_xticklabels(flank_ticklabels, fontsize=fig_fontsize)
        for _tidx, tick in enumerate(ax_profile[k].xaxis.get_majorticklabels()):
            if _tidx == 0:
                tick.set_horizontalalignment("left")
            elif _tidx == 2:
                tick.set_horizontalalignment("right")
            else:
                tick.set_horizontalalignment("center")

        # stackup axes controls ...
        ax_stackup[k].set_xticks(flank_ticks)
        ax_stackup[k].set_xticklabels(flank_ticklabels)
        ax_stackup[k].set_xticks(flank_ticks)
        ax_stackup[k].set_xticklabels(flank_ticklabels,fontsize=fig_fontsize)
        ax_stackup[k].tick_params(axis="x", length=6)
        ax_stackup[k].minorticks_off()
        ax_stackup[k].set_yticks([])
        ax_stackup[k].set_yticklabels([])
        for _tidx, tick in enumerate(ax_stackup[k].xaxis.get_majorticklabels()):
            if _tidx == 0:
                tick.set_horizontalalignment("left")
            elif _tidx == 2:
                tick.set_horizontalalignment("right")
            else:
                tick.set_horizontalalignment("center")

        # colorbar and itsd axes ...
        plt.colorbar(stack_hm, cax=ax_cbar[k], orientation="horizontal", ticks=vlims[k])
        ax_cbar[k].minorticks_off()
        ax_cbar[k].tick_params(axis="x", length=6)
        ax_cbar[k].set_xticklabels(vlims[k], fontsize=fig_fontsize)
        for _tidx, tick in enumerate(ax_cbar[k].xaxis.get_majorticklabels()):
            if _tidx == 0:
                tick.set_horizontalalignment("left")
            elif _tidx == 1:
                tick.set_horizontalalignment("right")


    return ax_xtra



def plot_stackups_sets(
    num_extra_plots,
    hms_dict_dict,  # heatmaps dict, that controls the order, the amount etc etc ...
    scales,
    vlims,
    titles,
    cmaps,
    binsizes,
    fillmissing=False,
    len_per_thousand=1.2,
    width_per_stack=3.5,
    extra_height=3.0,  # height that goes toward the profile and colorbar
    delta_h=.5,  # fixed distance between axes (vertically)
    delta_w=.2,  # fixed distance between axes (horizontally)
    # **plot_kwargs,
    fig_fontsize=25,
    **imshow_kwargs,
):
    """
    plot a buch of stackups ...
    """
    # rewrite everyhting assuming hms_dict_dict is a dict of stackup groups !
    # groups are plotted on top of each other ...

    if num_extra_plots:
        num_extra_plots = int(num_extra_plots)
    else:
        num_extra_plots = 0

    # inspect provided stacks and define figure with all of the panels ...
    num_stackup_groups, stackup_samples, stackup_group_heights, stack_width = _get_hms_nested_shape(hms_dict_dict)
    total_stack_height = sum(stackup_group_heights.values())
    num_cols = len(stackup_samples) + num_extra_plots
    stackup_height = total_stack_height*len_per_thousand / 1_000
    fig_height = stackup_height + extra_height
    fig_width = width_per_stack*num_cols
    # parameterize the hell out of it ...
    num_axes_vert = num_stackup_groups + 2  # + 2 is for average profile and a colorbar
    num_axes_horz = num_cols
    full_fig_width = fig_width + (num_axes_horz + 1)*delta_w
    full_fig_height = fig_height + (num_axes_vert + 1)*delta_h
    fig = plt.figure(
        figsize=(full_fig_width, full_fig_height),
        facecolor="white",
        layout="none",
    )
    gs = fig.add_gridspec(
        num_stackup_groups+2,
        num_cols,  # widsth ration are equal by default ...
        height_ratios = [
            0.95*extra_height,
            *[(_gh/total_stack_height)*stackup_height for _gh in stackup_group_heights.values()],
            0.05*extra_height,
        ],
        # horizontal spacings ...
        # left, right as a fraction of overall figure width
        left = delta_w/full_fig_width,
        right = 1. - delta_w/full_fig_width,
        # inter-axes spacing as a fraction of average axes height ...
        hspace = delta_h*num_axes_vert/fig_height,
        # vertical spacings ...
        # top, bottom as a fraction of overall figure height
        top = 1. - delta_h/full_fig_height,
        bottom = delta_h/full_fig_height,
        # inter-axes spacing as a fraction of average axes width ...
        wspace = delta_w*num_axes_horz/fig_width,
    )

    ##############
    ax_profile = {}
    ax_stackup = {}
    ax_xtra = []
    ax_cbar = {}
    # define nest dict of axes ...
    for jdx, sample in enumerate(stackup_samples):
        ax_profile[sample] = fig.add_subplot(gs[0, jdx])
        ax_cbar[sample] = fig.add_subplot(gs[-1, jdx])
        ax_stackup[sample] = {}
        for idx, group_k in enumerate(hms_dict_dict):
            ax_stackup[sample][group_k] = fig.add_subplot(gs[idx+1,jdx])
    # provide extra axes at the end ...
    for jdx in range(len(stackup_samples), num_cols):
        ax_xtra.append([fig.add_subplot(gs[_i+1,jdx]) for _i in range(num_stackup_groups)])

    # fill missing if needed and calculate profiles (per group) ...
    hms_dict_dict_copy = {}
    profile_hm = {}
    for group_k, hms_dict in hms_dict_dict.items():
        hms_dict_dict_copy[group_k] = _fillmissing_hms(hms_dict, how="col mean") if fillmissing else hms_dict
        # use modified stacks to calculate profiles ...
        profile_hm[group_k] = _get_profiles(hms_dict_dict_copy[group_k], scales)
    # replace hms_dict_dict with the updated copy ...
    hms_dict_dict = hms_dict_dict_copy
    # get norms - they are just per sample - regardless of the group ...
    norms = _get_norms(scales, vlims)

    last_group_k = list(hms_dict_dict.keys())[-1]
    first_sample = stackup_samples[0]

    # start plotting ...
    for group_k, hms_dict in hms_dict_dict.items():
        # we've checked that samples go in the same order ...
        for sample, hm in hms_dict.items():
            ax_profile[sample].plot(profile_hm[group_k][sample], linewidth=4)
            stack_hm = ax_stackup[sample][group_k].imshow(
                              hm,
                              norm=norms[sample],
                              aspect="auto",
                              cmap=cmaps[sample],
                              **imshow_kwargs,
            )
            # beautify ...
            ax_stackup[sample][group_k].set_xticks([])
            ax_stackup[sample][group_k].set_xticklabels([])
            ax_stackup[sample][group_k].set_yticks([])
            ax_stackup[sample][group_k].set_yticklabels([])
            ax_stackup[sample][group_k].minorticks_off()
            # #
            # if sample == first_sample:
            #     ax_stackup[sample][group_k].set_ylabel(group_k,fontsize=fig_fontsize)

            if group_k == last_group_k:
                # beautify ...
                # we have to do it for every samples - but not for every group ...
                first_bin = -.5
                center_bin = stack_width/2 - .5
                last_bin = stack_width - .5
                flank_in_kb = int( (center_bin+.5)*binsizes[sample]/1000 )
                flank_ticks = [first_bin, center_bin, last_bin]
                flank_ticklabels = [-flank_in_kb, "", flank_in_kb]
                #
                ax_profile[sample].set_title(titles[sample],fontsize=fig_fontsize)
                # ax_profile[sample].set_title(titles[sample])
                ax_profile[sample].minorticks_off()
                ax_profile[sample].set_xlim([first_bin, last_bin])
                ax_profile[sample].set_xticks(flank_ticks)
                ax_profile[sample].tick_params(axis="x", length=6)
                ax_profile[sample].set_xticklabels(flank_ticklabels,fontsize=fig_fontsize)
                for _tidx, tick in enumerate(ax_profile[sample].xaxis.get_majorticklabels()):
                    if _tidx == 0:
                        tick.set_horizontalalignment("left")
                    elif _tidx == 2:
                        tick.set_horizontalalignment("right")
                    else:
                        tick.set_horizontalalignment("center")
                ax_profile[sample].set_ylim(vlims[sample])
                ax_profile[sample].tick_params(axis="y", length=0, direction="in", pad=-5)
                ax_profile[sample].set_yticks(vlims[sample])
                ax_profile[sample].set_yticklabels(vlims[sample],fontsize=fig_fontsize)
                for _tidx, tick in enumerate(ax_profile[sample].yaxis.get_majorticklabels()):
                    tick.set_horizontalalignment("left")
                    if _tidx == 0:
                        tick.set_verticalalignment("bottom")
                    elif _tidx == 1:
                        tick.set_verticalalignment("top")

                # bottom one - show ticks for now ...
                ax_stackup[sample][group_k].set_xticks(flank_ticks)
                ax_stackup[sample][group_k].set_xticklabels(flank_ticklabels,fontsize=fig_fontsize)
                ax_stackup[sample][group_k].tick_params(axis="x", length=6)        
                ax_stackup[sample][group_k].set_yticks([])
                ax_stackup[sample][group_k].set_yticklabels([])
                for _tidx, tick in enumerate(ax_stackup[sample][group_k].xaxis.get_majorticklabels()):
                    if _tidx == 0:
                        tick.set_horizontalalignment("left")
                    elif _tidx == 2:
                        tick.set_horizontalalignment("right")
                    else:
                        tick.set_horizontalalignment("center")
                # colorbar - draw them one time per sample only !
                plt.colorbar(stack_hm,
                            cax=ax_cbar[sample],
                            orientation="horizontal",
                            ticks=vlims[sample])
                ax_cbar[sample].minorticks_off()
                ax_cbar[sample].tick_params(axis="x", length=6)
                ax_cbar[sample].set_xticklabels(vlims[sample],fontsize=fig_fontsize)
                for _tidx, tick in enumerate(ax_cbar[sample].xaxis.get_majorticklabels()):
                    if _tidx == 0:
                        tick.set_horizontalalignment("left")
                    elif _tidx == 1:
                        tick.set_horizontalalignment("right")

    return ax_xtra


# profile_height=0.35
# margin_h=0.2
# margin_v=0.2
# spacing_h=0.02
# spacing_v=0.15
# cbarh_spacing_v = 0.05
# fig_fontsize=6
# cbarh = 0.1
# from helper_func import _get_norms, _get_hms_nested_shape, _fillmissing_hms, _get_profiles

def plot_stackups_sets_new(
    num_extra_plots,
    hms_dict_dict,  # heatmaps dict, that controls the order, the amount etc etc ...
    scales,
    vlims,
    titles,
    cmaps,
    binsizes,
    fillmissing=False,
    extra_plots_position="left",
    len_per_thousand=0.75,
    width_per_stack=0.35,
    profile_height=0.35,
    cbar_height = 0.1,
    spacing_v=.5,  # fixed distance between axes (vertically)
    spacing_h=.2,  # fixed distance between axes (horizontally)
    # **plot_kwargs,
    fig_fontsize=6,
    **imshow_kwargs,
):
    """
    plot a buch of stackups ...
    """
    # rewrite everyhting assuming hms_dict_dict is a dict of stackup groups !
    # groups are plotted on top of each other ...

    if num_extra_plots:
        num_extra_plots = int(num_extra_plots)
    else:
        num_extra_plots = 0

    # inspect provided stacks and define figure with all of the panels ...
    num_stackup_groups, stackup_samples, stackup_group_heights, stack_width = _get_hms_nested_shape(hms_dict_dict)
    num_cols = len(stackup_samples) + num_extra_plots

    # in inches
    margin_h=0.2
    margin_v=0.2
    profile_spacing_v = 0.075
    cbarh_spacing_v = 0.05
    profile_color = "dimgray"
    imshow_kwargs = dict(interpolation="antialiased", interpolation_stage="data", filternorm=False)
    # imshow_kwargs = dict(interpolation="antialiased", interpolation_stage="rgba", filternorm=True)

    # horizontal splitting layout
    h_split = []
    h_split.append( Size.Fixed(margin_h) )
    h_split += [Size.Fixed(width_per_stack), Size.Fixed(spacing_h)]*(num_cols-1)
    h_split += [Size.Fixed(width_per_stack), Size.Fixed(margin_h)]

    # vertical splitting layout
    v_split = []
    v_split.append( Size.Fixed(margin_v) )
    v_split.append( Size.Fixed(cbar_height) )
    v_split.append( Size.Fixed(cbarh_spacing_v) )
    for _i, _num_row_per_stack in enumerate(stackup_group_heights.values()):
        v_split.append( Size.Fixed(len_per_thousand*(_num_row_per_stack/1_000)) )
        if _i < len(stackup_group_heights)-1:
            v_split.append( Size.Fixed(spacing_v) )
        else:
            v_split.append( Size.Fixed(spacing_v+profile_spacing_v) )
    profile_spacing_v
    v_split.append( Size.Fixed(profile_height) )
    v_split.append( Size.Fixed(margin_v) )


    # set figsize based on the tiling provided - i.e. post factum ...
    fig_width = sum(_h.fixed_size for _h in h_split)
    fig_height = sum(_v.fixed_size for _v in v_split)
    fig = plt.figure(
        figsize=(fig_width, fig_height),
        layout="none",
        # facecolor='lightblue'
    )
    print(f"figure overall is {fig_width=} {fig_height=}")

    divider = Divider(fig, (0, 0, 1, 1), h_split, v_split, aspect=False)
    _div_pos = divider.get_position()


    ax_profile = {}
    ax_stackup = {}
    ax_xtra = []
    ax_cbar = {}
    # define nest dict of axes ...
    # provide extra axes at the end ...
    if extra_plots_position == "left":
        # extra plots on the left ...
        for jdx in range(num_extra_plots):
            for idx, group_k in enumerate(hms_dict_dict):
                # ax_xtra.append([fig.add_subplot(gs[_i+1,jdx]) for _i in range(num_stackup_groups)])
                idx += 1  # adjust by 1, since there is a cbar at the bottom
                _stack_group_locator = divider.new_locator(nx=2*jdx+1, ny=2*idx+1)
                ax_xtra.append(fig.add_axes(_div_pos, axes_locator=_stack_group_locator))
        for jdx, sample in enumerate(stackup_samples):
            jdx += num_extra_plots  # adjust steps by the extract plots in front
            _cbar_locator = divider.new_locator(nx=2*jdx+1, ny=1)
            ax_cbar[sample] = fig.add_axes(_div_pos, axes_locator=_cbar_locator)
            ax_stackup[sample] = {}
            for idx, group_k in enumerate(hms_dict_dict):
                idx += 1  # adjust by 1, since there is a cbar at the bottom
                _stack_group_locator = divider.new_locator(nx=2*jdx+1, ny=2*idx+1)
                ax_stackup[sample][group_k] = fig.add_axes(_div_pos, axes_locator=_stack_group_locator)
            # profile ny is simply the next one :
            _profile_locator = divider.new_locator(nx=2*jdx+1, ny=2*(idx+1)+1)
            ax_profile[sample] = fig.add_axes(_div_pos, axes_locator=_profile_locator)
    # # provide extra axes at the end ...
    # if extra_plots_position == "right":
    else:  # RIGHT ...
        # start with the stacks
        for jdx, sample in enumerate(stackup_samples):
            _cbar_locator = divider.new_locator(nx=2*jdx+1, ny=1)
            ax_cbar[sample] = fig.add_axes(_div_pos, axes_locator=_cbar_locator)
            ax_stackup[sample] = {}
            for idx, group_k in enumerate(hms_dict_dict):
                idx += 1  # adjust by 1, since there is a cbar at the bottom
                _stack_group_locator = divider.new_locator(nx=2*jdx+1, ny=2*idx+1)
                ax_stackup[sample][group_k] = fig.add_axes(_div_pos, axes_locator=_stack_group_locator)
            # profile ny is simply the next one :
            _profile_locator = divider.new_locator(nx=2*jdx+1, ny=2*(idx+1)+1)
            ax_profile[sample] = fig.add_axes(_div_pos, axes_locator=_profile_locator)
        # add extra plots in the end (on the right) ...
        for jdx in range(len(stackup_samples), len(stackup_samples)+num_extra_plots):
            for idx, group_k in enumerate(hms_dict_dict):
                idx += 1  # adjust by 1, since there is a cbar at the bottom
                _stack_group_locator = divider.new_locator(nx=2*jdx+1, ny=2*idx+1)
                ax_xtra.append(fig.add_axes(_div_pos, axes_locator=_stack_group_locator))

    # fill missing if needed and calculate profiles (per group) ...
    hms_dict_dict_copy = {}
    profile_hm = {}
    for group_k, hms_dict in hms_dict_dict.items():
        hms_dict_dict_copy[group_k] = _fillmissing_hms(hms_dict, how="col mean") if fillmissing else hms_dict
        # use modified stacks to calculate profiles ...
        profile_hm[group_k] = _get_profiles(hms_dict_dict_copy[group_k], scales)
    # replace hms_dict_dict with the updated copy ...
    hms_dict_dict = hms_dict_dict_copy
    # get norms - they are just per sample - regardless of the group ...
    norms = {k: _get_norm(scales[k], vlims[k]) for k in stackup_samples}

    last_group_k = list(hms_dict_dict.keys())[-1]
    first_sample = stackup_samples[0]

    # start plotting ...
    for group_k, hms_dict in hms_dict_dict.items():
        # we've checked that samples go in the same order ...
        for sample, hm in hms_dict.items():
            ax_profile[sample].plot(profile_hm[group_k][sample], linewidth=1)
            stack_hm = ax_stackup[sample][group_k].imshow(
                              hm,
                              norm=norms[sample],
                              aspect="auto",
                              cmap=cmaps[sample],
                              **imshow_kwargs,
            )
            # beautify ...
            ax_stackup[sample][group_k].set_xticks([])
            ax_stackup[sample][group_k].set_xticklabels([])
            ax_stackup[sample][group_k].set_yticks([])
            ax_stackup[sample][group_k].set_yticklabels([])
            ax_stackup[sample][group_k].minorticks_off()
            # #
            # if sample == first_sample:
            #     ax_stackup[sample][group_k].set_ylabel(group_k,fontsize=fig_fontsize)

            if group_k == last_group_k:
                # beautify ...
                # we have to do it for every samples - but not for every group ...
                first_bin = -.5
                center_bin = stack_width/2 - .5
                last_bin = stack_width - .5
                _xticklength = 1
                flank_in_kb = int( (center_bin+.5)*binsizes[sample]/1000 )
                flank_ticks = [first_bin, center_bin, last_bin]
                flank_ticklabels = [-flank_in_kb, "", flank_in_kb]
                #
                ax_profile[sample].set_title(titles[sample],fontsize=fig_fontsize, pad=2.5)
                # ax_profile[sample].set_title(titles[sample])
                ax_profile[sample].minorticks_off()
                ax_profile[sample].set_xlim([first_bin, last_bin])
                ax_profile[sample].set_xticks(flank_ticks)
                ax_profile[sample].tick_params(axis="x", length=_xticklength, pad=0.5)
                ax_profile[sample].set_xticklabels(flank_ticklabels,fontsize=fig_fontsize)
                for _tidx, tick in enumerate(ax_profile[sample].xaxis.get_majorticklabels()):
                    if _tidx == 0:
                        tick.set_horizontalalignment("left")
                    elif _tidx == 2:
                        tick.set_horizontalalignment("right")
                    else:
                        tick.set_horizontalalignment("center")
                ax_profile[sample].set_ylim(vlims[sample])
                ax_profile[sample].tick_params(axis="y", length=0, direction="in", pad=-5)
                ax_profile[sample].set_yticks(vlims[sample])
                ax_profile[sample].set_yticklabels(vlims[sample],fontsize=fig_fontsize)
                for _tidx, tick in enumerate(ax_profile[sample].yaxis.get_majorticklabels()):
                    tick.set_horizontalalignment("left")
                    if _tidx == 0:
                        tick.set_verticalalignment("bottom")
                    elif _tidx == 1:
                        tick.set_verticalalignment("top")

                # # bottom one - show ticks for now ...
                # ax_stackup[sample][group_k].set_xticks(flank_ticks)
                # ax_stackup[sample][group_k].set_xticklabels(flank_ticklabels,fontsize=fig_fontsize)
                # ax_stackup[sample][group_k].tick_params(axis="x", length=6)
                # ax_stackup[sample][group_k].set_yticks([])
                # ax_stackup[sample][group_k].set_yticklabels([])
                # for _tidx, tick in enumerate(ax_stackup[sample][group_k].xaxis.get_majorticklabels()):
                #     if _tidx == 0:
                #         tick.set_horizontalalignment("left")
                #     elif _tidx == 2:
                #         tick.set_horizontalalignment("right")
                #     else:
                #         tick.set_horizontalalignment("center")
                # # colorbar - draw them one time per sample only !
                plt.colorbar(stack_hm,
                            cax=ax_cbar[sample],
                            orientation="horizontal",
                            ticks=vlims[sample])
                ax_cbar[sample].minorticks_off()
                ax_cbar[sample].tick_params(axis="x", length=_xticklength, pad=0.5)
                ax_cbar[sample].set_xticklabels(vlims[sample],fontsize=fig_fontsize)
                for _tidx, tick in enumerate(ax_cbar[sample].xaxis.get_majorticklabels()):
                    if _tidx == 0:
                        tick.set_horizontalalignment("left")
                    elif _tidx == 1:
                        tick.set_horizontalalignment("right")

    return ax_xtra