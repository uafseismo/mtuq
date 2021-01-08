
#
# graphics/uq/vw.py - uncertainty quantification on the v-w rectangle
#

import numpy as np

from matplotlib import pyplot
from pandas import DataFrame
from xarray import DataArray
from mtuq.graphics._gmt import read_cpt
from mtuq.graphics.uq._gmt import _nothing_to_plot
from mtuq.grid_search import MTUQDataArray, MTUQDataFrame
from mtuq.util import fullpath
from mtuq.util.math import closed_interval, open_interval
from os.path import exists


#
# For details about the  rectangle:`v-w` rectangle rectangle, see 
# Tape2015 - A uniform parameterization of moment tensors
# (https://doi.org/10.1093/gji/ggv262)
#

v_min = -1./3.
v_max = +1./3.
w_min = -3.*np.pi/8.
w_max = +3.*np.pi/8.
vw_area = (v_max-v_min)*(w_max-w_min)



def plot_misfit_vw(filename, ds, misfit_callback=None, title=''):
    """ Plots misfit values on `v-w` rectangle


    .. rubric :: Input arguments

    ``filename`` (`str`):
    Name of output image file

    ``ds`` (`DataArray` or `DataFrame`):
    data structure containing moment tensors and corresponding misfit values

    ``misfit_callback`` (func)
    User-supplied function applied to misfit values

    ``title`` (`str`):
    Optional figure title

    """
    _check(ds)
    ds = ds.copy()

    if issubclass(type(ds), DataArray):
        ds = ds.min(dim=('origin_idx', 'rho', 'kappa', 'sigma', 'h'))
        v = ds.coords['v']
        w = ds.coords['w']
        values = ds.values.transpose()

    elif issubclass(type(ds), DataFrame):
        ds = ds.reset_index()
        v, w, values = _bin(ds, lambda ds: ds.min())

    if misfit_callback:
        values = misfit_callback(values)

    _plot_misfit_vw(filename, v, w, values, title=title)


def plot_likelihood_vw(filename, ds, sigma=None, title=''):
    """ Plots maximum likelihoods on `v-w` rectangle


    .. rubric :: Input arguments

    ``filename`` (`str`):
    Name of output image file

    ``ds`` (`DataArray` or `DataFrame`):
    data structure containing moment tensors and corresponding misfit values

    ``title`` (`str`):
    Optional figure title

    """
    assert sigma is not None

    _check(ds)
    ds = ds.copy()

    if issubclass(type(ds), DataArray):
        ds.values = np.exp(-ds.values/(2.*sigma**2))
        ds.values /= ds.values.sum()
        ds = ds.max(dim=('origin_idx', 'rho', 'kappa', 'sigma', 'h'))
        v = ds.coords['v']
        w = ds.coords['w']
        values = ds.values.transpose()

    elif issubclass(type(ds), DataFrame):
        ds = np.exp(-ds/(2.*sigma**2))
        ds /= ds.sum()
        ds = ds.reset_index()
        v, w, values = _bin(ds, lambda ds: ds.max())

    values /= values.sum()
    values /= vw_area

    _plot_likelihood_vw(filename, v, w, values, title=title)


def plot_marginal_vw(filename, ds, sigma=None, title=''):
    """ Plots marginal likelihoods on `v-w` rectangle


    .. rubric :: Input arguments

    ``filename`` (`str`):
    Name of output image file

    ``ds`` (`DataArray` or `DataFrame`):
    data structure containing moment tensors and corresponding misfit values

    ``title`` (`str`):
    Optional figure title

    """
    assert sigma is not None

    _check(ds)
    ds = ds.copy()

    if issubclass(type(ds), DataArray):
        ds.values = np.exp(-ds.values/(2.*sigma**2))
        ds.values /= ds.values.sum()
        ds = ds.max(dim=('origin_idx', 'rho', 'kappa', 'sigma', 'h'))
        v = ds.coords['v']
        w = ds.coords['w']
        values = ds.values.transpose()

    elif issubclass(type(ds), DataFrame):
        ds = np.exp(-ds/(2.*sigma**2))
        ds /= ds.sum()
        ds = ds.reset_index()
        v, w, values = _bin(ds, lambda ds: ds.sum()/len(ds))

    values /= values.sum()
    values /= vw_area

    _plot_likelihood_vw(filename, v, w, values, title=title)



#
# pyplot wrappers
#

def _plot_misfit_vw(filename, v, w, values,
    colorbar_type=1, marker_type=1, title=''):

    if _nothing_to_plot(values):
        return

    _plot_vw(v, w, values, 
        colorbar_type=colorbar_type,
        cmap='hot',
        title=title)

    if marker_type:
        idx = np.unravel_index(values.argmin(), values.shape)
        coords = v[idx[1]], w[idx[0]]

        pyplot.scatter(*coords, s=333,
            marker='o',
            facecolors='none',
            edgecolors=[0,1,0],
            linewidths=1.75,
            )

    pyplot.savefig(filename)
    pyplot.close()


def _plot_likelihood_vw(filename, v, w, values,
    colorbar_type=1, marker_type=2, title=''):

    if _nothing_to_plot(values):
        return

    _plot_vw(v, w, values, 
        colorbar_type=colorbar_type,
        cmap='hot_r',
        title=title)

    if marker_type:
        idx = np.unravel_index(values.argmax(), values.shape)
        coords = v[idx[1]], w[idx[0]]

        pyplot.scatter(*coords, s=333,
            marker='o', 
            facecolors='none',
            edgecolors=[0,1,0],
            linewidths=1.75,
            clip_on=False,
            zorder=100,
            )

    pyplot.savefig(filename)


def _plot_vw(v, w, values, colorbar_type=0, cmap='hot_r', title=None):
    # create figure
    fig, ax = pyplot.subplots(figsize=(3., 8.), constrained_layout=True)

    # pcolor requires corners of pixels
    corners_v = _centers_to_edges(v)
    corners_w = _centers_to_edges(w)

    # `values` gets mapped to pixel colors
    pyplot.pcolor(corners_v, corners_w, values, cmap=cmap)

    # v and w have the following bounds
    # (see https://doi.org/10.1093/gji/ggv262)
    pyplot.xlim([-1./3., 1./3.])
    pyplot.ylim([-3./8.*np.pi, 3./8.*np.pi])

    pyplot.xticks([], [])
    pyplot.yticks([], [])

    if exists(_local_path(cmap)):
       cmap = read_cpt(_local_path(cmap))

    if colorbar_type:
        cbar = pyplot.colorbar(
            orientation='horizontal',
            pad=0.,
            )

        cbar.formatter.set_powerlimits((-2, 2))

    if title:
        fontdict = {'fontsize': 16}
        pyplot.title(title, fontdict=fontdict)


def _check(ds):
    """ Checks data structures
    """
    if type(ds) not in (DataArray, DataFrame, MTUQDataArray, MTUQDataFrame):
        raise TypeError("Unexpected grid format")


def _bin(df, handle, npts_v=20, npts_w=40):
    """ Bins DataFrame into rectangular cells
    """
    # define centers of cells
    centers_v = open_interval(-1./3., 1./3., npts_v)
    centers_w = open_interval(-3./8.*np.pi, 3./8.*np.pi, npts_w)

    # define corners of cells
    v = closed_interval(-1./3., 1./3., npts_v+1)
    w = closed_interval(-3./8.*np.pi, 3./8.*np.pi, npts_w+1)

    binned = np.empty((npts_w, npts_v))
    for _i in range(npts_w):
        for _j in range(npts_v):
            # which grid points lie within cell (i,j)?
            subset = df.loc[
                df['v'].between(v[_j], v[_j+1]) &
                df['w'].between(w[_i], w[_i+1])]

            binned[_i, _j] = handle(subset[0])

    return centers_v, centers_w, binned


def _centers_to_edges(v):
    if issubclass(type(v), DataArray):
        v = v.values.copy()
    else:
        v = v.copy()

    dv = (v[1]-v[0])
    v -= dv/2
    v = np.pad(v, (0, 1))
    v[-1] = v[-2] + dv

    return v


def _local_path(name):
    return fullpath('mtuq/graphics/_gmt/cpt', name+'.cpt')


