#
# graphics/uq/double_couple.py - uncertainty quantification of double couple sources
#

import numpy as np

from matplotlib import pyplot
from pandas import DataFrame
from xarray import DataArray
from mtuq.graphics._gmt import read_cpt
from mtuq.grid_search import MTUQDataArray, MTUQDataFrame
from mtuq.util import fullpath, warn
from mtuq.util.math import closed_interval, open_interval, to_delta, to_gamma
from os.path import exists


def plot_misfit_dc(filename, ds, title='',
    colorbar_type=1, marker_type=1):
    """ Plots misfit over strike, dip, and slip
    (matplotlib implementation)
    """
    _check(ds)
    ds = ds.copy()

    if issubclass(type(ds), DataArray):
        _plot_dc(filename, _squeeze(ds), cmap='viridis')
        
    elif issubclass(type(ds), DataFrame):
        warn('plot_misfit_dc not implemented for irregularly-spaced grids')


def plot_likelihood_dc(filename, ds, sigma=None, title=''):
    assert sigma is not None

    _check(ds)
    ds = ds.copy()

    if issubclass(type(ds), DataArray):
        ds.values = np.exp(-ds.values/(2.*sigma**2))
        _plot_dc(filename, _squeeze(ds), cmap='hot',
                 colorbar_type=colorbar_type, marker_type=marker_type)

    elif issubclass(type(ds), DataFrame):
        warn('plot_likelihood_dc not implemented for irregularly-spaced grids')


def plot_marginal_dc():
    raise NotImplementedError


def _squeeze(da):
    if 'origin_idx' in da.dims:
        da = da.max(dim='origin_idx')

    if 'rho' in da.dims:
        da = da.max(dim='rho')

    if 'v' in da.dims:
        assert len(da.coords['v'])==1
        da = da.squeeze(dim='v')

    if 'w' in da.dims:
        assert len(da.coords['w'])==1
        da = da.squeeze(dim='w')

    return da


def _check(ds):
    """ Checks data structures
    """
    if type(ds) not in (DataArray, DataFrame, MTUQDataArray, MTUQDataFrame):
        raise TypeError("Unexpected grid format")


#
# matplotlib backend
#

def _plot_dc(filename, da, colorbar_type=1, marker_type=1, cmap='hot', **kwargs):
    # FIXME: do labels correspond to the correct axes ?!

    # prepare axes
    fig, axes = pyplot.subplots(2, 2, 
        figsize=(8., 8.),
        )

    pyplot.subplots_adjust(
        wspace=0.4,
        hspace=0.4,
        )

    if exists(_local_path(cmap)):
       cmap = read_cpt(_local_path(cmap))

    # upper left panel
    marginal = da.min(dim=('sigma'))
    x = marginal.coords['h']
    y = marginal.coords['kappa']

    minmax1 = _minmax(x, y, marginal)

    axis = axes[0][0]

    _pcolor(axis, x, y, marginal.values, cmap, **kwargs)

    axis.set_xlabel('Dip', **axis_label_kwargs)
    axis.set_xticks(theta_ticks)
    axis.set_xticklabels(theta_ticklabels)

    axis.set_ylabel('Strike', **axis_label_kwargs)
    axis.set_yticks(kappa_ticks)
    axis.set_yticklabels(kappa_ticklabels)

    # upper right panel
    marginal = da.min(dim=('h'))
    x = marginal.coords['sigma']
    y = marginal.coords['kappa']

    minmax2 = _minmax(x, y, marginal)

    axis = axes[0][1]

    _pcolor(axis, x, y, marginal.values, cmap, **kwargs)

    axis.set_xlabel('Slip', **axis_label_kwargs)
    axis.set_xticks(sigma_ticks)
    axis.set_xticklabels(sigma_ticklabels)

    axis.set_ylabel('Strike', **axis_label_kwargs)
    axis.set_yticks(kappa_ticks)
    axis.set_yticklabels(kappa_ticklabels)

    # lower right panel
    marginal = da.min(dim=('kappa'))
    y = marginal.coords['h']
    x = marginal.coords['sigma']

    minmax3 = _minmax(x, y, marginal.T)

    axis = axes[1][1]

    _pcolor(axis, x, y, marginal.values.T, cmap, **kwargs)

    axis.set_xlabel('Slip', **axis_label_kwargs)
    axis.set_xticks(sigma_ticks)
    axis.set_xticklabels(sigma_ticklabels)

    axis.set_ylabel('Dip', **axis_label_kwargs)
    axis.set_yticks(theta_ticks)
    axis.set_yticklabels(theta_ticklabels)

    # lower left panel
    axes[1][0].axis('off')

    # optional markers
    if marker_type > 0:
        _add_marker(axes[0][0], minmax1[marker_type-1])
        _add_marker(axes[0][1], minmax2[marker_type-1])
        _add_marker(axes[1][1], minmax3[marker_type-1])

    pyplot.savefig(filename)


def _add_marker(axis, coords):
    axis.scatter(*coords, s=250,
        marker='o',
        facecolors='none',
        edgecolors=[0,1,0],
        linewidths=1.75,
        clip_on=False,
        zorder=100,
        )


def _minmax(x, y, values):
    x = x.values
    y = y.values
    iymin, ixmin = np.unravel_index(values.argmin(), values.shape)
    iymax, ixmax = np.unravel_index(values.argmax(), values.shape)
    xmin, ymin = x[ixmin], y[iymin]
    xmax, ymax = x[ixmax], y[iymax]
    return (xmin, ymin), (xmax, ymax)


axis_label_kwargs = {
    'fontsize': 14
    }


def _pcolor(axis, x, y, values, cmap, **kwargs):
    # workaround matplotlib compatibility issue
    try:
        axis.pcolor(x, y, values, cmap=cmap, shading='auto',  **kwargs)
    except:
        axis.pcolor(x, y, values, cmap=cmap, **kwargs)


kappa_ticks = [0, 45, 90, 135, 180, 225, 270, 315, 360]
kappa_ticklabels = ['0', '', '90', '', '180', '', '270', '', '360']

sigma_ticks = [-90, -67.5, -45, -22.5, 0, 22.5, 45, 67.5, 90]
sigma_ticklabels = ['-90', '', '-45', '', '0', '', '45', '', '90']

theta_ticks = [np.cos(np.radians(tick)) for tick in [0, 15, 30, 45, 60, 75, 90]]
theta_ticklabels = ['0', '', '30', '', '60', '', '90']


def _local_path(name):
    return fullpath('mtuq/graphics/_gmt/cpt', name+'.cpt')


