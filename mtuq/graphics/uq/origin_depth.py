
#
# graphics/uq/origin.py - uncertainty quantification of source origin
#

import numpy as np
import subprocess

from matplotlib import pyplot
from pandas import DataFrame
from xarray import DataArray
from mtuq.graphics.uq._gmt import exists_gmt, gmt_not_found_warning,\
    _parse_filetype, _parse_title
from mtuq.grid_search import MTUQDataArray, MTUQDataFrame
from mtuq.util import fullpath, warn
from mtuq.util.math import closed_interval, open_interval



def plot_misfit_depth(filename, ds, origins, sources, misfit_callback=None, title=''):
    """ Plots misfit versus depth


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
    depths = _get_depths(origins)

    _check(ds)
    ds = ds.copy()

    if issubclass(type(ds), DataArray):
        values, indices = _min_dataarray(ds)
        best_sources = _get_sources(sources, indices)

    elif issubclass(type(ds), DataFrame):
        values, indices = _min_dataframe(ds)
        best_sources = _get_sources(sources, indices)

    if misfit_callback:
        values = misfit_callback(values)

    _plot_depth(filename, depths, values, indices,
        title, xlabel='auto', ylabel='Misfit')



def plot_likelihood_depth(filename, ds, origins, sources, sigma=None, title=''):
    """ Plots marginal likelihood versus depth


    .. rubric :: Input arguments

    ``filename`` (`str`):
    Name of output image file

    ``ds`` (`DataArray` or `DataFrame`):
    data structure containing moment tensors and corresponding misfit values

    ``title`` (`str`):
    Optional figure title

    """
    assert sigma is not None

    depths = _get_depths(origins)

    _check(ds)
    ds = ds.copy()

    if issubclass(type(ds), DataArray):
        ds.values = np.exp(-ds.values/(2.*sigma**2))
        ds.values /= ds.values.sum()

        values, indices = _min_dataarray(ds)
        best_sources = _get_sources(sources, indices)

    elif issubclass(type(ds), DataFrame):
        ds = np.exp(-ds/(2.*sigma**2))
        ds /= ds.sum()

        values, indices = _min_dataframe(ds)
        best_sources = _get_sources(sources, indices)

    values /= values.sum()

    _plot_depth(filename, depths, values, indices, 
        title=title, xlabel='auto', ylabel='Likelihood')



def plot_marginal_depth(filename, ds, origins, sources, sigma=None, title=''):
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

    depths = _get_depths(origins)

    _check(ds)
    ds = ds.copy()

    if issubclass(type(ds), DataArray):
        ds = np.exp(-ds/(2.*sigma**2))
        ds /= ds.sum()

        values, indices = _max_dataarray(ds)
        best_sources = _get_sources(sources, indices)

    elif issubclass(type(ds), DataFrame):
        raise NotImplementedError
        ds = np.exp(-ds/(2.*sigma**2))
        ds /= ds.sum()

        values, indices = _min_dataframe(ds)
        best_sources = _get_sources(sources, indices)

    values /= values.sum()

    _plot_depth(filename, depths, values, indices, 
        title=title, xlabel='auto', ylabel='Likelihood')



#
# utility functions
#

def _check(ds):
    """ Checks data structures
    """
    if type(ds) not in (DataArray, DataFrame, MTUQDataArray, MTUQDataFrame):
        raise TypeError("Unexpected grid format")


def _get_depths(origins):
    depths = []
    for origin in origins:
        depths += [float(origin.depth_in_m)]
    return np.array(depths)


def _get_sources(sources, indices):
    return [sources.get(index) for index in indices]


def _min_dataarray(ds):
    values, indices = [], []
    for _i in range(ds.shape[-1]):
        sliced = ds[:,:,:,:,:,:,_i]
        values += [sliced.min()]
        indices += [int(sliced.argmin())]
    return np.array(values), indices


def _max_dataarray(ds):
    values, indices = [], []
    for _i in range(ds.shape[-1]):
        sliced = ds[:,:,:,:,:,:,_i]
        values += [sliced.max()]
        indices += [int(sliced.argmax())]
    return np.array(values), indices


def _sum_dataarray(ds):
    raise NotImplementedError

def _min_dataframe(ds):
    raise NotImplementedError

def _max_dataframe(ds):
    raise NotImplementedError

def _sum_dataframe(ds):
    raise NotImplementedError


#
# pyplot wrappers
#

def _plot_depth(filename, depths, values, best_sources, title='',
    xlabel='auto', ylabel='', show_magnitudes=False, show_beachballs=False,
    fontsize=16., normalize=False):

    if xlabel=='auto' and ((depths.max()-depths.min()) < 10000.):
       xlabel = 'Depth (m)'

    if xlabel=='auto' and ((depths.max()-depths.min()) >= 10000.):
       depths /= 1000.
       xlabel = 'Depth (km)'

    if normalize:
        values /= values.max()

    figsize = (6., 6.)
    pyplot.figure(figsize=figsize)

    pyplot.plot(depths, values, 'k-')

    if show_magnitudes:
        raise NotImplementedError

    if show_beachballs:
        raise NotImplementedError

    if xlabel:
         pyplot.xlabel(xlabel, fontsize=fontsize)

    if ylabel:
         pyplot.ylabel(ylabel, fontsize=fontsize)

    if title:
        pyplot.title(title, fontsize=fontsize)

    pyplot.savefig(filename)


