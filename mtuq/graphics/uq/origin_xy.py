
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


def plot_misfit_xy(filename, ds, origins, sources, title='', 
    labeltype='latlon', colorbar_type=0, marker_type=1):
    """ Plots misfit versus hypocenter position


    .. rubric :: Input arguments

    ``filename`` (`str`):
    Name of output image file

    ``ds`` (`DataArray` or `DataFrame`):
    data structure containing moment tensors and corresponding misfit values

    ``title`` (`str`):
    Optional figure title

    """
    x, y = _get_xy(origins)

    _check(ds)
    ds = ds.copy()

    if issubclass(type(ds), DataArray):
        values, indices = _min_dataarray(ds)
        best_sources = _get_sources(sources, indices)

    elif issubclass(type(ds), DataFrame):
        values, indices = _min_dataframe(ds)
        best_sources = _get_sources(sources, indices)

    _plot_misfit_xy(filename, x, y, values, title, labeltype)


def plot_mt_xy(filename, ds, origins, sources, title='', labeltype='latlon'):
    """ Plots focal mechanism versus hypocenter position


    .. rubric :: Input arguments

    ``filename`` (`str`):
    Name of output image file

    ``ds`` (`DataArray` or `DataFrame`):
    data structure containing moment tensors and corresponding misfit values

    ``title`` (`str`):
    Optional figure title

    """
    x, y = _get_xy(origins)

    _check(ds)
    ds = ds.copy()

    if issubclass(type(ds), DataArray):
        values, indices = _min_dataarray(ds)
        best_sources = _get_sources(sources, indices)

    elif issubclass(type(ds), DataFrame):
        values, indices = _min_dataframe(ds)
        best_sources = _get_sources(sources, indices)

    _plot_mt_xy_gmt(filename, x, y, best_sources, title, labeltype)


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


def _get_xy(origins):
    x, y = [], []
    for origin in origins:
        x += [float(origin.offset_x_in_m)]
        y += [float(origin.offset_y_in_m)]
    return np.array(x), np.array(y)


def _get_sources(sources, indices):
    return [sources.get(index) for index in indices]


def _min_dataarray(ds):
    values, indices = [], []
    for _i in range(ds.shape[-1]):
        sliced = ds[:,:,:,:,:,:,_i]#.squeeze()
        values += [sliced.values.min()]
        indices += [int(sliced.values.argmin())]
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

def _plot_misfit_xy(filename, x, y, values, title='', labeltype='latlon',
    colorbar_type=0, marker_type=1, cmap='hot'):

    xlabel, ylabel = _get_labeltype(x, y, labeltype)

    assert len(x)==len(y)==len(values), ValueError

    ux = np.unique(x)
    uy = np.unique(y)
    if len(ux)*len(uy)!=len(values):
        warn('Irregular x-y misfit grid')

    figsize = (6., 6.)
    pyplot.figure(figsize=figsize)

    pyplot.tricontourf(x, y, values, 100, cmap=cmap)

    if marker_type:
        idx = values.argmin()
        coords = x[idx], y[idx]

        pyplot.scatter(*coords, s=250,
            marker='o',
            facecolors='none',
            edgecolors=[0,1,0],
            linewidths=1.75,
            )

    if xlabel:
         pyplot.xlabel(xlabel)

    if ylabel:
         pyplot.ylabel(ylabel)

    if title:
        pyplot.title(title)

    pyplot.gca().axis('square')

    pyplot.savefig(filename)


#
# gmt wrappers
#

def _plot_mt_xy_gmt(filename, x, y, sources, title='',
    labeltype='latlon', show_magnitudes=False, show_beachballs=True):

    filetype = _parse_filetype(filename)
    title, subtitle = _parse_title(title)
    xlabel, ylabel = _get_labeltype(x, y, labeltype)

    xmin, xmax, ymin, ymax = _get_limits(x, y) 


    assert len(x)==len(y)==len(sources), ValueError

    ux = np.unique(x)
    uy = np.unique(y)
    if len(ux)*len(uy)!=len(sources):
        warn('Irregular x-y grid')

    mw_array = None
    if show_magnitudes:
        mw_array = np.zeros((len(sources), 3))
        for _i, source in enumerate(sources):
            mw_array[_i, 0] = x[_i]
            mw_array[_i, 1] = y[_i]
            mw_array = source.magnitude()

    mt_array = None
    if show_beachballs:
        mt_array = np.zeros((len(sources), 12))
        for _i, source in enumerate(sources):
            mt_array[_i, 0] = x[_i]
            mt_array[_i, 1] = y[_i]
            mt_array[_i, 3:9] = source.as_vector()

        mt_array[:,9] = 25.
        mt_array[:,10] = 0.
        mt_array[:,11] = 0.


    if mt_array is not None:
        mt_file = 'tmp_mt_'+filename+'.txt'
        np.savetxt(mt_file, mt_array)
    else:
        mt_file = "''"

    if mw_array is not None:
        mw_file = 'tmp_mw_'+filename+'.txt'
        np.savetxt(mw_file, mw_array)
    else:
        mw_file = "''"

    # call bash script
    if exists_gmt():
        subprocess.call("%s %s %s %s %s %s %f %f %f %f '%s' '%s' %s %s" %
           (fullpath('mtuq/graphics/_gmt/plot_mt_xy'),
            filename,
            filetype,
            mt_file,
            mw_file,
            '0',
            xmin, xmax,
            ymin, ymax,
            xlabel,
            ylabel,
            title,
            subtitle
            ),
            shell=True)
    else:
        gmt_not_found_warning(
            ascii_data)


def _get_labeltype(x,y,labeltype):
    if labeltype=='latlon':
       xlabel = None
       ylabel = None

    if labeltype=='offset' and ((x.max()-x.min()) >= 10000.):
       x /= 1000.
       y /= 1000.
       xlabel = 'E-W offset (km)'
       ylabel = 'N-S offset (km)'
    elif labeltype=='offset' and ((x.max()-x.min()) < 10000.):
       xlabel = 'E-W offset (m)'
       ylabel = 'N-S offset (m)'

    return xlabel,ylabel


def _get_limits(x,y):

    xmin, xmax = x.min(), x.max()
    ymin, ymax = y.min(), y.max()

    return xmin, xmax, ymin, ymax

