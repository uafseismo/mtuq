
#
# graphics/uq/force.py - uncertainty quantification of forces on the unit sphere
#

import numpy as np
import subprocess

from matplotlib import pyplot
from pandas import DataFrame
from xarray import DataArray

from mtuq.graphics.uq._gmt import exists_gmt, gmt_not_found_warning, \
    gmt_plot_misfit_force, gmt_plot_likelihood_force
from mtuq.grid_search import MTUQDataArray, MTUQDataFrame
from mtuq.util import fullpath
from mtuq.util.math import closed_interval, open_interval


def plot_misfit_force(filename, ds, misfit_callback=None, title='',
    colormap='viridis', colorbar_type=1, marker_type=1):
    """ Plots misfit values on the unit sphere


    .. rubric :: Input arguments

    ``filename`` (`str`):
    Name of output image file

    ``ds`` (`DataArray` or `DataFrame`):
    data structure containing forces and corresponding misfit values

    ``misfit_callback`` (func)
    User-supplied function applied to misfit values

    ``title`` (`str`):
    Optional figure title

    """
    _check(ds)
    ds = ds.copy()

    if issubclass(type(ds), DataArray):
        ds = ds.min(dim=('origin_idx', 'F0'))
        phi = ds.coords['phi']
        h = ds.coords['h']
        values = ds.values.transpose()

    elif issubclass(type(ds), DataFrame):
        ds = ds.reset_index()
        phi, h, values = _bin(ds, lambda ds: ds.min())

    if misfit_callback:
        values = misfit_callback(values)

    gmt_plot_misfit_force(filename, phi, h, values, title=title,
        colormap=colormap, colorbar_type=colorbar_type, marker_type=marker_type)


def plot_likelihood_force(filename, ds, sigma=None, title='',
    colormap='hot_r', colorbar_type=1, marker_type=2):

    """ Plots maximum likelihoods on the unit sphere


    .. rubric :: Input arguments

    ``filename`` (`str`):
    Name of output image file

    ``ds`` (`DataArray` or `DataFrame`):
    data structure containing forces and corresponding misfit values

    ``title`` (`str`):
    Optional figure title

    """
    _check(ds)
    ds = ds.copy()

    if issubclass(type(ds), DataArray):
        ds.values = np.exp(-ds.values/(2.*sigma**2))
        ds.values /= ds.values.sum()
        ds = ds.max(dim=('origin_idx', 'F0'))
        phi = ds.coords['phi']
        h = ds.coords['h']
        values = ds.values.transpose()

    elif issubclass(type(ds), DataFrame):
        ds[0] = np.exp(-ds[0]/(2.*sigma**2))
        ds = ds.reset_index()
        phi, h, values = _bin(ds, lambda ds: ds.max())

    values /= 4.*np.pi*values.sum()

    gmt_plot_likelihood_force(filename, phi, h, values, title=title,
        colormap=colormap, colorbar_type=colorbar_type, marker_type=marker_type)


def plot_marginal_force(filename, ds, sigma=None, title='',
    colormap='hot_r', colorbar_type=1, marker_type=2):
    """ Plots marginal likelihoods on the unit sphere


    .. rubric :: Input arguments

    ``filename`` (`str`):
    Name of output image file

    ``ds`` (`DataArray` or `DataFrame`):
    data structure containing forces and corresponding misfit values

    ``title`` (`str`):
    Optional figure title


    """
    _check(ds)
    ds = ds.copy()

    if issubclass(type(ds), DataArray):
        ds.values = np.exp(-ds.values/(2.*sigma**2))
        ds.values /= ds.values.sum()
        ds = ds.max(dim=('origin_idx', 'F0'))
        phi = ds.coords['phi']
        h = ds.coords['h']
        values = ds.values.transpose()

    elif issubclass(type(ds), DataFrame):
        ds = np.exp(-ds/(2.*sigma**2))
        #ds /= ds.sum()
        ds = ds.reset_index()
        phi, h, values = _bin(ds, lambda ds: ds.sum()/len(ds), normalize=True)

    values /= 4.*np.pi*values.sum()

    gmt_plot_likelihood_force(filename, phi, h, values, title=title,
        colormap=colormap, colorbar_type=colorbar_type, marker_type=marker_type)

def plot_force_amplitude(filename, ds, source_dict, title='',
    colormap='viridis', colorbar_type=1, marker_type=3):
    """ Plots force amplitude values on the unit sphere


    .. rubric :: Input arguments

    ``filename`` (`str`):
    Name of output image file

    ``ds`` (`DataArray`):
    Data structure containing forces and corresponding misfit values.
    Only supporting DataArray at the moment

    ``source_dict`` (`DataArray`):
    Solution's source dictionary used to define the optimal F0 value.

    ``title`` (`str`):
    Optional figure title

    """
    _check(ds)
    ds = ds.copy()

    if issubclass(type(ds), DataArray):
        phi = ds.coords['phi']
        h = ds.coords['h']

    elif issubclass(type(ds), DataFrame):
        raise NotImplementedError

    force_map = np.empty((ds.phi.size, ds.h.size))
    for iphi in range(ds.phi.size):
        for ih in range(ds.h.size):
            idx = ds[:,iphi,ih,:].values.argmin()
            force_map[iphi, ih] = ds.coords['F0'][idx]
    normalized_force_map=np.log(force_map / source_dict['F0']).transpose()

    #
    # Computes global minimum map coordinates
    #
    global_min_lon = source_dict['phi']
    if global_min_lon >= 360:
        global_min_lon -= 360
    global_min_lon += 90
    global_min_lat = (np.rad2deg(np.arccos(source_dict['h']))-90)*-1

    gmt_plot_misfit_force(filename, phi, h, normalized_force_map,
        title=title, colormap=colormap, colorbar_type=colorbar_type,
        marker_type=marker_type, global_min_lon=global_min_lon, global_min_lat=global_min_lat)

def _check(ds):
    """ Checks data structures
    """
    if type(ds) not in (DataArray, DataFrame, MTUQDataArray, MTUQDataFrame):
        raise TypeError("Unexpected grid format")


def _bin(df, handle, npts_phi=60, npts_h=30):
    """ Bins DataFrame into rectangular cells
    """
    # define centers of cells
    centers_phi = open_interval(0., 360., npts_phi)
    centers_h = open_interval(-1., +1., npts_h)

    # define corners of cells
    phi = closed_interval(0., 360, npts_phi+1)
    h = closed_interval(-1., +1., npts_h+1)

    binned = np.empty((npts_h, npts_phi))
    for _i in range(npts_h):
        for _j in range(npts_phi):
            # which grid points lie within cell (i,j)?
            subset = df.loc[
                df['phi'].between(phi[_j], phi[_j+1]) &
                df['h'].between(h[_i], h[_i+1])]

            if len(subset)==0:
                print("Encountered empty bin\n"
                      "phi: %f, %f\n"
                      "h: %f, %f\n" %
                      (phi[_j], phi[_j+1], h[_i], h[_i+1]) )

            binned[_i, _j] = handle(subset[0])

    return centers_phi, centers_h, binned
