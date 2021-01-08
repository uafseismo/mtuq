
#
# graphics/uq/lune.py - uncertainty quantification on the eigenvalue lune
#
# For details about the eigenvalue lune, see
# Tape2012 - A geometric setting for moment tensors
# (https://doi.org/10.1111/j.1365-246X.2012.05491.x)
#

import numpy as np

from matplotlib import pyplot
from pandas import DataFrame
from xarray import DataArray
from mtuq.graphics.uq._gmt import exists_gmt, gmt_not_found_warning, \
    gmt_plot_misfit_lune, gmt_plot_likelihood_lune, gmt_plot_misfit_mt_lune
from mtuq.grid_search import MTUQDataArray, MTUQDataFrame
from mtuq.util.math import lune_det, to_gamma, to_delta, to_v, to_w, semiregular_grid, to_mij, to_Mw


def plot_misfit_lune(filename, ds, misfit_callback=None, title='',
    colormap='viridis', colormap_reverse=False, colorbar_type=1, marker_type=1):

    """ Plots misfit values on eigenvalue lune (requires GMT)


    .. rubric :: Input arguments

    ``filename`` (`str`):
    Name of output image file

    ``ds`` (`DataArray` or `DataFrame`):
    Data structure containing moment tensors and corresponding misfit values

    ``misfit_callback`` (func)
    User-supplied function applied to misfit values

    ``title`` (`str`):
    Optional figure title

    """
    _check(ds)
    ds = ds.copy()

    if issubclass(type(ds), DataArray):
        ds = ds.min(dim=('origin_idx', 'rho', 'kappa', 'sigma', 'h'))
        gamma = to_gamma(ds.coords['v'])
        delta = to_delta(ds.coords['w'])
        values = ds.values.transpose()

    elif issubclass(type(ds), DataFrame):
        ds = ds.reset_index()
        gamma, delta, values = _bin(ds, lambda ds: ds.min())

    if misfit_callback:
        values = misfit_callback(values)

    gmt_plot_misfit_lune(filename, gamma, delta, values, title=title,
        colormap=colormap, colorbar_type=colorbar_type, marker_type=marker_type)


def plot_misfit_mt_lune(filename, ds, misfit_callback=None, title='',
    colormap='viridis', colormap_reverse=False, colorbar_type=1, marker_type=1):

    """ Plots misfit values and best fitting moment tensors on eigenvalue lune
    (requires GMT)


    .. rubric :: Input arguments

    ``filename`` (`str`):
    Name of output image file

    ``ds`` (`DataArray`):
    Data structure containing moment tensors and corresponding misfit values

    ``misfit_callback`` (func)
    User-supplied function applied to misfit values

    ``title`` (`str`):
    Optional figure title

    """
    _check(ds)
    ds = ds.copy()
    ds_for_plotting = ds.copy()

    if issubclass(type(ds), DataArray):
        ds = ds.min(dim=('origin_idx', 'rho', 'kappa', 'sigma', 'h'))
        gamma = to_gamma(ds.coords['v'])
        delta = to_delta(ds.coords['w'])
        values = ds.values.transpose()

    elif issubclass(type(ds), DataFrame):
        raise NotImplementedError


    _write_gmt_mt_params('tmp_FMT_coords.txt', ds_for_plotting, values)


    gmt_plot_misfit_mt_lune(filename, gamma, delta, values, title=title,
        colormap=colormap, colorbar_type=colorbar_type, marker_type=marker_type)

def plot_magnitude_lune(filename, ds, source_dict, misfit_callback=None, title='',
    colormap='viridis', colormap_reverse=False, colorbar_type=1, marker_type=3):

    """ Plots magnitude values on eigenvalue lune (requires GMT)


    .. rubric :: Input arguments

    ``filename`` (`str`):
    Name of output image file

    ``ds`` (`DataArray`):
    Data structure containing moment tensors and corresponding misfit values

    ``misfit_callback`` (func)
    User-supplied function applied to misfit values

    ``title`` (`str`):
    Optional figure title

    """
    _check(ds)
    ds = ds.copy()

    if issubclass(type(ds), DataArray):
        values = _extract_magnitude_map(ds)
        gamma = to_gamma(ds.coords['v'])
        delta = to_delta(ds.coords['w'])


    elif issubclass(type(ds), DataFrame):
        raise NotImplementedError

    if misfit_callback:
        values = misfit_callback(values)

    global_min_lon = to_gamma(source_dict['v'])
    global_min_lat = to_delta(source_dict['w'])

    gmt_plot_misfit_lune(filename, gamma, delta, values, title=title,
        colormap=colormap, colorbar_type=colorbar_type, marker_type=marker_type, global_min_lon=global_min_lon, global_min_lat=global_min_lat)


def plot_likelihood_lune(filename, ds, sigma=None, title='',
    colormap='hot_r', colorbar_type=1, marker_type=2):

    """ Plots maximum likelihoods on eigenvalue lune (requires GMT)


    .. rubric :: Input arguments

    ``filename`` (`str`):
    Name of output image file

    ``ds`` (`DataArray` or `DataFrame`):
    data structure containing moment tensors and corresponding misfit values

    ``sigma`` (`float`):
    Standard deviation applied to misfit values

    ``title`` (`str`):
    Optional figure title

    """
    assert sigma is not None

    _check(ds)
    ds = ds.copy()

    if issubclass(type(ds), DataArray):
        ds.values = np.exp(-ds.values/(2.*sigma**2))
        ds = ds.max(dim=('origin_idx', 'rho', 'kappa', 'sigma', 'h'))
        gamma = to_gamma(ds.coords['v'])
        delta = to_delta(ds.coords['w'])
        values = ds.values.transpose()

    elif issubclass(type(ds), DataFrame):
        ds = np.exp(-ds/(2.*sigma**2))
        ds = ds.reset_index()
        gamma, delta, values = _bin(ds, lambda ds: ds.max())

    #values /= lune_det(delta, gamma)

    area = (2./3.)*np.pi
    values /= area*values.sum()

    gmt_plot_likelihood_lune(filename, gamma, delta, values, title=title,
        colormap=colormap, colorbar_type=colorbar_type, marker_type=marker_type)


def plot_marginal_lune(filename, ds, sigma=None, title='',
    colormap='hot_r', colorbar_type=1, marker_type=2):

    """ Plots marginal likelihoods on eigenvalue lune (requires GMT)


    .. rubric :: Input arguments

    ``filename`` (`str`):
    Name of output image file

    ``ds`` (`DataArray` or `DataFrame`):
    data structure containing moment tensors and corresponding misfit values

    ``sigma`` (`float`):
    Standard deviation applied to misfit values

    ``title`` (`str`):
    Optional figure title

    """
    assert sigma is not None

    _check(ds)
    ds = ds.copy()

    if issubclass(type(ds), DataArray):
        ds.values = np.exp(-ds.values/(2.*sigma**2))
        ds = ds.sum(dim=('origin_idx', 'rho', 'kappa', 'sigma', 'h'))
        gamma = to_gamma(ds.coords['v'])
        delta = to_delta(ds.coords['w'])
        values = ds.values.transpose()

    elif issubclass(type(ds), DataFrame):
        ds = np.exp(-ds/(2.*sigma**2))
        ds = ds.reset_index()
        gamma, delta, values = _bin(ds, lambda ds: ds.sum()/len(ds), normalize=True)

    #values /= lune_det(delta, gamma)

    area = (2./3.)*np.pi
    values /= area*values.sum()

    gmt_plot_likelihood_lune(filename, gamma, delta, values, title=title,
        colormap=colormap, colorbar_type=colorbar_type, marker_type=marker_type)



def _check(ds):
    """ Checks data structures
    """
    if type(ds) not in (DataArray, DataFrame, MTUQDataArray, MTUQDataFrame):
        raise TypeError("Unexpected grid format")


def _bin(df, handle, npts_v=20, npts_w=40, tightness=0.6, normalize=False):
    """ Bins DataFrame into rectangular cells
    """
    # at which points will we plot values?
    centers_v, centers_w = semiregular_grid(
        npts_v, npts_w, tightness=tightness)

    # what cell edges correspond to the above centers?
    centers_gamma = to_gamma(centers_v)
    edges_gamma = np.array(centers_gamma[:-1] + centers_gamma[1:])/2.
    edges_v = to_v(edges_gamma)

    centers_delta = to_delta(centers_w)
    edges_delta = np.array(centers_delta[:-1] + centers_delta[1:])/2.
    edges_w = to_w(edges_delta)

    edges_v = np.pad(edges_v, 1)
    edges_v[0] = -1./3.
    edges_v[-1] = +1./3.

    edges_w = np.pad(edges_w, 1)
    edges_w[0] = -3.*np.pi/8.
    edges_w[-1] = +3.*np.pi/8


    # bin grid points into cells
    binned = np.empty((npts_w, npts_v))
    binned[:] = np.nan
    for _i in range(npts_w):
        for _j in range(npts_v):
            # which grid points lie within cell (i,j)?
            subset = df.loc[
                df['v'].between(edges_v[_j], edges_v[_j+1]) &
                df['w'].between(edges_w[_i], edges_w[_i+1])]

            if len(subset)==0:
                print("Encountered empty bin")

            binned[_i, _j] = handle(subset[0])

            if normalize:
              # normalize by area of cell
              binned[_i, _j] /= edges_v[_j+1] - edges_v[_j]
              binned[_i, _j] /= edges_w[_i+1] - edges_w[_i]

    return to_gamma(centers_v), to_delta(centers_w), binned


def _write_gmt_mt_params(filename, ds_for_plotting, misfit_values):
    """ Write full moment tensor parameters for GMT psmeca function in a
    temporary file. This temporary file is used in the plot_lune_mt GMT script
    before being removed from the working directory.
    """

    normalized_values = misfit_values.T.flatten() - np.min(misfit_values)
    normalized_values /= np.max(normalized_values)

    nv, nw = len(ds_for_plotting.coords['v']), len(ds_for_plotting.coords['w'])
    best_orientation=np.empty((nv*nw, 12))
    id = 0
    for iv in range(len(ds_for_plotting.coords['v'])):
        for iw in range(len(ds_for_plotting.coords['w'])):
            idx = np.unravel_index(np.argmin(ds_for_plotting[:,iv,iw,:,:,:,0].values, axis=None), np.shape(ds_for_plotting[:,iv,iw,:,:,:,0]))
            best_orientation[id, 0] = to_gamma(ds_for_plotting.coords['v'][iv])
            best_orientation[id, 1] = to_delta(ds_for_plotting.coords['w'][iw])
            best_orientation[id, 2] = normalized_values[id]
            rho, v, w, kappa, sigma, h = ds_for_plotting['rho'][idx[0]],\
                                        ds_for_plotting['v'][iv],\
                                        ds_for_plotting['w'][iw],\
                                        ds_for_plotting['kappa'][idx[1]],\
                                        ds_for_plotting['sigma'][idx[2]],\
                                        ds_for_plotting['h'][idx[3]]

            # Adding a random negative perturbation to the dip, to avoid GMT plotting bug.
            random_dip_perturbation = np.random.uniform(0.2,0.4)
            if sigma > (-90.0 + 0.4):
                sigma -= random_dip_perturbation
            mt = to_mij(rho, v, w, kappa, sigma, h)
            exponent = np.max([int('{:.2e}'.format(mt[i]).split('e+')[1]) for i in range(len(mt))])
            scaled_mt = mt/10**(exponent)

            best_orientation[id, 3:9] = scaled_mt
            best_orientation[id, 9] = exponent+7
            best_orientation[id, 10:12] = 0, 0
            id += 1
    np.savetxt(filename, best_orientation[:,:], fmt='%.4f')


def _extract_magnitude_map(ds):
    """ Extract and returns the optimal magnitude map to be plotted in
    "plot_magnitude_lune". For each location on the lune, the displayed
    magnitude corresponds to the best fitting moment tensor at this location.
    """
    M0 = to_Mw(ds.idxmin()['rho'].values)
    nv, nw = len(ds.coords['v']), len(ds.coords['w'])
    best_magnitude_map=np.empty((nv,nw))
    for iv in range(len(ds.coords['v'])):
        for iw in range(len(ds.coords['w'])):
            magnitude_idx = np.unravel_index(np.argmin(ds[:,iv,iw,:,:,:,0].values, axis=None), np.shape(ds[:,iv,iw,:,:,:,0]))[0]
            best_magnitude_map[iv, iw] = to_Mw(ds['rho'][magnitude_idx].values) - M0
    return(best_magnitude_map.T)
