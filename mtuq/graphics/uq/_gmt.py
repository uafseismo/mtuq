
import numpy as np
import shutil
import subprocess

from mtuq.graphics._gmt import exists_gmt, gmt_not_found_warning, gmt_version,\
    gmt_formats
from mtuq.util import fullpath, warn
from mtuq.util.math import wrap_180
from os.path import basename, exists, splitext



def gmt_plot_misfit_lune(filename, lon, lat, values, colormap='panoply', **kwargs):

    if _nothing_to_plot(values):
        return

    lon, lat =  _parse_lonlat(lon,lat)
    values, minval, maxval, exp = _parse_values(values)

    _call(fullpath('mtuq/graphics/uq/_gmt/plot_lune'),
        filename,
        lon, lat, values,
        z_min=minval,
        z_max=maxval,
        z_exp=exp,
        cpt_name=colormap,
        cpt_step=(maxval-minval)/20.,
        **kwargs)


def gmt_plot_misfit_mt_lune(filename, lon, lat, values, colormap='panoply', **kwargs):

    if _nothing_to_plot(values):
        return

    lon, lat =  _parse_lonlat(lon,lat)
    values, minval, maxval, exp = _parse_values(values)

    _call(fullpath('mtuq/graphics/uq/_gmt/plot_lune_mt'),
        filename,
        lon, lat, values,
        z_min=minval,
        z_max=maxval,
        z_exp=exp,
        cpt_name=colormap,
        cpt_step=(maxval-minval)/20.,
        **kwargs)

def gmt_plot_likelihood_lune(filename, lon, lat, values, colormap='hot', **kwargs):

    if _nothing_to_plot(values):
        return

    lon, lat =  _parse_lonlat(lon,lat)
    values, minval, maxval, exp = _parse_values(values)

    _call(fullpath('mtuq/graphics/uq/_gmt/plot_lune'),
        filename,
        lon, lat, values,
        z_min=minval,
        z_max=maxval,
        z_exp=exp,
        cpt_name=colormap,
        cpt_step=(maxval-minval)/20.,
        **kwargs)


def gmt_plot_misfit_force(filename, phi, h, values, colormap='panoply', **kwargs):

    if _nothing_to_plot(values):
        return

    lat = np.degrees(np.pi/2 - np.arccos(h))
    lon = wrap_180(phi + 90.)

    lon, lat =  _parse_lonlat(lon,lat)
    values, minval, maxval, exp = _parse_values(values)

    _call(fullpath('mtuq/graphics/uq/_gmt/plot_force'),
        filename,
        lon, lat, values,
        z_min=minval,
        z_max=maxval,
        z_exp=exp,
        cpt_name=colormap,
        cpt_step=(maxval-minval)/20.,
        **kwargs)


def gmt_plot_likelihood_force(filename, phi, h, values, colormap='hot', **kwargs):

    if _nothing_to_plot(values):
        return

    lat = np.degrees(np.pi/2 - np.arccos(h))
    lon = wrap_180(phi + 90.)

    lon, lat =  _parse_lonlat(lon,lat)
    values, minval, maxval, exp = _parse_values(values)

    _call(fullpath('mtuq/graphics/uq/_gmt/plot_force'),
        filename,
        lon, lat, values,
        z_min=minval,
        z_max=maxval,
        z_exp=exp,
        cpt_name=colormap,
        cpt_step=(maxval-minval)/20.,
        **kwargs)


def _call(shell_script, filename, lon, lat, values,
    z_min=None, z_max=None, z_exp=0,
    cpt_name='panoply', cpt_step=None, cpt_reverse=False,
    colorbar_type=0, marker_type=0, title='', global_min_lon=None, global_min_lat=None):

    print('  calling GMT script: %s' % basename(shell_script))

    filetype = _parse_filetype(filename)
    title, subtitle = _parse_title(title)

    # write lon,lat,val ASCII table
    ascii_data = 'tmp_'+filename+'.txt'
    _savetxt(ascii_data, lon, lat, values)

    cpt_local = fullpath('mtuq/graphics/_gmt/cpt', cpt_name+'.cpt')

    if exists(cpt_local):
       cpt_name = cpt_local

    # call bash script
    if exists_gmt():
        subprocess.call("%s %s %s %s %e %e %d %e %s %d %d %d %s %s %s %s" %
           (shell_script,
            ascii_data,
            filename,
            filetype,
            z_min,
            z_max,
            z_exp,
            cpt_step,
            cpt_name,
            int(bool(cpt_reverse)),
            int(colorbar_type),
            int(marker_type),
            title,
            subtitle,
            global_min_lon,
            global_min_lat
            ),
            shell=True)
    else:
        gmt_not_found_warning(
            ascii_data)



#
# utility functions
#

def _nothing_to_plot(values):
    mask = np.isnan(values)
    if np.all(mask):
        warn(
            "Nothing to plot: all values are NaN",
            Warning)
        return True

    masked = np.ma.array(values, mask=mask)
    minval = masked.min()
    maxval = masked.max()

    if minval==maxval:
        warn(
            "Nothing to plot: all values are identical",
            Warning)
        return True


def _parse_lonlat(lon, lat):

    lon, lat = np.meshgrid(lon, lat)
    lat = lat.flatten()
    lon = lon.flatten()
    return lon, lat


def _parse_values(values):

    values = values.flatten()
    masked = np.ma.array(values, mask=np.isnan(values))

    minval = masked.min()
    maxval = masked.max()
    exp = np.floor(np.log10(np.max(np.abs(masked))))

    if -1 <= exp <= 2:
        return masked, minval, maxval, 0

    else:
        masked /= 10**exp
        minval /= 10**exp
        maxval /= 10**exp
        return masked, minval, maxval, exp



def _parse_title(title):

    try:
        parts=title.split('\n')
    except:
        parts=[]

    if len(parts) >= 2:
        title = "'%s'" % parts[0]
        subtitle = "'%s'" % parts[1]
    elif len(parts) == 1:
        title = "'%s'" % parts[0]
        subtitle = "''"
    else:
        title = "''"
        subtitle = "''"

    return title, subtitle


def _parse_filetype(filename):

    parts = splitext(filename)
    name, ext = parts[0], parts[1].lstrip('.')

    if ext.upper() in ['PS']:
        return 'EPS'

    elif ext.upper() in ['JPG', 'JPEG']:
        return 'JPEG'

    elif ext.upper() in ['TIF', 'TIFF']:
        return 'TIFF'

    elif ext.upper() in gmt_formats:
        return ext.upper()

    else:
        warn('Unrecognized extension: defaulting to PNG')
        return 'PNG'



def _savetxt(filename, gamma, delta, values):
    # FIXME: can GMT accept virtual files?
    np.savetxt(filename, np.column_stack([gamma, delta, values]))

#
