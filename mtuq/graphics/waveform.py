
import numpy as np
import matplotlib.pyplot as pyplot
import warnings
from mtuq.graphics.header import attach_header, generate_header


#
# functions that generate entire figures
#

def plot_data_synthetics(filename, 
        data_bw,
        data_sw,
        synthetics_bw,
        synthetics_sw,
        stations,
        mt=None,
        total_misfit_bw=1., 
        total_misfit_sw=1., 
        normalize='maximum_amplitude',
        title=None, 
        header=None,
        station_labels=True, 
        trace_labels=True):

    """ Creates CAP-style data/synthetics figure
    """

    # how many stations have at least one trace?
    nstations = _count([data_bw, data_sw])

    # dimensions of subplot array
    nrow = nstations
    ncol = 5
    irow = 0

    # figure dimensions in inches
    height = 1.4*nrow
    width = 16.

    margin_bottom = 0.25
    margin_top = 0.25
    margin_left = 0.25
    margin_right = 0.25

    if station_labels:
        margin_left += 0.75

    height += margin_bottom
    height += margin_top
    width += margin_left
    width += margin_right

    if not title:
        title = filename.split('.')[0]

    if header:
        # attach optional CAP-style header
        header_height = 2.5
        height += header_height
        fig = pyplot.figure(figsize=(width, height))
        attach_header(title, header, mt, header_height)

    else:
        header_height = 0.
        fig = pyplot.figure(figsize=(width, height))


    # adjust subplot spacing
    pyplot.subplots_adjust(
        left=margin_left/width,
        right=1.-margin_right/width,
        bottom=margin_bottom/height,
        top=1.-(margin_top+header_height)/height,
        wspace=0.,
        hspace=0.,
        )

    # determine maximum amplitudes
    max_amplitude_bw = 0.
    if data_bw.max() > max_amplitude_bw:
        max_amplitude_bw = data_bw.max()
    if synthetics_bw.max() > max_amplitude_bw:
        max_amplitude_bw = synthetics_bw.max()

    max_amplitude_sw = 0.
    if data_sw.max() > max_amplitude_sw:
        max_amplitude_sw = data_sw.max()
    if synthetics_sw.max() > max_amplitude_sw:
        max_amplitude_sw = synthetics_sw.max()


    #
    # loop over stations
    #

    for _i in range(len(stations)):

        # skip empty stations
        if len(data_bw[_i])==len(data_sw[_i])==0:
            continue

        # add station labels
        if station_labels:
            meta = stations[_i]
            pyplot.subplot(nrow, ncol, ncol*irow+1)
            add_station_labels(meta)


        #
        # plot body wave traces
        #

        stream_dat = data_bw[_i]
        stream_syn = synthetics_bw[_i]

        for dat, syn in zip(stream_dat, stream_syn):
            component = dat.stats.channel[-1].upper()
            weight = getattr(dat, 'weight', 1.)

            # skip bad traces
            if component != syn.stats.channel[-1].upper():
                warnings.warn('Mismatched components, skipping...')
                continue
            elif weight==0.:
                continue

            # plot traces
            if component=='Z':
                pyplot.subplot(nrow, ncol, ncol*irow+1)
                plot(dat, syn)
            elif component=='R':
                pyplot.subplot(nrow, ncol, ncol*irow+2)
                plot(dat, syn)
            else:
                continue

            # normalize amplitudes
            if normalize=='trace_amplitude':
                max_trace = _max(dat, syn)
                ylim = [-2*max_trace, +2*max_trace]
                pyplot.ylim(*ylim)
            elif normalize=='maximum_amplitude':
                ylim = [-2*max_amplitude_bw, +2*max_amplitude_bw]
                pyplot.ylim(*ylim)

            if trace_labels:
                add_trace_labels(dat, syn, total_misfit_bw)


        #
        # plot surface wave traces
        #

        stream_dat = data_sw[_i]
        stream_syn = synthetics_sw[_i]

        for dat, syn in zip(stream_dat, stream_syn):
            component = dat.stats.channel[-1].upper()
            weight = getattr(dat, 'weight', 1.)

            # skip bad traces
            if component != syn.stats.channel[-1].upper():
                warnings.warn('Mismatched components, skipping...')
                continue
            elif weight==0.:
                continue

            # plot traces
            if component=='Z':
                pyplot.subplot(nrow, ncol, ncol*irow+3)
                plot(dat, syn)
            elif component=='R':
                pyplot.subplot(nrow, ncol, ncol*irow+4)
                plot(dat, syn)
            elif component=='T':
                pyplot.subplot(nrow, ncol, ncol*irow+5)
                plot(dat, syn)
            else:
                continue

            # amplitude normalization
            if normalize=='trace_amplitude':
                max_trace = _max(dat, syn)
                ylim = [-max_trace, +max_trace]
                pyplot.ylim(*ylim)
            elif normalize=='maximum_amplitude':
                ylim = [-max_amplitude_sw, +max_amplitude_sw]
                pyplot.ylim(*ylim)

            if trace_labels:
                add_trace_labels(dat, syn, total_misfit_sw)

        irow += 1

    pyplot.savefig(filename)



def plot_data_greens(filename, 
        data_bw,
        data_sw,
        greens_bw,
        greens_sw,
        process_bw,
        process_sw,
        misfit_bw,
        misfit_sw,
        stations,
        origin,
        mt,
        **kwargs):

    """ Creates CAP-style data/synthetics figure

    Similar to plot_data_synthetics, except provides different input argument
    syntax
    """
    event_name = filename.split('.')[0]

    greens_bw = greens_bw.select(origin)
    greens_sw = greens_sw.select(origin)
    synthetics_bw = greens_bw.get_synthetics(mt)
    synthetics_sw = greens_sw.get_synthetics(mt)
    total_misfit_bw = misfit_bw(data_bw, greens_bw, mt)
    total_misfit_sw = misfit_sw(data_sw, greens_sw, mt)

    header = generate_header(event_name,
        process_bw, process_sw, misfit_bw, misfit_bw,
        greens_bw[0].model, 'syngine', mt, origin.depth_in_m)

    plot_data_synthetics(filename,
            data_bw, data_sw, synthetics_bw, synthetics_sw, stations,
            total_misfit_bw=total_misfit_bw, total_misfit_sw=total_misfit_sw,
            mt=mt, header=header, **kwargs)



#
# functions that act on individual axes
#

def plot(dat, syn, label=None):
    """ Plots data and synthetics time series on current axes
    """
    t1,t2,nt,dt = _time_stats(dat)

    start = getattr(syn, 'start', 0)
    stop = getattr(syn, 'stop', len(syn.data))

    meta = dat.stats
    d = dat.data
    s = syn.data

    ax = pyplot.gca()

    t = np.linspace(0,t2-t1,nt,dt)
    ax.plot(t, d, 'k')
    ax.plot(t, s[start:stop], 'r')

    _hide_axes(ax)


def add_station_labels(meta):
    """ Displays station id, distance, and azimuth to the left of current axes
    """
    ax = pyplot.gca()

    # display station name
    label = '.'.join([meta.network, meta.station])
    pyplot.text(-0.25,0.45, label, fontsize=12, transform=ax.transAxes)

    # display distance
    distance = '%d km' % round(meta.preliminary_distance_in_m/1000.)
    pyplot.text(-0.25,0.30, distance, fontsize=12, transform=ax.transAxes)

    # display azimuth
    azimuth =  '%d%s' % (round(meta.preliminary_azimuth), u'\N{DEGREE SIGN}')
    pyplot.text(-0.25,0.15, azimuth, fontsize=12, transform=ax.transAxes)

    _hide_axes(ax)



def add_trace_labels(dat, syn, total_misfit=1.):
    """ Adds CAP-style annotations to current axes
    """
    ax = pyplot.gca()
    ymin = ax.get_ylim()[0]

    s = syn.data
    d = dat.data

    # display cross-correlation time shift
    time_shift = getattr(syn, 'time_shift', np.nan)
    pyplot.text(0.,(1/4.)*ymin, '%.2f' %time_shift, fontsize=12)

    # display maximum cross-correlation coefficient
    Ns = np.dot(s,s)**0.5
    Nd = np.dot(d,d)**0.5
    if Ns*Nd > 0.:
        max_cc = np.correlate(s, d, 'valid').max()
        max_cc /= (Ns*Nd)
        pyplot.text(0.,(2/4.)*ymin, '%.2f' %max_cc, fontsize=12)
    else:
        max_cc = np.nan
        pyplot.text(0.,(2/4.)*ymin, '%.2f' %max_cc, fontsize=12)

    # display percent of total misfit
    misfit = getattr(syn, 'misfit', np.nan)
    misfit /= total_misfit
    if misfit >= 0.1:
        pyplot.text(0.,(3/4.)*ymin, '%.1f' %(100.*misfit), fontsize=12)
    else:
        pyplot.text(0.,(3/4.)*ymin, '%.2f' %(100.*misfit), fontsize=12)



#
# utility functions
#

def _time_stats(trace):
    # returns time scheme
    return (
        float(trace.stats.starttime),
        float(trace.stats.endtime),
        trace.stats.npts,
        trace.stats.delta,
        )


def _count(datasets):
    # counts number of nonempty streams in dataset(s)
    count = 0
    for streams in zip(*datasets):
        for stream in streams:
            if len(stream) > 0:
                count += 1
                break
    return count


def _max(dat, syn):
    # maximum amplitude of two traces
    return max(
        abs(dat.max()),
        abs(syn.max()))


def _hide_axes(ax):
    # hides axes lines, ticks, and labels
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])


