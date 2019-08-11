
import os
import numpy as np

from mtuq import read, open_db, download_greens_tensors
from mtuq.graphics import plot_data_greens, plot_beachball
from mtuq.grid import DoubleCoupleGridRegular
from mtuq.grid_search import grid_search
from mtuq.misfit import Misfit
from mtuq.process_data import ProcessData
from mtuq.util import fullpath
from mtuq.util.cap import Trapezoid



if __name__=='__main__':
    #
    # Grid search integration test
    #
    # This script is similar to examples/SerialGridSearch.DoubleCouple.py,
    # except here we use a coarser grid, and at the end we assert that the test
    # result equals the expected result
    #
    # The compare against CAP/FK:
    #
    # cap.pl -H0.02 -P1/15/60 -p1 -S2/10/0 -T15/150 -D1/1/0.5 -C0.1/0.333/0.025/0.0625 -Y1 -Zweight_test.dat -Mscak_34 -m4.5 -I1/1/10/10/10 -R0/0/0/0/0/360/0/90/-180/180 20090407201255351
    #
    # Note however that CAP uses a different method for defining regular grids
    #


    # by default, the script runs with figure generation and error checking
    # turned on
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--no_checks', action='store_true')
    parser.add_argument('--no_figures', action='store_true')
    args = parser.parse_args()
    run_checks = (not args.no_checks)
    run_figures = (not args.no_figures)


    path_greens=  fullpath('data/tests/benchmark_cap/greens/scak')
    path_data=    fullpath('data/examples/20090407201255351/*.[zrt]')
    path_weights= fullpath('data/examples/20090407201255351/weights_all.dat')
    event_name=   '20090407201255351'
    model=        'scak'


    process_bw = ProcessData(
        filter_type='Bandpass',
        freq_min= 0.1,
        freq_max= 0.333,
        pick_type='from_fk_metadata',
        fk_database=path_greens,
        window_type='cap_bw',
        window_length=15.,
        padding_length=0.,
        weight_type='cap_bw',
        cap_weight_file=path_weights,
        )

    process_sw = ProcessData(
        filter_type='Bandpass',
        freq_min=0.025,
        freq_max=0.0625,
        pick_type='from_fk_metadata',
        fk_database=path_greens,
        window_type='cap_sw',
        window_length=150.,
        padding_length=0.,
        weight_type='cap_sw',
        cap_weight_file=path_weights,
        )


    misfit_bw = Misfit(
        norm='L2',
        time_shift_max=0.,
        time_shift_groups=['Z','R'],
        )

    misfit_sw = Misfit(
        norm='L2',
        time_shift_max=0.,
        time_shift_groups=['Z','R','T'],
        )


    #
    # Following obspy, we use the variable name "source" for the mechanism of
    # an event and "origin" for the location of an event
    #

    sources = DoubleCoupleGridRegular(
        npts_per_axis=1,
        magnitude=4.5)

    wavelet = Trapezoid(
        magnitude=4.5)


    #
    # The main I/O work starts now
    #

    print 'Reading data...\n'
    data = read(path_data, format='sac',
        event_id=event_name,
        tags=['units:cm', 'type:velocity']) 

    data.sort_by_distance()

    stations = data.get_stations()
    origin = data.get_origins()[0]


    print 'Processing data...\n'
    data_bw = data.map(process_bw)
    data_sw = data.map(process_sw)

    print 'Reading Greens functions...\n'
    db = open_db(path_greens, format='FK', model=model)
    greens = db.get_greens_tensors(stations, origin)

    print 'Processing Greens functions...\n'
    greens.convolve(wavelet)
    greens_bw = greens.map(process_bw)
    greens_sw = greens.map(process_sw)


    #
    # The main computational work starts nows
    #

    print 'Evaluating body wave misfit...\n'

    results_bw_0 = misfit_bw(
        data_bw, greens_bw, sources, optimization_level=0)

    results_bw_1 = misfit_bw(
        data_bw, greens_bw, sources, optimization_level=1)

    results_bw_2 = misfit_bw(
        data_bw, greens_bw, sources, optimization_level=2)

    print '\n\n\n'
    print results_bw_0.min()
    print results_bw_1.min()
    print results_bw_2.min()

    print 'Evaluating surface wave misfit...\n'

    results_sw_0 = misfit_sw(
        data_sw, greens_sw, sources, optimization_level=0)

    results_sw_1 = misfit_sw(
        data_sw, greens_sw, sources, optimization_level=1)

    results_sw_2 = misfit_sw(
        data_sw, greens_sw, sources, optimization_level=2)


    print '\n\n\n'
    print results_sw_0.min()
    print results_sw_1.min()
    print results_sw_2.min()
