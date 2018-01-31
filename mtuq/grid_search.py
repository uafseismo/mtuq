

import numpy as np


def grid_search(data, greens, misfit, grid):
    """ Grid search over moment tensor parameters only
    """
    for _i in grid.size:
        # generate moment tensor 
        mt = grid.index2mt(_i)

        # generate_synthetics
        categories = data.keys()
        synthetics = {}
        for key in categories:
            synthetics[key] = greens[key].combine(mt)

        sum_misfit = 0.
        for key in categories:
            sum_misfit += misfit(data[key], synthetics[key])


def grid_search_mpi(data, greens, misfit, grid):
    raise NotImplementedError




### grid iterators

class BaseGridRegular(object):
    """ Base class
    """
    def __init__(self, bounds, type):
        if not hasattr(maps, type):
            raise ValueError

        # parameter names
        keys = sorted(bounds.keys())
        ndim = len(keys)

        # grid shape
        shape = []
        for key in range(ndim):
            shape += bounds[_i,2]

        self.keys = keys
        self.bounds = bounds

        self.ndim = ndim
        self.shape = shape
        self.size = np.product(shape)

        # define coordinate vectors
        self.define_coords()

        # define parameter conversion 
        self.map = getattr(maps, type)


    def define_coords(self):
        # define coordinate vectors for each parameter
        self._coords = {}
        for key in keys:
            self._vectors[key] = [np.linspace(bounds[key][0], bounds[key][1], nsamples+2)[1:-1]]


    def index2mt(self, index):
        t = index2tuple(index)
        m = {}
        for key in self.keys():
            m[key] = self._coords[key][t[key]]

        # convert to mt vector (gcmt convention)
        return self.map(m)


    def index2tuple(self, index):
       t = {}
       for ikey, key in enumerate(self.keys):
           t[key] = index % self.shape[ikey]
       return t


class BaseGridRandom(BaseGridRegular):
    """ Base class
    """
    def define_coords(self):
        self._coords = {}
        for key in keys:
            self._vectors[key] = [np.random.randvec(*bounds[key])]



class MTGridRandom(BaseGridRandom):
    """ Full moment tensor grid
    """
    def __init__(self):
        super(MTGridRandom, self).__init__()


class MTGridRegular(BaseGridRegular):
    """ Full moment tensor grid
    """
    def __init__(self):
        super(MTGridRregular, self).__init__()


class DTGridRandom:
    """ Double-couple moment tensor grid
    """
    def __init__(self):
        raise NotImplementedError


class DTGridRegular:
    """ Double-couple moment tensor grid
    """
    def __init__(self):
        raise NotImplementedError


