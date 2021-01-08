
import numpy as np
from scipy.signal import fftconvolve


#
# numerical
#

def isclose(X, Y):
    EPSVAL = 1.e-6
    X = np.array(X)
    Y = np.array(Y)
    return bool(
        np.linalg.norm(X-Y) < EPSVAL)


def correlate(v1, v2):
    """ Fast cross-correlation function

    Correlates unpadded array v1 and padded array v2, producing result of
    shape ``len(v2) - len(v1)``
    """
    n1, n2 = len(v1), len(v2)

    if n1>2000 or n2-n1>200:
        # for long traces, frequency-domain implementation is usually faster
        return fftconvolve(v1, v2[::-1], 'valid')
    else:
        # for short traces, time-domain implementation is usually faster
        return np.correlate(v1, v2, 'valid')

def wrap_180(angle_in_deg):
    """ Wraps angle to (-180, 180)
    """
    angle_in_deg %= 360.
    idx = np.where(angle_in_deg > 180.)
    angle_in_deg[idx] -= 360.
    return angle_in_deg



#
# set theoretic
#


def list_intersect(a, b):
    """ Intersection of two lists
    """
    return list(set(a).intersection(set(b)))


def list_intersect_with_indices(a, b):
    intersection = list(set(a).intersection(set(b)))
    indices = [a.index(item) for item in intersection]
    return intersection, indices


def open_interval(x1, x2, N):
    """ Covers the open interval (x1, x2) with N regularly-spaced points
    """

    # NOTE: np.linspace(x1, x2, N)[1:-1] would be slightly simpler
    # but not as readily used by matplotlib.pyplot.pcolor

    return np.linspace(x1, x2, 2*N+1)[1:-1:2]


def closed_interval(x1, x2, N):
    """ Covers the closed interval [x1, x2] with N regularly-spaced points
    """
    return np.linspace(x1, x2, N)


def tight_interval(x1,x2,N,tightness=0.999):
    # tightness (float)
    # 0. reduces to ``open_intervel``, 1. reduces to ``closed_intervel``
    Lo = open_interval(x1,x2,N)
    Lc = closed_interval(x1,x2,N)
    return Lo*(1.-tightness) + Lc*tightness



#
# moment tensor and force
#

def lune_det(delta, gamma):
    """ Determinant of lune mapping as function of lune coordinates
    """
    delta, gamma = np.meshgrid(np.deg2rad(delta), np.deg2rad(gamma))
    beta = np.pi/2. - delta
    return 4./np.pi * np.sin(beta)**3 * np.cos(3.*gamma)


def to_mij(rho, v, w, kappa, sigma, h):
    """ Converts from lune parameters to moment tensor parameters
    (up-south-east convention)
    """
    kR3 = np.sqrt(3.)
    k2R6 = 2.*np.sqrt(6.)
    k2R3 = 2.*np.sqrt(3.)
    k4R6 = 4.*np.sqrt(6.)
    k8R6 = 8.*np.sqrt(6.)

    m0 = rho/np.sqrt(2.)

    delta, gamma = to_delta_gamma(v, w)
    beta = 90. - delta

    gamma = np.deg2rad(gamma)
    beta = np.deg2rad(90. - delta)
    kappa = np.deg2rad(kappa)
    sigma = np.deg2rad(sigma)
    theta = np.arccos(h)

    Cb  = np.cos(beta)
    Cg  = np.cos(gamma)
    Cs  = np.cos(sigma)
    Ct  = np.cos(theta)
    Ck  = np.cos(kappa)
    C2k = np.cos(2.0*kappa)
    C2s = np.cos(2.0*sigma)
    C2t = np.cos(2.0*theta)

    Sb  = np.sin(beta)
    Sg  = np.sin(gamma)
    Ss  = np.sin(sigma)
    St  = np.sin(theta)
    Sk  = np.sin(kappa)
    S2k = np.sin(2.0*kappa)
    S2s = np.sin(2.0*sigma)
    S2t = np.sin(2.0*theta)

    mt0 = m0 * (1./12.) * \
        (k4R6*Cb + Sb*(kR3*Sg*(-1. - 3.*C2t + 6.*C2s*St*St) + 12.*Cg*S2t*Ss))

    mt1 = m0* (1./24.) * \
        (k8R6*Cb + Sb*(-24.*Cg*(Cs*St*S2k + S2t*Sk*Sk*Ss) + kR3*Sg * \
        ((1. + 3.*C2k)*(1. - 3.*C2s) + 12.*C2t*Cs*Cs*Sk*Sk - 12.*Ct*S2k*S2s)))

    mt2 = m0* (1./6.) * \
        (k2R6*Cb + Sb*(kR3*Ct*Ct*Ck*Ck*(1. + 3.*C2s)*Sg - k2R3*Ck*Ck*Sg*St*St +
        kR3*(1. - 3.*C2s)*Sg*Sk*Sk + 6.*Cg*Cs*St*S2k +
        3.*Ct*(-4.*Cg*Ck*Ck*St*Ss + kR3*Sg*S2k*S2s)))

    mt3 = m0* (-1./2.)*Sb*(k2R3*Cs*Sg*St*(Ct*Cs*Sk - Ck*Ss) +
        2.*Cg*(Ct*Ck*Cs + C2t*Sk*Ss))

    mt4 = -m0* (1./2.)*Sb*(Ck*(kR3*Cs*Cs*Sg*S2t + 2.*Cg*C2t*Ss) +
        Sk*(-2.*Cg*Ct*Cs + kR3*Sg*St*S2s))

    mt5 = -m0* (1./8.)*Sb*(4.*Cg*(2.*C2k*Cs*St + S2t*S2k*Ss) +
        kR3*Sg*((1. - 2.*C2t*Cs*Cs - 3.*C2s)*S2k + 4.*Ct*C2k*S2s))

    if type(mt0) is np.ndarray:
        return np.column_stack([mt0, mt1, mt2, mt3, mt4, mt5])
    else:
        return np.array([mt0, mt1, mt2, mt3, mt4, mt5])


def to_xyz(F0, phi, h):
    """ Converts from spherical to Cartesian coordinates (east-north-up)
    """
    # spherical coordinates in "physics convention"
    r = F0
    phi = np.radians(phi)
    theta = np.arccos(h)

    x = F0*np.sin(theta)*np.cos(phi)
    y = F0*np.sin(theta)*np.sin(phi)
    z = F0*np.cos(theta)

    if type(F0) is np.ndarray:
        return np.column_stack([x, y, z])
    else:
        return np.array([x, y, z])


def to_rtp(F0, phi, h):
    """ Converts from spherical to Cartesian coordinates (up-south-east)
    """
    # spherical coordinates in "physics convention"
    r = F0
    phi = np.radians(phi)
    theta = np.arccos(h)

    x = F0*np.sin(theta)*np.cos(phi)
    y = F0*np.sin(theta)*np.sin(phi)
    z = F0*np.cos(theta)

    if type(F0) is np.ndarray:
        return np.column_stack([z, -y, x,])
    else:
        return np.array([z, -y, x])


def to_delta_gamma(v, w):
    """ Converts from Tape2015 parameters to lune coordinates
    """
    return to_delta(w), to_gamma(v)


def to_gamma(v):
    """ Converts from Tape2015 parameter v to lune longitude
    """
    gamma = (1./3.)*np.arcsin(3.*v)
    return np.rad2deg(gamma)


def to_delta(w):
    """ Converts from Tape2015 parameter w to lune latitude
    """
    beta0 = np.linspace(0, np.pi, 100)
    u0 = 0.75*beta0 - 0.5*np.sin(2.*beta0) + 0.0625*np.sin(4.*beta0)
    beta = np.interp(3.*np.pi/8. - w, u0, beta0)
    delta = np.rad2deg(np.pi/2. - beta)
    return delta


def to_v_w(delta, gamma):
    """ Converts from lune coordinates to Tape2015 parameters
    """
    return to_v(gamma), to_w(delta)


def to_v(gamma):
    """ Converts from lune longitude to Tape2015 parameter v
    """
    v = (1./3.)*np.sin(3.*np.deg2rad(gamma))
    return v


def to_w(delta):
    """ Converts from lune latitude to Tape2015 parameter w
    """
    beta = np.deg2rad(90. - delta)
    u = (0.75*beta - 0.5*np.sin(2.*beta) + 0.0625*np.sin(4.*beta))
    w = 3.*np.pi/8. - u
    return w


def to_M0(Mw):
    """ Converts from moment magnitude to scalar moment
    """
    return 10.**(1.5*float(Mw) + 9.1)


def to_rho(Mw):
    """ Converts from moment magnitude to Tape2012 magnitude parameter
    """
    return to_M0(Mw)*np.sqrt(2.)

def to_Mw(rho):
    """ Converts from Tape2012 magnitude parameter to moment magnitude
    """
    return ((np.log10(rho/np.sqrt(2))-9.1)/1.5)


def semiregular_grid(npts_v, npts_w, tightness=0.5):
    """ Semiregular moment tensor grid

    For tightness~0, grid will be regular in Tape2012 parameters delta, gamma.
    For tightness~1, grid will be regular in Tape2015 parameters v, w.
    For intermediate values, the grid will be "semiregular" in the sense of
    a linear interpolation between the above cases.
    """
    assert 0. <= tightness < 1.,\
        Exception("Allowable range: 0. <= tightness < 1.")

    v1 = open_interval(-1./3., 1./3., npts_v)
    w1 = open_interval(-3./8.*np.pi, 3./8.*np.pi, npts_w)

    gamma1 = to_gamma(v1)
    delta1 = to_delta(w1)

    gamma2 = np.linspace(-29.5, 29.5, npts_v)
    delta2 = np.linspace(-87.5, 87.5, npts_w)

    delta = delta1*(1.-tightness) + delta2*tightness
    gamma = gamma1*(1.-tightness) + gamma2*tightness

    return to_v(gamma), to_w(delta)
