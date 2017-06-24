import numpy as np
from scipy.optimize import leastsq
import scipy

def leastsq_helper(p, y, lsfunc, x, *args):
    return y - lsfunc(p, x, *args)


def fboltz_up(p, x):
    """
    Boltzmann function (upwards from 0 to 1)

    Parameters
    ----------
    p : numpy.ndarray
        :math:`p_0`, location of half maximum; :math:`p_1`, slope
    x : numpy.ndarray
        Dependent variable

    Returns
    -------
    boltzmann : numpy.ndarray
        :math:`y = 1 - \\frac{1}{1 + e^\\left(\\left(x-p_0\\right)/p_1\\right)}`
    """
    return 1.0 - 1.0/(1.0+np.exp((x-p[0])/p[1]))


def fexp(p, x):
    """
    Exponential function

    Parameters
    ----------
    p : numpy.ndarray
        :math:`p_0`, amplitude; :math:`p_1`, :math:`\tau`; :math:`p_2`, offset
    x : numpy.ndarray
        Dependent variable

    Returns
    -------
    boltzmann : numpy.ndarray
        :math:`y = p_0 \left( e^{\frac{-x}{\tau}} \right) + p_2`
    """
    amp = p[0]
    tau = p[1]
    offset = p[2]
    return amp*(-np.exp(-x/tau)) + offset


def gv(i, v, erev):
    """
    Fit Boltzmann function to normalized conductance values

    Parameters
    ----------
    i : numpy.ndarray
        Peak current values
    v : numpy.ndarray
        Command voltages
    erev : float
        Reversal potential

    Returns
    -------
    g : numpy.ndarray
        Normalized conductances
    gfit : numpy.ndarray
        Half-maximal voltage and slope of best-fit Boltzmann function
    """
    g = i / (v-erev)
    g /= g.max()

    v50_init = 0.0
    slope_init = 1.0
    gfit = leastsq(
        leastsq_helper, (v50_init, slope_init), args=(g, fboltz_up, v))[0]

    return g, gfit


def tri_norm(x, *args):
    m1, m2, m3, s1, s2, s3, k1, k2, k3 = args
    ret = k1*scipy.stats.norm.pdf(x, loc=m1 ,scale=s1)
    ret += k2*scipy.stats.norm.pdf(x, loc=m2 ,scale=s2)
    ret += k3*scipy.stats.norm.pdf(x, loc=m3 ,scale=s3)
    return ret


def bi_norm(x, *args):
    m1, m2, s1, s2, k1, k2 = args
    ret = k1*scipy.stats.norm.pdf(x, loc=m1 ,scale=s1)
    ret += k2*scipy.stats.norm.pdf(x, loc=m2 ,scale=s2)
    return ret


def single_norm(x, *args):
    m1, s1, k1 = args
    return k1*scipy.stats.norm.pdf(x, loc=m1 ,scale=s1)
