import numpy as np
from scipy.optimize import leastsq

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


