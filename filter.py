from scipy.signal import bessel, lfilter


def lowpass(x, dt, f_c):
    """
    Lowpass filter

    Parameters
    ----------
    x : stfio_plot.Timeseries
        Input data
    dt : float
        Sampling interval in ms
    f_c : float
        Cutoff frequency in kHz (-3 dB)

    Returns
    -------
    x convolved with a Gaussian filter kernel.
    """
    fs = 1.0/dt
    cutoff = f_c
    B, A = bessel(1, cutoff / (fs / 2), btype='low') # 1st order Butterworth low-pass
    return lfilter(B, A, x, axis=0)


def highpass(x, dt, f_c):
    """
    Highpass filter

    Parameters
    ----------
    x : stfio_plot.Timeseries
        Input data
    dt : float
        Sampling interval in ms
    f_c : float
        Cutoff frequency in kHz (-3 dB)

    Returns
    -------
    x convolved with a Gaussian filter kernel.
    """
    fs = 1.0/dt
    cutoff = f_c
    B, A = bessel(1, cutoff / (fs / 2), btype='high') # 1st order Butterworth low-pass
    return lfilter(B, A, x, axis=0)


def bandpass(x, dt, f_lo, f_hi):
    """
    Highpass filter

    Parameters
    ----------
    x : stfio_plot.Timeseries
        Input data
    dt : float
        Sampling interval in ms
    f_c : float
        Cutoff frequency in kHz (-3 dB)

    Returns
    -------
    x convolved with a Gaussian filter kernel.
    """
    fs = 1.0/dt
    B, A = bessel(1, [f_lo / (fs / 2), f_hi / (fs / 2)], btype='bandpass') # 1st order Butterworth low-pass
    return lfilter(B, A, x, axis=0)
