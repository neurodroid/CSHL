import sys

import numpy as np

from mpl_toolkits.axes_grid.axislines import SubplotZero, Subplot
import matplotlib.gridspec as gridspec

from stfio import plot as stfio_plot


class ZeroAxis(Subplot):
    def __init__(self, *args, **kwargs):
        kwargs['frameon'] = False

        super(ZeroAxis, self).__init__(*args, **kwargs)

        args[0].add_axes(self)

        for direction in ["right", "top"]:
            self.axis[direction].set_visible(False)
           

def plot_iv(i, v, iunits, vunits, fig, subplot, sharex=None):

    ax = ZeroAxis(fig, subplot, sharex=sharex)
    ax.plot(v, i, '-o', color='k', ms=12, mec='none')
    ax.spines['left'].set_position('center')
    ax.spines['bottom'].set_position('center')

    ax.set_xlabel("$V$ ({0})".format(vunits))
    ax.set_ylabel("$I$ ({0})".format(iunits))

    return ax


def plot_gv(g, v, vunits, gfit, fig, subplot, sharex=None):

    ax = ZeroAxis(fig, subplot, sharex=sharex)
    ax.plot(v, g, 'o', color='k', ms=12, mec='none')
    xfit = np.arange(v.min(), v.max(), 0.1)
    ax.plot(xfit, fboltz_up(gfit, xfit), '-k')

    ax.spines['left'].set_position('center')
    ax.spines['bottom'].set_position('center')

    ax.set_xlabel("$V$ ({0})".format(vunits))
    ax.set_ylabel("$g/g_{max}$")

    return ax


def plot_traces(plotwindow=None, ichannel=0, vchannel=1):
    """
    Show traces in a figure

    Parameters
    ----------
    plotwindow : (float, float), optional
        Plot window (in ms from beginning of trace)
        None for whole trace. Default: None
    ichannel : int, optional
        current channel number. Default: 0
    vchannel : int, optional
        voltage channel number. Default: 1
    """

    import stf
    if not stf.check_doc():
        return None

    nchannels = stf.get_size_recording()
    if nchannels < 2:
        sys.stderr.write(
            "Function requires 2 channels (0: current; 1: voltage)\n")
        return

    dt = stf.get_sampling_interval()

    fig = stf.mpl_panel(figsize=(12, 8)).fig
    fig.clear()
    gs = gridspec.GridSpec(4, 1)
    ax_currents = stfio_plot.StandardAxis(
        fig, gs[:3, 0], hasx=False, hasy=False)
    ax_voltages = stfio_plot.StandardAxis(
        fig, gs[3:, 0], hasx=False, hasy=False, sharex=ax_currents)
    if plotwindow is not None:
        istart = int(plotwindow[0]/dt)
        istop = int(plotwindow[1]/dt)
    else:
        istart = 0
        istop = None

    for ntrace in range(stf.get_size_channel()):
        stf.set_trace(ntrace)
        stf.set_channel(ichannel)
        trace = stf.get_trace()[istart:istop]

        ax_currents.plot(np.arange(len(trace))*dt, trace)

        # Measure pulse amplitude
        stf.set_channel(vchannel)
        trace = stf.get_trace()[istart:istop]
        ax_voltages.plot(np.arange(len(trace))*dt, trace)

    # Reset active channel
    stf.set_channel(ichannel)

    stfio_plot.plot_scalebars(
        ax_currents, xunits=stf.get_xunits(), yunits=stf.get_yunits(channel=0))
    stfio_plot.plot_scalebars(
        ax_voltages, xunits=stf.get_xunits(), yunits=stf.get_yunits(channel=1))


def timeconstants(fitwindow, pulsewindow, ichannel=0, vchannel=1):
    """
    Compute and plot decay time constants

    Parameters
    ----------
    fitwindow : (float, float), optional
        Window for fitting time constant (time in ms from beginning of sweep)
        None for current cursor settings. Default: None
    pulsewindow : (float, float), optional
        Window for voltage pulse measurement (time in ms from beginning of sweep)
        None for current cursor settings. Default: None
    ichannel : int, optional
        current channel number. Default: 0
    vchannel : int, optional
        voltage channel number. Default: 1

    Returns
    -------
    v_commands : numpy.ndarray
        Command voltages
    taus : numpy.ndarray
        Time constants
    """

    import stf
    if not stf.check_doc():
        return None

    nchannels = stf.get_size_recording()
    if nchannels < 2:
        sys.stderr.write(
            "Function requires 2 channels (0: current; 1: voltage)\n")
        return

    dt = stf.get_sampling_interval()

    v_commands = []
    taus = []

    fig = stf.mpl_panel(figsize=(12, 8)).fig
    fig.clear()
    gs = gridspec.GridSpec(4, 8)
    ax_currents = stfio_plot.StandardAxis(
        fig, gs[:3, :4], hasx=False, hasy=False)
    ax_voltages = stfio_plot.StandardAxis(
        fig, gs[3:, :4], hasx=False, hasy=False, sharex=ax_currents)
    for ntrace in range(stf.get_size_channel()):
        stf.set_trace(ntrace)
        stf.set_channel(ichannel)
        trace = stf.get_trace()

        ax_currents.plot(np.arange(len(trace))*dt, trace)

        if fitwindow is not None:
            stf.fit.cursor_time = fitwindow
        res = stf.leastsq(0, False)
        taus.append(res['Tau_0'])

        # Measure pulse amplitude
        stf.set_channel(vchannel)
        trace = stf.get_trace()
        ax_voltages.plot(np.arange(len(trace))*dt, trace)

        stf.set_peak_direction("up")
        stf.set_peak_mean(-1)
        if pulsewindow is not None:
            stf.peak.cursor_time = pulsewindow
        stf.measure()
        v_commands.append(stf.peak.value)

    stfio_plot.plot_scalebars(
        ax_currents, xunits=stf.get_xunits(),
        yunits=stf.get_yunits(channel=ichannel))
    stfio_plot.plot_scalebars(
        ax_voltages, xunits=stf.get_xunits(),
        yunits=stf.get_yunits(channel=vchannel))

    v_commands = np.array(v_commands)
    taus = np.array(taus)

    ax_taus = plot_iv(
        taus, v_commands, "ms",
        stf.get_yunits(channel=vchannel), fig, 122)

    # Reset peak computation to single sampling point
    stf.set_peak_mean(1)

    # Reset active channel
    stf.set_channel(ichannel)

    # Compute conductances:
    stf.show_table_dictlist({
        "Voltage ({0})".format(
            stf.get_yunits(channel=vchannel)): v_commands.tolist(),
        "Taus (ms)": taus.tolist(),
    })

    return v_commands, taus

    
def iv(peakwindow=None, basewindow=None, pulsewindow=None,
       erev=None, peakmode="both", ichannel=0, vchannel=1,
       exclude=None):
    """
    Compute and plot an IV curve for currents

    Parameters
    ----------
    peakwindow : (float, float), optional
        Window for peak measurement (time in ms from beginning of sweep)
        None for current cursor settings. Default: None
    basewindow : (float, float), optional
        Window for baseline measurement (time in ms from beginning of sweep)
        None for current cursor settings. Default: None
    pulsewindow : (float, float), optional
        Window for voltage pulse measurement (time in ms from beginning of sweep)
        None for current cursor settings. Default: None
    erev : float, optional
        End of v clamp pulse in ms or None to determine automatically.
        Default: None
    peakmode : string, optional
        Peak direction - one of "up", "down", "both" or "mean". Default: "up"
    ichannel : int, optional
        current channel number. Default: 0
    vchannel : int, optional
        voltage channel number. Default: 1
    exclude : list of ints, optional
        List of trace indices to be excluded from the analysis. Default: None

    Returns
    -------
    v_commands : numpy.ndarray
        Command voltages
    ipeaks : numpy.ndarray
        Peak currents
    gpeaks : numpy.ndarray
        Peak normalized conductances
    g_fit : numpy.ndarray
        Half-maximal voltage and slope of best-fit Boltzmann function
    """

    import stf
    if not stf.check_doc():
        return None

    nchannels = stf.get_size_recording()
    if nchannels < 2:
        sys.stderr.write(
            "Function requires 2 channels (0: current; 1: voltage)\n")
        return

    dt = stf.get_sampling_interval()
    olddirection = stf.get_peak_direction()

    v_commands = []
    ipeaks = []
    if basewindow is not None:
        stf.base.cursor_time = basewindow

    fig = stf.mpl_panel(figsize=(12, 8)).fig
    fig.clear()
    gs = gridspec.GridSpec(4, 8)
    ax_currents = stfio_plot.StandardAxis(
        fig, gs[:3, :4], hasx=False, hasy=False)
    ax_voltages = stfio_plot.StandardAxis(
        fig, gs[3:, :4], hasx=False, hasy=False, sharex=ax_currents)
    for ntrace in range(stf.get_size_channel()):
        if exclude is not None:
            if ntrace in exclude:
                continue

        stf.set_trace(ntrace)
        stf.set_channel(ichannel)
        trace = stf.get_trace()

        ax_currents.plot(np.arange(len(trace))*dt, trace)

        # Measure only downward peaks (inward currents)
        if peakmode is "mean":
            stf.set_peak_direction("up")
            stf.set_peak_mean(-1)
        else:
            stf.set_peak_direction(peakmode)
            # Set peak computation to single sampling point
            stf.set_peak_mean(1)

        if peakwindow is not None:
            stf.peak.cursor_time = peakwindow
        stf.measure()
        if basewindow is not None:
            ipeaks.append(stf.peak.value-stf.base.value)
        else:
            ipeaks.append(stf.peak.value)

        # Measure pulse amplitude
        stf.set_channel(vchannel)
        trace = stf.get_trace()
        ax_voltages.plot(np.arange(len(trace))*dt, trace)

        stf.set_peak_direction("up")
        stf.set_peak_mean(-1)
        if pulsewindow is not None:
            stf.peak.cursor_time = pulsewindow
        stf.measure()
        v_commands.append(stf.peak.value)

    stfio_plot.plot_scalebars(
        ax_currents, xunits=stf.get_xunits(), yunits=stf.get_yunits(channel=0))
    stfio_plot.plot_scalebars(
        ax_voltages, xunits=stf.get_xunits(), yunits=stf.get_yunits(channel=1))

    v_commands = np.array(v_commands)
    ipeaks = np.array(ipeaks)

    if erev is None:
        # Find first zero crossing in ipeaks:
        for npulse in range(ipeaks.shape[0]-1):
            if np.sign(ipeaks[npulse]) != np.sign(ipeaks[npulse+1]):
                # linear interpolation
                m1 = (ipeaks[npulse+1]-ipeaks[npulse]) / (
                    v_commands[npulse+1]-v_commands[npulse])
                c1 = ipeaks[npulse] - m1*v_commands[npulse]
                erev = -c1/m1
                break
        if erev is None:
            sys.stderr.write(
                "Could not determine reversal potential. Aborting now\n")
            return None

    # Reset peak computation to single sampling point
    stf.set_peak_mean(1)
    stf.set_peak_direction(olddirection)

    # Reset active channel
    stf.set_channel(ichannel)

    # Compute conductances:
    gpeaks, g_fit = gv(ipeaks, v_commands, erev)

    ax_ipeaks = plot_iv(
        ipeaks, v_commands, stf.get_yunits(channel=ichannel),
        stf.get_yunits(channel=1), fig, 222)

    ax_ipeaks.set_title("Peak current")

    ax_gpeaks = plot_gv(
        gpeaks, v_commands, stf.get_yunits(channel=vchannel),
        g_fit, fig, 224)
    ax_gpeaks.set_title("Peak conductance")

    stf.show_table_dictlist({
        "Voltage ({0})".format(
            stf.get_yunits(channel=vchannel)): v_commands.tolist(),
        "Peak current ({0})".format(
            stf.get_yunits(channel=ichannel)): ipeaks.tolist(),
        "Peak conductance (g/g_max)": gpeaks.tolist(),
    })

    return v_commands, ipeaks, gpeaks, g_fit


def fi(pulsewindow=None, vthreshold=-45.0, vchannel=0, ichannel=1):
    """
    Compute and plot an f-I curve for current clamp recordings

    Parameters
    ----------
    pulsewindow : (float, float), optional
        Window for current pulse measurement (time in ms from beginning of sweep)
        None for current cursor settings. Default: None
    vthreshold : float, optional
        Voltage threshold of action potentials. Default: -45.0
    vchannel : int, optional
        voltage channel number. Default: 0
    ichannel : int, optional
        current channel number. Default: 1

    Returns
    -------
    i_commands : numpy.ndarray
        Command current pulses
    f : numpy.ndarray
        Action potential frequencies
    """
    f = None
    i_commands = None

    return i_commands, f
