import sys

import numpy as np
from scipy.optimize import leastsq

from mpl_toolkits.axes_grid.axislines import SubplotZero
import matplotlib.gridspec as gridspec

from stfio import plot as stfio_plot


class ZeroAxis(SubplotZero):
    def __init__(self, *args, **kwargs):
        kwargs['frameon'] = False

        super(ZeroAxis, self).__init__(*args, **kwargs)

        args[0].add_axes(self)

        for direction in ["xzero", "yzero"]:
            self.axis[direction].set_visible(True)

        for direction in ["left", "right", "bottom", "top"]:
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


def leastsq_helper(p, y, lsfunc, x, *args):
    return y - lsfunc(p, x, *args)


def fboltz_up(p, x):
    """
    Boltzmann function (upwards from 0 to 1)

    Parameters
    ==========
    p -- numpy.ndarray
        p[0]: location of half maximum
        p[1]: slope
    x -- numpy.ndarray
        Dependent variable

    Returns
    =======
    $y = 1 - \frac{1}{1 + e^{x-p_0}/p_1}
    """
    return 1.0 - 1.0/(1.0+np.exp((x-p[0])/p[1]))


def iv(window, erev, peakmode="up"):
    """
    Compute and plot an IV curve for currents

    Parameters
    ==========
    window -- (float, float)
        Window for peak measurement, starting from beginning of pulse
    erev -- float
        End of v clamp pulse in ms or None to determine automatically
    """
    import stf
    nchannels = stf.get_size_recording()
    if nchannels < 2:
        sys.stderr.write(
            "Function requires 2 channels (0: current; 1: voltage)\n")
        return

    dt = stf.get_sampling_interval()

    pulse = stf.get_trace(trace=stf.get_size_channel()-1, channel=1)

    start_pulse = np.argmax(np.diff(pulse))*dt
    end_pulse = np.argmin(np.diff(pulse))*dt

    v_commands = []
    ipeaks = []
    stf.base.cursor_time = (start_pulse-20.0, start_pulse-10.0)

    fig = stf.mpl_panel(figsize=(12, 8)).fig
    fig.clear()
    gs = gridspec.GridSpec(4, 8)
    ax_currents = stfio_plot.StandardAxis(
        fig, gs[:3, :4], hasx=False, hasy=False)
    ax_voltages = stfio_plot.StandardAxis(
        fig, gs[3:, :4], hasx=False, hasy=False, sharex=ax_currents)
    for ntrace in range(stf.get_size_channel()):
        stf.set_trace(ntrace)
        stf.set_channel(0)
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

        stf.peak.cursor_time = (
            start_pulse+window[0], start_pulse+window[1])
        stf.measure()
        ipeaks.append(stf.peak.value)

        # Measure pulse amplitude
        stf.set_channel(1)
        trace = stf.get_trace()
        ax_voltages.plot(np.arange(len(trace))*dt, trace)

        stf.set_peak_direction("up")
        stf.set_peak_mean(-1)
        stf.peak.cursor_time = (end_pulse-20.0, end_pulse-10.0)
        stf.measure()
        v_commands.append(stf.peak.value)

    stfio_plot.plot_scalebars(
        ax_currents, xunits=stf.get_xunits(), yunits=stf.get_yunits(channel=0))
    stfio_plot.plot_scalebars(
        ax_voltages, xunits=stf.get_xunits(), yunits=stf.get_yunits(channel=1))

    v_commands = np.array(v_commands)
    ipeaks = np.array(ipeaks)

    # Reset peak computation to single sampling point
    stf.set_peak_mean(1)

    # Reset active channel
    stf.set_channel(0)

    # Compute conductances:
    gpeaks, g_fit = gv(ipeaks, v_commands, erev)

    ax_ipeaks = plot_iv(
        ipeaks, v_commands, stf.get_yunits(channel=0),
        stf.get_yunits(channel=1), fig, 222)

    ax_ipeaks.set_title("Peak current")

    ax_gpeaks = plot_gv(
        gpeaks, v_commands, stf.get_yunits(channel=1),
        g_fit, fig, 224)
    ax_gpeaks.set_title("Peak conductance")

    stf.show_table_dictlist({
        "Voltage ({0})".format(
            stf.get_yunits(channel=1)): v_commands.tolist(),
        "Peak current ({0})".format(
            stf.get_yunits(channel=0)): ipeaks.tolist(),
        "Peak conductance (g/g_max)": gpeaks.tolist(),
    })

    return v_commands, ipeaks, gpeaks, g_fit


def gv(i, v, erev):
    g = i / (v-erev)
    g /= g.max()

    v50_init = 0.0
    slope_init = 1.0
    gfit = leastsq(
        leastsq_helper, (v50_init, slope_init), args=(g, fboltz_up, v))[0]

    return g, gfit
