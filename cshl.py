import numpy as np
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


def whole_cell_iv(start_pulse=None, end_pulse=None):
    """
    Compute and plot a whole-cell IV curve for an open file with Stimfit.

    Parameters
    ==========
    start_pulse -- float
        Start of v clamp pulse in ms or None to determine automatically
    end_pulse -- float
        End of v clamp pulse in ms or None to determine automatically
    """
    import stf
    dt = stf.get_sampling_interval()

    pulse = stf.get_trace(trace=stf.get_size_channel()-1, channel=1)
    if start_pulse is None:
        start_pulse = np.argmax(np.diff(pulse))*dt
    if end_pulse is None:
        end_pulse = np.argmin(np.diff(pulse))*dt
        
    v_commands = []
    inward_peaks, outward_peaks = [], []
    stf.base.cursor_time = (start_pulse-20.0, start_pulse-10.0)

    fig = stf.mpl_panel(figsize=(12, 8)).fig
    fig.clear()
    gs = gridspec.GridSpec(4, 5)
    ax_currents = stfio_plot.StandardAxis(fig, gs[:3, :3], hasx=False, hasy=False)
    ax_voltages = stfio_plot.StandardAxis(
        fig, gs[3:, :3], hasx=False, hasy=False, sharex=ax_currents)
    for ntrace in range(stf.get_size_channel()):
        stf.set_trace(ntrace)
        stf.set_channel(0)
        trace = stf.get_trace()

        ax_currents.plot(np.arange(len(trace))*dt, trace)

        # Measure only downward peaks (inward currents)
        stf.set_peak_direction("down")
        # Set peak computation to single sampling point
        stf.set_peak_mean(1)
        stf.peak.cursor_time = (start_pulse+0.4, start_pulse+20.0)
        stf.measure()
        inward_peaks.append(stf.peak.value)

        stf.set_peak_direction("up")
        # Set peak computation to mean of peak window
        stf.set_peak_mean(-1)
        stf.peak.cursor_time = (end_pulse-20.0, end_pulse-10.0)
        stf.measure()
        outward_peaks.append(stf.peak.value)

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

    # Reset peak computation to single sampling point
    stf.set_peak_mean(1)

    ax_inward = plot_iv(
        inward_peaks, v_commands, stf.get_yunits(channel=0),
        stf.get_yunits(channel=1), fig, 222)
    ax_outward = plot_iv(
        outward_peaks, v_commands, stf.get_yunits(channel=0),
        stf.get_yunits(channel=1), fig, 224, sharex=ax_inward)

    ax_inward.set_title("Peak inward current")
    ax_outward.set_title("Steady-state outward current")    
