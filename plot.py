import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

def plot_traces(mode, xlabel='Time (ms)'):
    gs = gridspec.GridSpec(4, 1)
    fig = plt.figure()
    axi = fig.add_subplot(gs[:3, 0])
    axv = fig.add_subplot(gs[3:, 0], sharex = axi)
    axi.spines['right'].set_color('none')
    axi.spines['top'].set_color('none')
    axi.spines['bottom'].set_color('none')
    plt.setp(axi.get_xticklabels(), visible=False)
    plt.setp(axi.get_xticklines(), visible=False)
    axv.spines['right'].set_color('none')
    axv.spines['top'].set_color('none')
    axi.spines['left'].set_smart_bounds(True)
    axv.spines['left'].set_smart_bounds(True)
    axv.spines['bottom'].set_smart_bounds(True)
    if mode is 'vclamp':
        axi.set_ylabel('Current (pA)')
        axv.set_ylabel('Voltage (mV)')
    else:
        axv.set_ylabel('Current (pA)')
        axi.set_ylabel('Voltage (mV)')

    axv.set_xlabel(xlabel)

    return fig, axi, axv


def plot_iv(xlabel='Voltage (mV)', ylabel='Current (pA)', zeroori=True):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    if zeroori:
        ax.spines['left'].set_position('zero')
        ax.spines['bottom'].set_position('zero')
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.spines['left'].set_smart_bounds(True)
    ax.spines['bottom'].set_smart_bounds(True)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    return fig, ax

def plot_fi(ylabel='Frequency (Hz)', zeroori=True):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    if zeroori:
        ax.spines['left'].set_position('zero')
        ax.spines['bottom'].set_position('zero')
    ax.spines['left'].set_smart_bounds(True)
    ax.spines['bottom'].set_smart_bounds(True)
    ax.set_xlabel('Current (pA)')
    ax.set_ylabel(ylabel)

    return fig, ax


def openfile_dialog():
    from PyQt4 import QtGui
    app = QtGui.QApplication([dir])
    fname = QtGui.QFileDialog.getOpenFileName(None, "Select a file...", '.', filter="All files (*)")
    return str(fname)
