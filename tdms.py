import numpy as np
import sys
try:
    from nptdms import TdmsFile
except ImportError:
    sys.stderr.write("nptdms module unavailable")

def read_tdms(fn):

    tdms_file = TdmsFile(fn)

    try:
        times = np.array(
            [[channel.data
              for channel in tdms_file.group_channels(group)
              if channel.data is not None]
             for group in tdms_file.groups()
             if group.lower() == "time"][0][0])
        dt = np.mean(np.diff(times))
    except IndexError:
        if not "Sampling Rate" in tdms_file.object().properties.keys():
            if not "Sampling Rate(AI)" in tdms_file.object().properties.keys():
                dt = 1.0
            else:
                sr = float(tdms_file.object().properties['Sampling Rate(AI)'])
                if sr > 0:
                    dt = 1e3/sr
                else:
                    dt = 1.0/25.0
        else:
            sr = float(tdms_file.object().properties['Sampling Rate'])
            if sr > 0:
                dt = 1e3/sr
            else:
                dt = 1.0/25.0

    yunits = tdms_file.object().properties['Units']
    try:
        meta = tdms_file.group_channels('Meta')
    except:
        meta = ''

    recording = {group: [
        channel.data
        for channel in tdms_file.group_channels(group)
        if channel.data is not None]
                 for group in tdms_file.groups()}
    recording["dt"] = dt
    recording["yunits"] = yunits
    recording["holding"] = meta

    return recording
