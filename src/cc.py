from obspy import read, UTCDateTime, signal, Stream
from scipy import signal
import numpy as np
import os
import csv

startday = UTCDateTime(2004, 4, 1)
endday = UTCDateTime(2022, 3, 31)

stations = ["N.KISF"]
components = ["UB", "NB", "EB"]
fn_head = "./FNET/"
format = "sac"
dt = 0.05
dt_dec = 1.0
twin = 300
ntwin = int(twin/dt_dec)
dt_cc = 10
dir_out = "./cc/"

time_cc = np.arange(0, 86400, dt_cc)
b, a = signal.butter(2, [0.02, 0.05], "bandpass", fs=int(1/dt))
for station in stations:
    os.makedirs(dir_out + station, exist_ok=True)
    day = startday
    while (day <= endday):
        date = str(day.year)[-2:] + str(day.month).zfill(2) + str(day.day).zfill(2)
        day += 86400

        st = Stream()
        for component in components:
            fn = fn_head + date + "/" + station + "." + component
            if os.path.exists(fn):
                st += read(fn, format=format)
        if (len(st) < 3):
            continue
        if (abs(st[0].stats.npts - int(86400/dt)) > int(1/dt)) or \
            (abs(st[1].stats.npts - int(86400/dt)) > int(1/dt)) or \
            (abs(st[2].stats.npts - int(86400/dt)) > int(1/dt)):
            continue
        if (st[0].stats.npts != st[1].stats.npts) or \
            (st[0].stats.npts != st[2].stats.npts):
            continue
        st_hf = st.copy().filter(type="bandpass", freqmin=2, freqmax=8, \
                corners=2, zerophase=True)
        hf_sq = st_hf[0].data**2 + st_hf[1].data**2 + st_hf[2].data**2
        hf_sq_bp = signal.filtfilt(b, a, hf_sq)
        st_lf = st.copy().detrend("linear").filter(type="bandpass", freqmin=0.02, \
                freqmax=0.05, corners=2, zerophase=True).integrate().\
                detrend("linear")
        lf = st_lf[0].data

        hf_sq_bp = signal.decimate(hf_sq_bp, round(dt_dec/dt))
        lf = signal.decimate(lf, round(dt_dec/dt))
        cc = np.zeros(int(86400/dt_cc))
        for i in range(int(twin/2/dt_cc), int((86400-twin/2)/dt_cc)):
            i_wave = int(i*dt_cc/dt_dec)
            i_wave_l = i_wave - int(ntwin/2)
            i_wave_r = i_wave + int(ntwin/2)
            a1 = hf_sq_bp[i_wave_l:i_wave_r] - np.mean(hf_sq_bp[i_wave_l:i_wave_r])
            a2 = lf[i_wave_l:i_wave_r] - np.mean(lf[i_wave_l:i_wave_r])
            #cc[i] = np.dot(a1, a2) / np.linalg.norm(a2)
            cc[i] = np.dot(a1, a2) / \
                    (np.linalg.norm(a1, ord=2)*np.linalg.norm(a2, ord=2))

        np.savetxt(dir_out + station + "/" + date, np.stack([time_cc, cc], axis=1),\
                newline="\n", fmt=['%8.2f', '%10.5e'])
