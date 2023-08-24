from obspy import read, UTCDateTime, Stream
import numpy as np
import csv
import os

stations = ["N.KISF"]
components = ["UB", "EB", "NB"]
years = list(range(2004, 2022+1))
format = "sac"
dt = 0.05
fn_head = "./FNET/"
dir_out = "./rms/"

rms_hf = [0, 0, 0]
rms_lf = [0, 0, 0]
for year in years:
    for station in stations:
        f_out = open(dir_out+"/"+station+"_"+str(year), mode="w")
        day = UTCDateTime(year, 1, 1)
        while (day  <= UTCDateTime(year, 12, 31)):
            date = str(day.year)[-2:] + str(day.month).zfill(2) + str(day.day).zfill(2)
            day += 86400
    
            st = Stream()
            for component in components:
                fn = fn_head + date + "/" + station + "." + component
                if os.path.exists(fn):
                    st += read(fn, format=format)
            if (len(st) < 3):
                for ihour in range(24):
                    f_out.write(",,,,,\n")
                continue
            if ((abs(st[0].stats.npts - int(86400/dt)) > int(1/dt)) or
                (abs(st[1].stats.npts - int(86400/dt)) > int(1/dt)) or
                (abs(st[2].stats.npts - int(86400/dt)) > int(1/dt))):
                for ihour in range(24):
                    #f_out.write(",,,,,\n")
                    f_out.write(",,,,,,,,,\n")
                continue
            if (st[0].stats.npts != st[1].stats.npts) or (st[0].stats.npts != st[2].stats.npts):
                for ihour in range(24):
                    f_out.write(",,,,,\n")
                continue
            st_hf = st.copy().filter(type="bandpass", freqmin=2, freqmax=8, corners=2, zerophase=True)
            st_lf = st.copy().detrend("linear").filter(type="bandpass", freqmin=0.02, freqmax=0.05, corners=2, zerophase=True)
            for ihour in range(24):
                ipts0 = int(ihour*3600/dt)
                ipts1 = min(int((ihour+1)*3600/dt), st[0].stats.npts)
                for j in range(3):
                    rms_hf[j] = np.sqrt(np.mean(st_hf[j].data[ipts0:ipts1]**2))
                    rms_lf[j] = np.sqrt(np.mean(st_lf[j].data[ipts0:ipts1]**2))
                rms_hfhor = np.sqrt(np.mean(st_hf[1].data[ipts0:ipts1]**2+st_hf[2].data[ipts0:ipts1]**2))
                rms_hfall = np.sqrt(np.mean(st_hf[0].data[ipts0:ipts1]**2+st_hf[1].data[ipts0:ipts1]**2+st_hf[2].data[ipts0:ipts1]**2))
                rms_lfhor = np.sqrt(np.mean(st_lf[1].data[ipts0:ipts1]**2+st_lf[2].data[ipts0:ipts1]**2))
                rms_lfall = np.sqrt(np.mean(st_lf[0].data[ipts0:ipts1]**2+st_lf[1].data[ipts0:ipts1]**2+st_lf[2].data[ipts0:ipts1]**2))
                f_out.write(f"{rms_hf[0]:.3e},{rms_hf[1]:.3e},{rms_hf[2]:.3e},{rms_hfhor:.3e},{rms_hfall:.3e},{rms_lf[0]:.3e},{rms_lf[1]:.3e},{rms_lf[2]:.3e},{rms_lfhor:.3e},{rms_lfall:.3e}\n")
            
        f_out.close()
