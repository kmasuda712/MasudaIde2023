from obspy import read, UTCDateTime
import numpy as np
import os
import csv

stations = ["N.KISF"]
startday = UTCDateTime(2004, 1, 1)
endday = UTCDateTime(2022, 3, 31)
threshold = 4
fn_cc_head = "./cc/"
fn_out_head = "./detection/"
dt_cc = 10
twin_mvave = 10000

arr_mvave = np.ones(int(twin_mvave/dt_cc)) / int(twin_mvave/dt_cc)
for station in stations:
    day = startday
    cc = np.array([])
    cc_std = np.array([])
    nday1 = 0
    nday2 = 0
    while (day <= endday):
        nday1 += 1
        date = str(day.year)[-2:] + str(day.month).zfill(2) + str(day.day).zfill(2)
        day += 86400
        
        fn = fn_cc_head + station + "/" + date
        if os.path.exists(fn):
            cc_tmp = np.loadtxt(fn)
            cc = np.concatenate([cc, cc_tmp[:,1]])
            cc_std = np.concatenate([cc_std, cc_tmp[:,1]])
            nday2 += 1
        else:
            cc = np.concatenate([cc, np.zeros(int(86400/dt_cc))])
    
    np.nan_to_num(cc, copy=False, nan=0.0, posinf=0.0, neginf=0.0)
    np.nan_to_num(cc_std, copy=False, nan=0.0, posinf=0.0, neginf=0.0)
    
    
    cc_mvave = np.convolve(cc, arr_mvave, mode="same")
    posinega = True
    if (len(cc_std) > 0):
        cc_std_mvave = np.convolve(cc_std, arr_mvave, mode="same")
        std = np.std(cc_std_mvave)
        if (np.count_nonzero(cc_std_mvave>std) <
                np.count_nonzero(cc_std_mvave<-std)):
            posinega = False
            cc = -cc
            cc_mvave = -cc_mvave
    else:
        std = 1
    
    pflag = False
    psacday = UTCDateTime(0)
    periods = []
    ndetect = 0
    sec_detects = 0
    for itime in range(cc.size):
        flag = cc_mvave[itime] > threshold*std
        if not pflag and flag:
            itime_start = itime
        if pflag and not flag:
            itime_end = itime
            time_start = startday + dt_cc*itime_start
            time_end = startday + dt_cc*itime_end
            periods.append([time_start.year, time_start.month, time_start.day, time_start.hour, time_start.minute, time_start.second,
                            time_end.year, time_end.month, time_end.day, time_end.hour, time_end.minute, time_end.second])
            ndetect += 1
            sec_detects += time_end - time_start
        pflag = flag
    info = [nday1, nday2, ndetect, sec_detects, std, posinega, startday.year, startday.month, startday.day, endday.year, endday.month, endday.day]
    
    with open(fn_out_head + station + ".csv", 'w') as f:
        writer = csv.writer(f)
        writer.writerow(info)
        writer.writerows(periods)
