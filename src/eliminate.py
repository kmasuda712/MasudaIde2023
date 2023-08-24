from obspy import read, UTCDateTime
import csv

stations = ["N.KISF"]
rms_th = 1*10**-7
years = list(range(2004, 2022+1))
fn_detection_head = "./detection/"
fn_rms_head = "./rms/"
fn_out_head = "./detection_rms/"

for station in stations:
    starttime = UTCDateTime(years[0], 1, 1, 0, 0, 0)
    endtime = UTCDateTime(years[-1], 12, 31, 23, 59, 59)

    rms = []
    for year in years:
        with open(fn_rms_head + station + "_" + str(year), mode="r") as f:
            for row in f:
                if row[:-1].split(",")[8] == "":
                    rms.append(True)
                    continue
                if float(row[:-1].split(",")[8]) > rms_th:
                    rms.append(True)
                    continue
                rms.append(False)

    periods = []
    ndetect = 0
    sec_detects = 0
    with open(fn_detection_head + station + ".csv", mode="r") as f:
        header = f.readline()[:-1].split(",")
        for row in f:
            eventinfo = row[:-1].split(",")
            y1, mo1, d1, h1, mi1, s1 = [int(eventinfo[i]) for i in range(6)]
            y2, mo2, d2, h2, mi2, s2 = [int(eventinfo[i]) for i in range(6, 12)]
            t1 = UTCDateTime(y1, mo1, d1, h1, mi1, s1)
            t2 = UTCDateTime(y2, mo2, d2, h2, mi2, s2)
            
            irms1 = int((t1-starttime)/3600)
            irms2 = int((t2-starttime)/3600)
            if sum(rms[irms1:irms2+1]) == 0:
                periods.append([y1, mo1, d1, h1, mi1, s1, y2, mo2, d2, h2, mi2, s2])
                ndetect += 1
                sec_detects += t2 - t1
    header[2] = ndetect
    header[3] = sec_detects
    with open(fn_out_head + station + ".csv", mode="w") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(periods)
