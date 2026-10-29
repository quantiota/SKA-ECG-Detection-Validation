#!/usr/bin/env python3
"""Why the R detector misses beats on QTDB: the QRS energy of each missed beat against the threshold.

The R detector of ska_complex_events.r_peaks (unchanged) keeps an energy peak if it is above 40 % of
the maximum of its 5-s window and at least 300 ms from the previous one. The same energy is
recomputed here (same filters), and for every annotated QRS the new rule missed (qtdb_check.py):
    energy_ratio      = the beat's energy peak / the threshold (< 1: below the 40 % threshold)
    gap_prev_ms       = distance to the previous kept peak (< 300: removed by the refractory distance)
    window_max_ratio  = the window maximum / the median energy of the record's beats (a large beat
                        or an artefact raising the threshold)

    python qtdb_missed.py   -> results/qtdb_missed.csv
"""
import pathlib

import numpy as np
import pandas as pd
from scipy.signal import butter, filtfilt, find_peaks

import detect_qs
from qtdb_check import read_signal

HERE = pathlib.Path(__file__).resolve().parent
sce = detect_qs.sce


def energy(x, fs):
    y = sce.edge_signal(x, fs, 50.0)
    b, a = butter(3, [5.0, 20.0], btype="band", fs=fs)
    y = filtfilt(b, a, y)
    return np.convolve(np.gradient(y) ** 2, np.ones(int(0.15 * fs)) / (0.15 * fs), "same")


def main():
    d = pd.read_csv(HERE / "results" / "qtdb_qs.csv")
    miss = d[(d.detector == "new") & ~d.matched]
    rows = []
    for rec, g in miss.groupby("record"):
        x, fs, lead = read_signal(rec)
        fs = int(fs)
        e = energy(x, fs)
        win = 5 * fs
        local_max = np.maximum.reduceat(e, np.arange(0, len(e), win))
        thr = 0.4 * np.repeat(local_max, win)[:len(e)]
        kept, _ = find_peaks(e, height=thr, distance=int(0.3 * fs))
        allpk, _ = find_peaks(e)
        med = np.median(e[kept])
        for a in g.itertuples():
            c = (a.onset + a.end) // 2
            near = allpk[np.abs(allpk - c) <= 0.1 * fs]
            p = near[np.argmax(e[near])] if len(near) else c
            prev = kept[kept < p]
            rows.append(dict(record=rec, lead=lead, qrs_sample=c, energy_ratio=e[p] / thr[p],
                             gap_prev_ms=(p - prev[-1]) * 1000 / fs if len(prev) else np.nan,
                             window_max_ratio=local_max[p // win] / med))
    r = pd.DataFrame(rows)
    r.to_csv(HERE / "results" / "qtdb_missed.csv", index=False, float_format="%.2f")
    print(r.groupby(["record", "lead"]).agg(missed=("qrs_sample", "size"),
                                            below_threshold=("energy_ratio", lambda v: int((v < 1).sum())),
                                            median_energy_ratio=("energy_ratio", "median"),
                                            median_gap_prev_ms=("gap_prev_ms", "median"),
                                            median_window_max_ratio=("window_max_ratio", "median")).round(2).to_string())


if __name__ == "__main__":
    main()
