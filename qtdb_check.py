#!/usr/bin/env python3
"""Q and S timing against the cardiologists' marks of QTDB: old rule and new rule (detect_qs.py).

Data   QTDB (PhysioNet qtdb 1.0.0): 105 records, 15 min, 250 Hz, 2 leads; the first lead is used.
       .q1c = the first cardiologist's manual marks: QRS onset "(", peak N, QRS end ")" on selected beats.
For every annotated QRS with both its onset and its end, the detected beat whose R is nearest
(within 150 ms) is compared:
    Q inside [QRS onset, R]  and  Q − QRS onset (ms)
    S inside [R, QRS end]    and  S − QRS end (ms)
Detectors  old = ska_complex_events.detect unchanged; new = detect_qs.detect.

    python qtdb_check.py   -> results/qtdb_qs.csv, results/qtdb_summary.csv, figures/qtdb_qs.png
                              (data in data/qtdb, see README)
"""
import pathlib

import matplotlib
import numpy as np
import pandas as pd

import detect_qs

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = pathlib.Path(__file__).resolve().parent
DATA = HERE / "data" / "qtdb"
OUT = HERE
TOL = 0.150
SYM = {1: "N", 39: "(", 40: ")"}


def read_signal(rec):
    """First lead of a format-212 record, in mV, and its sampling rate."""
    lines = (DATA / f"{rec}.hea").read_text().splitlines()
    n, fs = int(lines[0].split()[1]), float(lines[0].split()[2].split("/")[0])
    f = lines[1].split()
    assert f[1] == "212", f"{rec}: format {f[1]}"
    gain, base = (f[2].split("(") + ["0)"])[:2]
    b = np.fromfile(DATA / f"{rec}.dat", dtype=np.uint8)
    b = b[: len(b) // 3 * 3].reshape(-1, 3).astype(np.int32)
    s = np.empty(2 * len(b), np.int32)
    s[0::2] = b[:, 0] | ((b[:, 1] & 0x0F) << 8)
    s[1::2] = b[:, 2] | ((b[:, 1] >> 4) << 8)
    s[s > 2047] -= 4096
    x = s[: len(s) // n * n].reshape(-1, n)[:, 0]
    g = float(gain) or 200.0                                  # WFDB: a gain of 0 means the default, 200
    return (x - float(base.rstrip(")"))) / g, fs, f[-1]


def read_annotations(rec):
    """MIT-format annotation file -> list of (sample, symbol)."""
    w = np.fromfile(DATA / f"{rec}.q1c", dtype="<u2")
    out, t, i = [], 0, 0
    while i < len(w):
        code, val = int(w[i]) >> 10, int(w[i]) & 0x3FF
        if code == 0 and val == 0:
            break
        if code == 59:
            t += (int(w[i + 1]) << 16) | int(w[i + 2]); i += 3; continue
        if code in (60, 61, 62):
            i += 1; continue
        if code == 63:
            i += 1 + (val + 1) // 2; continue
        t += val
        out.append((t, SYM.get(code, "?")))
        i += 1
    return out


def qrs(ann):
    """Annotated QRS complexes with onset, peak and end."""
    return [dict(onset=ann[i - 1][0], peak=t, end=ann[i + 1][0]) for i, (t, s) in enumerate(ann)
            if s == "N" and 0 < i < len(ann) - 1 and ann[i - 1][1] == "(" and ann[i + 1][1] == ")"]


def main():
    recs = (DATA / "RECORDS").read_text().split()
    rows = []
    for rec in recs:
        x, fs, lead = read_signal(rec)
        ann = qrs(read_annotations(rec))
        if not ann:
            continue
        dets = {"old": [ev for ev, _ in detect_qs._detect_old(x, int(fs))],
                "new": [ev for ev, _ in detect_qs.detect(x, int(fs))]}
        ms = 1000 / fs
        for name, beats in dets.items():
            R = np.array([b["R"] for b in beats])
            for a in ann:
                j = int(np.argmin(np.abs(R - a["peak"]))) if len(R) else -1
                ok = j >= 0 and abs(R[j] - a["peak"]) <= TOL * fs
                row = dict(record=rec, lead=lead, detector=name, onset=a["onset"], end=a["end"], matched=ok)
                if ok:
                    b = beats[j]
                    for v, ref in (("Q", a["onset"]), ("S", a["end"])):
                        p = b[v]
                        if p is None or (isinstance(p, float) and np.isnan(p)):
                            row[f"{v}_found"] = False
                            continue
                        row[f"{v}_found"] = True
                        row[f"{v}_inside"] = (a["onset"] <= p <= b["R"]) if v == "Q" else (b["R"] <= p <= a["end"])
                        row[f"{v}_err_ms"] = (p - ref) * ms
                rows.append(row)
        print(rec, lead, len(ann), "annotated QRS", flush=True)

    d = pd.DataFrame(rows)
    (OUT / "results").mkdir(parents=True, exist_ok=True)
    (OUT / "figures").mkdir(parents=True, exist_ok=True)
    d.to_csv(OUT / "results" / "qtdb_qs.csv", index=False, float_format="%.1f")
    s = []
    for name, g in d.groupby("detector", sort=False):
        m = g[g.matched]
        row = dict(detector=name, annotated_qrs=len(g), records=g.record.nunique(), matched_pct=100 * len(m) / len(g))
        for v in "QS":
            f = m[m[f"{v}_found"] == True]                                    # noqa: E712
            row[f"{v}_found_pct"] = 100 * len(f) / len(m)
            row[f"{v}_inside_qrs_pct"] = 100 * f[f"{v}_inside"].astype(bool).mean()
            row[f"{v}_median_ms_from_{'onset' if v == 'Q' else 'end'}"] = f[f"{v}_err_ms"].median()
        s.append(row)
    s = pd.DataFrame(s)
    s.to_csv(OUT / "results" / "qtdb_summary.csv", index=False, float_format="%.1f")
    print(s.round(1).to_string(index=False))

    cols = {"old": "#9a9893", "new": "#3987e5"}
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.2))
    for a, v, lab in ((ax[0], "Q", "Q − QRS onset (ms)"), (ax[1], "S", "S − QRS end (ms)")):
        for name, c in cols.items():
            e = d[(d.detector == name) & d.matched][f"{v}_err_ms"].dropna().clip(-120, 120)
            a.hist(e, bins=np.arange(-120, 124, 4), histtype="step", lw=2, color=c, label=name, density=True)
        a.axvline(0, color="#52514e", lw=0.6)
        a.set_xlabel(lab)
        a.set_yticks([])
        a.spines[["top", "right", "left"]].set_visible(False)
    ax[0].legend(frameon=False)
    fig.suptitle(f"QTDB, {d.record.nunique()} records, first lead: Q and S against the cardiologist's QRS onset and end")
    fig.tight_layout()
    fig.savefig(OUT / "figures" / "qtdb_qs.png", dpi=120)


if __name__ == "__main__":
    main()
