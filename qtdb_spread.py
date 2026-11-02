#!/usr/bin/env python3
"""T peak and the spread of the wave intervals against the cardiologists' marks of QTDB (old rule).

Data   QTDB (PhysioNet qtdb 1.0.0): 105 records, 15 min, 250 Hz, 2 leads; the first lead is used.
       .q1c = the first cardiologist's manual marks on selected beats: P peak "p", R peak (beat label),
       T peak "t". A beat's P is the last "p" before its R, its T the first "t" after it.
For every marked beat, the detected beat whose R is nearest (within 150 ms) is compared:
    T error   detected T peak − marked T peak (ms), split by the sign of the marked T on the signal
              minus the detector's baseline (median of [R − 300, R − 220] ms)
    intervals P→R, R→T and T→P (T to the P of the next beat, when the next beat is also marked and
              detected), marked and detected on the same beats
    spread    per record and interval (at least 10 beats): SD and IQR, marked and detected
R→T and T→P are taken on the beats whose marked T is positive (the detector takes T as the highest
point of [R + 120 ms, R + 0.6 RR]); P→R on all beats.

    python qtdb_spread.py   -> results/qtdb_spread_beats.csv, results/qtdb_spread.csv,
                               results/qtdb_spread_summary.csv, figures/qtdb_spread.png
                               (data in data/qtdb, see README)
"""
import matplotlib
import numpy as np
import pandas as pd

import detect_qs
import qtdb_check as q

matplotlib.use("Agg")
import matplotlib.pyplot as plt

SYM = {**{c: "N" for c in range(1, 14)}, 24: "p", 27: "t"}      # beat labels 1-13 = R; p, t = P and T peaks
INTERVALS = ["P→R", "R→T", "T→P"]


def marked_beats(ann):
    """Marked beats: R with the last P before it and the first T after it (None if not marked)."""
    N = [i for i, (_, s) in enumerate(ann) if s == "N"]
    out = []
    for k, i in enumerate(N):
        lo, hi = N[k - 1] if k else -1, N[k + 1] if k + 1 < len(N) else len(ann)
        p = [t for t, s in ann[lo + 1:i] if s == "p"]
        t = [t for t, s in ann[i + 1:hi] if s == "t"]
        out.append(dict(P=p[-1] if p else None, R=ann[i][0], T=t[0] if t else None))
    return out


def main():
    q.SYM.clear(); q.SYM.update(SYM)
    rows = []
    for rec in (q.DATA / "RECORDS").read_text().split():
        x, fs, lead = q.read_signal(rec)
        beats = marked_beats(q.read_annotations(rec))
        if not beats:
            continue
        D = [ev for ev, _ in detect_qs._detect_old(x, int(fs))]
        DR = np.array([d["R"] for d in D])
        ms = 1000 / fs
        for k, b in enumerate(beats):
            j = int(np.argmin(np.abs(DR - b["R"])))
            if abs(DR[j] - b["R"]) > q.TOL * fs:
                continue
            d, row = D[j], dict(record=rec, lead=lead, R=b["R"])
            if b["P"] is not None:
                row.update(PR_marked=(b["R"] - b["P"]) * ms, PR_detected=(d["R"] - d["P"]) * ms)
            if b["T"] is not None:
                base = np.median(x[b["R"] - int(0.3 * fs): b["R"] - int(0.22 * fs)])
                row.update(T_sign=int(np.sign(x[b["T"]] - base)), T_err_ms=(d["T"] - b["T"]) * ms,
                           RT_marked=(b["T"] - b["R"]) * ms, RT_detected=(d["T"] - d["R"]) * ms)
                nb, dn = (beats[k + 1], D[j + 1]) if k + 1 < len(beats) and j + 1 < len(D) else (None, None)
                if nb and nb["P"] is not None and abs(dn["R"] - nb["R"]) <= q.TOL * fs:
                    row.update(TP_marked=(nb["P"] - b["T"]) * ms, TP_detected=(dn["P"] - d["T"]) * ms)
            rows.append(row)
        print(rec, lead, len(beats), "marked beats", flush=True)

    d = pd.DataFrame(rows)
    d.to_csv(q.OUT / "results" / "qtdb_spread_beats.csv", index=False, float_format="%.1f")

    print("\nT error (detected − marked), ms")
    for sign, g in d.dropna(subset=["T_err_ms"]).groupby("T_sign"):
        e = g.T_err_ms
        print(f"  marked T {'positive' if sign > 0 else 'negative' if sign < 0 else 'on the baseline'}: {len(e)} beats, "
              f"{g.record.nunique()} records, "
              f"median {e.median():.0f}, within 20 ms {100 * (e.abs() <= 20).mean():.1f} %")

    iqr = lambda v: np.subtract(*np.percentile(v, [75, 25]))
    s = []
    for rec, g in d.groupby("record", sort=False):
        for name, c in zip(INTERVALS, ["PR", "RT", "TP"]):
            m = (g if c == "PR" else g[g.T_sign > 0])[[f"{c}_marked", f"{c}_detected"]].dropna()
            if len(m) < 10:
                continue
            a, b = m.iloc[:, 0], m.iloc[:, 1]
            s.append(dict(record=rec, interval=name, beats=len(m), median_marked=a.median(), median_detected=b.median(),
                          sd_marked=a.std(), sd_detected=b.std(), iqr_marked=iqr(a), iqr_detected=iqr(b)))
    s = pd.DataFrame(s)
    s.to_csv(q.OUT / "results" / "qtdb_spread.csv", index=False, float_format="%.2f")
    summ = []
    for name in INTERVALS:
        g = s[s.interval == name]
        summ.append(dict(interval=name, records=len(g), beats_median=g.beats.median(),
                         **{c: g[c].median() for c in ["median_marked", "median_detected", "sd_marked", "sd_detected",
                                                        "iqr_marked", "iqr_detected"]},
                         sd_correlation=np.corrcoef(g.sd_marked, g.sd_detected)[0, 1]))
    summ = pd.DataFrame(summ)
    summ.to_csv(q.OUT / "results" / "qtdb_spread_summary.csv", index=False, float_format="%.2f")
    print("\n" + summ.round(2).to_string(index=False))

    fig, ax = plt.subplots(1, 3, figsize=(13, 4.2))
    for a, name in zip(ax, INTERVALS):
        g = s[s.interval == name]
        lo, top = 1, 1.3 * max(g.sd_marked.max(), g.sd_detected.max())         # log axes: a few records reach
        a.plot([lo, top], [lo, top], color="#bbbbbb", lw=0.8)                  # hundreds of ms
        a.scatter(g.sd_marked.clip(lower=lo), g.sd_detected.clip(lower=lo), s=14, color="#2a78d6")
        a.set_xscale("log"); a.set_yscale("log")
        a.set_xlim(lo, top); a.set_ylim(lo, top)
        a.set_title(f"{name}: {len(g)} records", fontsize=10, loc="left")
        a.set_xlabel("SD of the marked interval, ms")
    ax[0].set_ylabel("SD of the detected interval, ms")
    fig.suptitle("QTDB: spread of the wave intervals per record, cardiologist marks vs SKA detector (same beats)",
                 fontsize=10, x=0.01, ha="left")
    fig.tight_layout()
    fig.savefig(q.OUT / "figures" / "qtdb_spread.png", dpi=120)
    print("saved results/qtdb_spread_beats.csv, results/qtdb_spread.csv, results/qtdb_spread_summary.csv, "
          "figures/qtdb_spread.png")


if __name__ == "__main__":
    main()
