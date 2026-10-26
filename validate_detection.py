#!/usr/bin/env python3
"""Validate the P, Q, R, S, T detection of complex-events/ska_complex_events.py against the
cardiologists' annotations of the PhysioNet LUDB (lead II, 500 Hz, 10 s, 200 records).

LUDB marks, per wave, its onset '(' , its peak (p = P, N = QRS peak, t = T) and its end ')'.
For every annotated wave, the detected peak of the same beat (the one nearest to it, within
150 ms) is compared:
    P, R, T : error = detected peak − annotated peak (ms)
    Q       : inside [QRS onset, R peak]       S : inside [R peak, QRS end]
False detections: inside the annotated span of a record (first to last annotation), a detected
P, R or T peak with no annotated wave of the same kind within 150 ms.
The detector is imported from SKA-ECG-Heart-Rate-Variability unchanged.

    python validate_detection.py
        -> results/ludb_detection.csv, results/ludb_summary.csv, figures/ludb_detection_errors.png
"""
import pathlib
import sys

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "SKA-ECG-Heart-Rate-Variability" / "complex-events"))
from ska_complex_events import detect                                    # noqa: E402

DATA, LEAD, TOL = HERE / "data", "ii", 0.150          # lead II, match within 150 ms
SYM = {1: "N", 24: "p", 27: "t", 39: "(", 40: ")"}     # MIT annotation codes used by LUDB


def read_signal(rec):
    lines = (DATA / f"{rec}.hea").read_text().splitlines()
    n, fs = int(lines[0].split()[1]), float(lines[0].split()[2])
    names = [l.split()[-1] for l in lines[1:n + 1]]
    j = names.index(LEAD)
    gain, base = lines[1 + j].split()[2].split("/")[0].rstrip(")").split("(")
    raw = np.fromfile(DATA / f"{rec}.dat", dtype="<i2").reshape(-1, n)[:, j]
    return (raw - float(base)) / float(gain), fs


def read_annotations(rec):
    """MIT-format annotation file -> list of (sample, symbol)."""
    w = np.fromfile(DATA / f"{rec}.{LEAD}", dtype="<u2")
    out, t, i = [], 0, 0
    while i < len(w):
        code, val = int(w[i]) >> 10, int(w[i]) & 0x3FF
        if code == 0 and val == 0:
            break
        if code == 59:                                   # SKIP: 32-bit interval follows
            t += (int(w[i + 1]) << 16) | int(w[i + 2]); i += 3; continue
        if code in (60, 61, 62):                         # NUM, SUB, CHN
            i += 1; continue
        if code == 63:                                   # AUX: val bytes follow
            i += 1 + (val + 1) // 2; continue
        t += val
        out.append((t, SYM.get(code, "?")))
        i += 1
    return out


def waves(ann):
    """Annotated waves: dicts with kind (P, QRS, T), onset, peak, end (samples)."""
    res = []
    for i, (t, s) in enumerate(ann):
        if s in "pNt" and 0 < i < len(ann) - 1 and ann[i - 1][1] == "(" and ann[i + 1][1] == ")":
            res.append(dict(kind={"p": "P", "N": "QRS", "t": "T"}[s], onset=ann[i - 1][0], peak=t, end=ann[i + 1][0]))
    return res


def main():
    recs = [r.split("/")[-1] for r in (DATA / "RECORDS").read_text().split()]
    rows, det = [], []
    for rec in recs:
        if not (DATA / f"{rec}.{LEAD}").exists():
            continue
        x, fs = read_signal(rec)
        beats = [ev for ev, _ in detect(x, int(fs))]
        tol = TOL * fs
        ann = read_annotations(rec)
        wvs = waves(ann)
        lo, hi = ann[0][0], ann[-1][0]
        for kind, key in (("P", "P"), ("QRS", "R"), ("T", "T")):
            peaks = np.array([w["peak"] for w in wvs if w["kind"] == kind])
            for ev in beats:
                if lo <= ev[key] <= hi:
                    det.append(dict(record=rec, wave=key, annotated_in_record=bool(len(peaks)),
                                    false=not len(peaks) or np.min(np.abs(peaks - ev[key])) > tol))
        for wv in wvs:
            if wv["kind"] == "QRS":
                ref = np.array([ev["R"] for ev in beats])
            else:
                ref = np.array([ev[wv["kind"]] for ev in beats])
            if not len(ref) or np.min(np.abs(ref - wv["peak"])) > tol:
                rows.append(dict(record=rec, wave=wv["kind"], matched=False)); continue
            ev = beats[int(np.argmin(np.abs(ref - wv["peak"])))]
            ms = 1000 / fs
            if wv["kind"] == "QRS":
                rows.append(dict(record=rec, wave="R", matched=True, error_ms=(ev["R"] - wv["peak"]) * ms))
                rows.append(dict(record=rec, wave="Q", matched=True, inside=wv["onset"] <= ev["Q"] <= ev["R"],
                                 error_ms=(ev["Q"] - wv["onset"]) * ms))
                rows.append(dict(record=rec, wave="S", matched=True, inside=ev["R"] <= ev["S"] <= wv["end"],
                                 error_ms=(ev["S"] - wv["end"]) * ms))
            else:
                rows.append(dict(record=rec, wave=wv["kind"], matched=True, error_ms=(ev[wv["kind"]] - wv["peak"]) * ms,
                                 inside=wv["onset"] <= ev[wv["kind"]] <= wv["end"]))
    d = pd.DataFrame(rows)
    (HERE / "results").mkdir(exist_ok=True); (HERE / "figures").mkdir(exist_ok=True)
    d.to_csv(HERE / "results" / "ludb_detection.csv", index=False)

    det = pd.DataFrame(det)
    s = []
    for wave in "PQRST":
        g = d[d.wave == wave] if wave in "QS" else d[d.wave.isin([wave, "QRS" if wave == "R" else wave])]
        m = g[g.matched]
        e = m.error_ms.abs()
        s.append(dict(wave=wave, annotated=len(g) if wave not in "QS" else len(m), matched=len(m),
                      matched_pct=round(100 * len(m) / max(len(g), 1), 1),
                      median_abs_error_ms=None if wave in "QS" else round(e.median(), 1),
                      within_20ms_pct=None if wave in "QS" else round(100 * (e <= 20).mean(), 1),
                      inside_wave_pct=None if wave == "R" else round(100 * m.inside.astype(bool).mean(), 1),
                      detected_in_span=None if wave in "QS" else int((det.wave == wave).sum()),
                      false_detections=None if wave in "QS" else int(det[det.wave == wave].false.sum()),
                      false_pct=None if wave in "QS" else round(100 * det[det.wave == wave].false.mean(), 1)))
    s = pd.DataFrame(s)
    s.to_csv(HERE / "results" / "ludb_summary.csv", index=False)
    print(f"{len(set(d.record))} records, lead {LEAD}")
    print(s.to_string(index=False))
    for wave in "PRT":
        g = det[det.wave == wave]
        a, b = g[g.annotated_in_record], g[~g.annotated_in_record]
        print(f"{wave} false detections: {int(a.false.sum())}/{len(a)} in records where {wave} is annotated, "
              f"{int(b.false.sum())}/{len(b)} in the {b.record.nunique()} records where it is not")

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.5), sharey=True)
    for a, wave in zip(ax, "PRT"):
        e = d[(d.wave == wave) & d.matched].error_ms
        a.hist(e, bins=np.arange(-150, 152, 4), color="#3987e5")
        a.axvline(0, color="black", lw=.6); a.set_title(f"{wave}: detected − annotated peak")
        a.set_xlabel("ms")
    ax[0].set_ylabel("waves")
    fig.suptitle(f"LUDB lead II, {len(set(d.record))} records: SKA wave detection vs the cardiologists' annotations")
    fig.tight_layout()
    fig.savefig(HERE / "figures" / "ludb_detection_errors.png", dpi=110)
    print(f"saved {HERE / 'figures' / 'ludb_detection_errors.png'}")


if __name__ == "__main__":
    main()
