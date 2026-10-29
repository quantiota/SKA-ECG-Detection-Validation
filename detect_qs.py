"""New Q and S rule for the SKA detector, applied in memory (ska_complex_events.py stays unchanged).

Old rule   Q = lowest point in [R − 60, R] ms,  S = lowest point in [R, R + 60] ms
New rule   on the same baseline-removed signal xb, walking from R (at most 120 ms):
    Q  walk back while the signal falls; Q = the foot of the R upstroke, the first minimum.
       No baseline condition: the PR segment can sit below the baseline.
    S  walk forward; S = the first minimum below the baseline. If none is below the baseline
       within 120 ms, S = the first minimum (the sign check of complex_events then marks it absent).
    A minimum ends when the signal rises from it for at least 3 ms or by more than 0.005 mV,
    so a one- or two-sample ripple on the slope of R does not stop the walk.
    The search for a minimum starts only once the walk is below half of the R height, so a notch
    on the top of R is not taken for Q or S. If the walk never goes below half of R within 120 ms,
    the old position is kept.

    import detect_qs; detect_qs.patch()   # ska_complex_events.detect now uses the new rule
"""
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "SKA-ECG-Heart-Rate-Variability" / "complex-events"))
import ska_complex_events as sce                                   # noqa: E402

LIMIT_MS, RISE_MS, RISE_MV, HALF = 120, 3, 0.005, 0.5
_detect_old = sce.detect


def minima(xb, r, step, fs):
    """Successive minima met walking from r in direction step (−1 back, +1 forward), within LIMIT_MS."""
    end = r + step * int(LIMIT_MS * fs / 1000)
    end = min(max(end, 0), len(xb) - 1)
    need = max(1, round(RISE_MS * fs / 1000))
    i = r                                             # skip the top of R (notches): start below half of R
    while i != end and xb[i] > HALF * xb[r]:
        i += step
    if i == end:
        return []                                     # R never halves: keep the old position
    m = M = i
    n, falling, out = 0, True, []
    while i != end:
        i += step
        if falling:                                   # tracking a minimum
            if xb[i] < xb[m]:
                m, n = i, 0
            else:
                n += 1
                if xb[i] - xb[m] > RISE_MV or n >= need:   # a real rise: the minimum is confirmed
                    out.append(m)
                    falling, M, n = False, i, 0
        else:                                         # rising again: wait for the next fall
            if xb[i] > xb[M]:
                M, n = i, 0
            else:
                n += 1
                if xb[M] - xb[i] > RISE_MV or n >= need:
                    falling, m, n = True, i, 0
    return out or [m]


def q_peak(xb, r, fs, old):
    ms = minima(xb, r, -1, fs)
    return ms[0] if ms else old


def s_peak(xb, r, fs, old):
    ms = minima(xb, r, +1, fs)
    below = [m for m in ms if xb[m] < 0]
    return below[0] if below else ms[0] if ms else old


def detect(x, fs=1000):
    """ska_complex_events.detect with the new Q and S; P, R, T and the baseline unchanged."""
    beats = _detect_old(x, fs)
    for ev, base in beats:
        xb = x - base
        ev["Q"], ev["S"] = q_peak(xb, ev["R"], fs, ev["Q"]), s_peak(xb, ev["R"], fs, ev["S"])
    return beats


def patch():
    sce.detect = detect
