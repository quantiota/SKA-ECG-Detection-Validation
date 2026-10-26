# SKA ECG Detection Validation

The SKA ECG learner reads one step per wave, P → Q → R → S → T. The waves are found in the raw
ECG by [`complex-events/ska_complex_events.py`](https://github.com/quantiota/SKA-ECG-Heart-Rate-Variability/tree/main/complex-events)
of [SKA-ECG-Heart-Rate-Variability](https://github.com/quantiota/SKA-ECG-Heart-Rate-Variability).
This repository checks that detection against waves marked by cardiologists.

## Data

[LUDB](https://physionet.org/content/ludb/1.0.1/), the Lobachevsky University Electrocardiography
Database (PhysioNet): 200 records, 10 s, 500 Hz, 12 leads. For each wave the cardiologists marked
its onset, peak and end (P, QRS, T). Lead II is used, as in the SKA ECG repositories.

The detector is used unchanged.

## Method

For every annotated wave, the detected wave of the same kind nearest to it (within 150 ms) is compared:

- **P, R, T**: error = detected peak − annotated peak.
- **Q**: is it inside [QRS onset, R peak]?
- **S**: is it inside [R peak, QRS end]?
- **False detection**: inside the annotated span of a record, a detected P, R or T with no annotated wave of the same kind within 150 ms.

## Result

| wave | annotated | matched | median error | within 20 ms | inside the annotated wave | false detections |
|---|---|---|---|---|---|---|
| P | 1,401 | 96.5 % | 2 ms | 97.2 % | 99.0 % | 3 / 1,355 * |
| R | 1,817 | 95.0 % | 2 ms | 99.0 % | — | 5 / 1,729 * |
| T | 1,642 | 90.1 % | 4 ms | 90.7 % | 94.2 % | 86 / 1,565 |
| Q | 1,726 QRS | — | — | — | 82.1 % | — |
| S | 1,726 QRS | — | — | — | 81.1 % | — |

\* In records where the wave is annotated. In 23 records no P is annotated (atrial fibrillation,
atrial flutter): there the detector still reports a P in each cycle (204). In one record no QRS
is annotated (9 R).

![Detected − annotated peak](figures/ludb_detection_errors.png)

## Files

- `validate_detection.py`: the validation
- `results/ludb_detection.csv`: one row per annotated wave
- `results/ludb_summary.csv`: the table above
- `figures/ludb_detection_errors.png`
- `logs/validate_detection.log`

## Run

```bash
git clone https://github.com/quantiota/SKA-ECG-Heart-Rate-Variability.git
git clone https://github.com/quantiota/SKA-ECG-Detection-Validation.git
cd SKA-ECG-Detection-Validation
pip install numpy pandas scipy matplotlib

# LUDB, lead II only: .hea, .dat and the .ii annotations of the 200 records
mkdir -p data && cd data
B=https://physionet.org/files/ludb/1.0.1
wget -q $B/RECORDS $B/SHA256SUMS.txt
for r in $(sed 's#data/##' RECORDS); do for e in hea dat ii; do wget -q -c $B/data/$r.$e; done; done
# check every file against PhysioNet's checksums: prints only the files that fail (download them again)
grep -E 'data/[0-9]+\.(hea|dat|ii)$' SHA256SUMS.txt | sed 's# data/# #' | sha256sum -c --quiet -
cd ..

python validate_detection.py
```

## Reference

A. I. Kalyakulina et al., "LUDB: A New Open-Access Validation Tool for Electrocardiogram
Delineation Algorithms," IEEE Access, vol. 8, pp. 186181–186190, 2020,
doi: [10.1109/ACCESS.2020.3029211](https://doi.org/10.1109/ACCESS.2020.3029211).
Data: PhysioNet, ODC Attribution License.

## License

MIT
