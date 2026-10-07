# IF-detection-power-quality
Code and data that reproduce every signal, table and figure of the paper:

> J. M. Sanabria Villamizar, I. López García, M. Bueno López, C. S. López-Monsalvo and F. Beltrán Carbajal,
> **"Revisiting Instantaneous Frequency Detection in Power Quality: When Does Adaptive Decomposition Actually Help?"**,
> submitted to *Algorithms* (MDPI), 2026.

The paper separates the two stages of Hilbert–Huang-type pipelines (decomposition and instantaneous-frequency estimation) and compares 19 estimators used in measurement practice over ten Monte Carlo scenarios, a parametric sweep of the amplitude–frequency plane and a validation with a measured electric-vehicle-charger background.

<!-- After archiving the repository in Zenodo, replace this line with the DOI badge. -->

---

## Contents

The code is distributed as a single compressed archive in this repository. After extracting it:

```
├── LEEME_CODIGO.txt        Detailed usage guide (in Spanish)
├── requirements.txt        Exact library versions used for the paper
├── methods/                Estimators and test signals
├── experiments/journal/    Experiments, figure scripts and result files
├── experiments/            IHT experiments (preliminary version)
├── m1m2/raw_inputs.py      Reader for the raw oscilloscope records
├── data/                   Measured waveforms and raw records
└── figures/                Figures of the paper (overwritten when regenerated)
```

The result files of every experiment (`res_*.json`, `*.jsonl`) are included, so the tables of the paper can be checked without rerunning the experiments.

## Installation

Python 3.10 or later (the results were produced with Python 3.12.3).

```bash
pip install -r requirements.txt
```

| Library | Version | License | Used for |
|---|---|---|---|
| NumPy | 2.4.4 | BSD-3 | Numerical computing |
| SciPy | 1.17.1 | BSD-3 | Filters, signal processing |
| Matplotlib | 3.10.8 | Matplotlib (PSF-based) | Figures |
| EMD-signal (PyEMD) | 1.10.0 | Apache-2.0 | EMD, EEMD, CEEMDAN |
| vmdpy | 0.2 | MIT | VMD |
| PyWavelets | 1.10.0 | MIT | Continuous wavelet transform |

**Use these exact versions to reproduce the published values.** In a few borderline cases (errors very close to the 1 % success threshold), a different library version can change the outcome of individual realizations.

## Quick start

All scripts are run from `experiments/journal/`:

```bash
cd experiments/journal
python3 fig_key_rho.py        # Figure 2 (seconds)
python3 mc3_growing.py 5      # Scenario G with 5 realizations (quick check)
```

Most experiments take the number of Monte Carlo realizations as their first argument (default: 100). A small number is useful to check the installation; the paper uses the defaults.

> **Note.** Rerunning an experiment overwrites its result file. Copy the original results elsewhere first if you want to compare.

## Reproducing the paper

| Result in the paper | Script(s) | Result file |
|---|---|---|
| Figure 1 (signal model) | `fig_signal_model.py` | — |
| Figure 2 (slope criterion) | `fig_key_rho.py` | — |
| Table 3 (Cramér–Rao bound) | `mc12_crlb.py`, `analyze_round6.py` | `crlb.jsonl`, `res_round6.json` |
| Table 4 (Scenario S) | `mc1_stationary.py 60 15` or `mc1_parte.py <SNR> <node>`; `mc8_new_methods.py` (Matrix Pencil, Prony); IHT rows: `../exp6_iht_single_pass.py`, `../exp7_iht_deflation.py` | `res_mc1_stationary.json`, `res_mc8_new_methods.json` |
| Figure 3 and Proposition 1 (Scenario P) | `mc6_sweep.py`, `fig_sweep.py` | `res_mc6_sweep.json` |
| Table 5, Figures 4–5 (Scenario R) | `mc2_rocof.py`, `mc13_rocof_stats.py`, `mc10r_resumable.py`, `fig_rocof_example.py`, `fig_rocof_mc.py` | `res_mc2_rocof.json`, `rstats.jsonl`, `res_mc10_ceemdan100.json` |
| Table 6 (Scenario G) | `mc3_growing.py` | `res_mc3_growing.json` |
| Table 7, Figure 6 (Scenario D) | `mc4_transient.py`, `fig_D_U5.py` | `res_mc4_transient.json` |
| Table 8, Figure 7 (Scenarios U1, U2) | `mc5_unknown.py`, `mc10r_resumable.py`, `fig_unknown_rho.py` | `res_mc5_unknown.json`, `U1_ce_raw.jsonl` |
| Table 9, Figure 8 (Scenarios U3–U5) | `mc9r_resumable.py`, `fig_D_U5.py` | `res_mc9_new_scenarios.json`, `U3_raw.jsonl`, `U4_raw.jsonl` |
| Table 10 (measured background) | `mc15_semisynthetic.py`, `mc17_ipdft_fix.py`, `agg_semi.py` (see below) | `semi_R_raw.jsonl`, `semi_G_raw.jsonl`, `real_f_raw.jsonl` |
| Table 11 (computational cost) | `timing.py` | `res_timing.json` |
| Supplementary S1 (60 Hz grid) | `mc7_rocof60.py` | `res_mc7_rocof60.json` |
| Supplementary S4 (sensitivity) | `mc14_sensitivity.py` | `sens.jsonl` |
| Supplementary S8 (PMU compliance) | `mc11_pmu_compliance.py` | `res_mc11_pmu_compliance.json` |

Table 12 (overall comparison) summarizes the assessments of the previous tables and has no script of its own.

### Validation with the measured background (Table 10)

```bash
# Linux / macOS
REALBG_RAW=1 RLIB_SUFFIX=_raw python3 mc15_semisynthetic.py
REALBG_RAW=1 RLIB_SUFFIX=_raw python3 mc17_ipdft_fix.py
RLIB_SUFFIX=_raw python3 agg_semi.py
```

```powershell
# Windows (PowerShell)
$env:REALBG_RAW="1"; $env:RLIB_SUFFIX="_raw"
python mc15_semisynthetic.py
python mc17_ipdft_fix.py
python agg_semi.py
```

### Long runs

Methods based on CEEMDAN take several seconds per realization. The scripts ending in `_resumable.py`, as well as `mc12_crlb.py`, `mc13_rocof_stats.py`, `mc14_sensitivity.py`, `mc15_semisynthetic.py` and `mc17_ipdft_fix.py`, save each realization as it finishes (`.jsonl` files). If a run is interrupted, launching the same command again resumes it. Delete the `.jsonl` file to start from scratch.

## Data

| File | Content |
|---|---|
| `data/raw/V_NoLISN_2V.mat` | EV-charger voltage measured without a LISN |
| `data/raw/V_SiLISN_2V.mat` | EV-charger voltage measured with a LISN |
| `data/t_*.npy`, `data/v_*.npy` | Waveform exports of the same records (used in an earlier version) |

The raw records were acquired with a PicoScope oscilloscope at 1.25 MHz (about 2.5 s) and are resampled to 4 kHz. `m1m2/raw_inputs.py` reads them from `data/raw/`; set the environment variable `M1M2_DATA` to use a different folder.

## Reproducibility notes

- **Random seeds** are fixed in every experiment: each run produces identical results with the library versions above.
- **Success criterion:** an estimator succeeds in a realization if its relative frequency error is below 1 %.
- **EMD** is computed with the EMD-signal library. A sanity check: a pure sinusoid must return a single IMF and a negligible residue.
- **Interpolated DFT** uses the *periodic* (DFT-even) Hann window assumed by the two-point interpolation formula. The symmetric window returned by default by many libraries biases the estimate by tens of millihertz.

## Citation

If you use this code, please cite the paper:

```bibtex
@article{SanabriaVillamizar2026revisiting,
  author  = {Sanabria Villamizar, Johinner Mauricio and L{\'o}pez Garc{\'i}a, Irvin and
             Bueno L{\'o}pez, Maximiliano and L{\'o}pez-Monsalvo, C{\'e}sar Sim{\'o}n and
             Beltr{\'a}n Carbajal, Francisco},
  title   = {Revisiting Instantaneous Frequency Detection in Power Quality:
             When Does Adaptive Decomposition Actually Help?},
  journal = {Algorithms},
  year    = {2026},
  note    = {Submitted}
}
```

## License

The code written by the authors is released under the MIT License (see `LICENSE`). Third-party libraries keep their own licenses (table above).

## Contact

Corresponding author: Irvin López García, Universidad Autónoma Metropolitana, Unidad Azcapotzalco, Mexico City, Mexico — ilg@azc.uam.mx

First author: Johinner Mauricio Sanabria Villamizar — al2212801004@azc.uam.mx
