# IF-detection-power-quality
Code and data that reproduce every signal, table and figure of the paper:

> J. M. Sanabria Villamizar, I. López García, M. Bueno López, C. S. López-Monsalvo and F. Beltran-Carbajal,
> **"Revisiting Instantaneous Frequency Detection in Power Quality: When Does Adaptive Decomposition Actually Help?"**,
> submitted to *Electric Power Systems Research* (Elsevier), 2026.

The paper separates the two stages of Hilbert–Huang-type pipelines (decomposition and instantaneous-frequency estimation) and compares 19 estimators used in measurement practice over ten Monte Carlo scenarios, a parametric sweep of the amplitude–frequency plane and a validation with a measured electric-vehicle-charger background.

Each release is archived on Zenodo with its own DOI:

- **v1.1.1**, the release that corresponds to the manuscript submitted to *Electric Power Systems Research*: latest version under <https://doi.org/10.5281/zenodo.23214230>
- v1.1.0: <https://doi.org/10.5281/zenodo.24139294>. Same code and results as v1.1.1; in v1.1.0 the file `analisis_pendientes/d_fase_aleatoria.json` still carried a note written before that run was adopted in the paper (Section 6.4 and Section S7), and its full-HHT row was labelled `HHT, naive IMF rule` instead of `HHT, band-occupancy rule`. No number changed.
- v1.0.0: <https://doi.org/10.5281/zenodo.23214231>
- all versions: <https://doi.org/10.5281/zenodo.23214230>

Table and figure numbers below are those of the submitted manuscript.

---

## Contents

The code is distributed as a single compressed archive (`4_Codigo_M3.zip`). After extracting it:

```
├── LEEME.md                Detailed usage guide (in Spanish)
├── requirements.txt        Exact library versions used for the paper
├── methods/                Estimators and test signals
├── experiments/journal/    Monte Carlo experiments, figure scripts and result files
├── experiments/            IHT experiments (preliminary version)
├── analisis_pendientes/    Analyses added in v1.1.0 (see below), with their results
├── m1m2/raw_inputs.py      Reader for the raw oscilloscope records
├── data/                   Measured waveforms and raw records
├── figures/                Figures of the paper (overwritten when regenerated)
└── extras/                 Earlier figure script and two one-off maintenance scripts
```

The folder `analisis_pendientes/` is also available uncompressed at the top level of this repository, for browsing.

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
python3 fig_key_rho.py        # Figure 1 (seconds)
python3 mc3_growing.py 5      # Scenario G with 5 realizations (quick check)
```

Most experiments take the number of Monte Carlo realizations as their first argument (default: 100). A small number is useful to check the installation; the paper uses the defaults.

> **Note.** Rerunning an experiment overwrites its result file. Copy the original results elsewhere first if you want to compare.

## Reproducing the paper

Scripts without a folder are in `experiments/journal/`; `ap/` stands for `analisis_pendientes/`. "S" numbers refer to the Supplementary Material.

| Result in the paper | Script(s) | Result file |
|---|---|---|
| Figure 1 (slope criterion) | `fig_key_rho.py` | — |
| Table 3 (regimes of each scenario) | no script: computed from Eq. (2) and Table 2 | — |
| Table 4, Tables S5–S6 (Scenario S, Cramér–Rao bound) | `mc1_stationary.py 60 15` or `mc1_parte.py <SNR> <node>`; `mc8_new_methods.py` (Matrix Pencil, Prony); `ceemdan_tabla4.py`, `ceemdan_uno.py`; IHT rows: `../exp6_iht_single_pass.py`, `../exp7_iht_deflation.py`; `mc12_crlb.py`, `analyze_round6.py` | `res_mc1_stationary.json`, `res_mc8_new_methods.json`, `res_ceemdan_tabla4.json`, `crlb.jsonl`, `res_round6.json` |
| Figure 2 and Section 6.2 (Scenario P) | `mc6_sweep.py` (all methods except CEEMDAN); `ap/p_ceemdan30.py` (CEEMDAN, 30 trials, 20 realizations per cell); `ap/fig_sweep_p30.py`; `ap/p_frontera_rf.py` (boundaries of Rilling and Flandrin) | `res_mc6_sweep.json`, `ap/p_ceemdan30.json`, `ap/p_frontera_rf.json` |
| Table 5, Table S8 (Scenario G) | `mc3_growing.py`, `mc10r_resumable.py` (CEEMDAN, FN-EMD), `mc8_new_methods.py` (FSST, Matrix Pencil, Prony) | `res_mc3_growing.json`, `res_mc10_ceemdan100.json`, `res_mc8_new_methods.json` |
| Table 6, Table S9, Figure S4 (Scenario D) | `mc4_transient.py`, `mc10r_resumable.py`, `mc8_new_methods.py`, `fig_scenarioD.py` | `res_mc4_transient.json`, `D_ce_raw.jsonl`, `res_mc8_new_methods.json` |
| Section 6.4 and S7 (Scenario D: phase of the transient, EMD seed) | `ap/d_fases.py` (72 phases), `ap/d_fase_aleatoria.py` (drawn phase), `ap/d_semilla_sin_emd.py` (same fit without EMD), `ap/s7.py` (IMF-by-IMF analysis) | `ap/d_fases.json`, `ap/d_fase_aleatoria.json`, `ap/d_semilla_sin_emd.json`, `ap/s7.json` |
| Table 7, Tables S3 and S10, Figure S5 (Scenarios U1–U4) | `mc5_unknown.py`, `mc10r_resumable.py`, `mc9r_resumable.py`, `mc8_new_methods.py`, `fig_unknown_rho.py` | `res_mc5_unknown.json`, `U1_ce_raw.jsonl`, `U2_ce_raw.jsonl`, `U3_raw.jsonl`, `U4_raw.jsonl`, `res_mc9_new_scenarios.json` |
| Table 8, Figure 3, Tables S2, S3 and S7, Figures S1–S3 (Scenarios R and N) | `mc2_rocof.py`, `mc13_rocof_stats.py`, `mc10r_resumable.py`, `mc8_new_methods.py` (P and M class, FSST), `mc17_ipdft_fix.py` (IpDFT, HT-IpDFT), `mc9_new_scenarios.py` (N), `fig_rocof_mc.py`, `fig_signal_model.py`, `fig_rocof_example.py`, `fig_scenarioU5.py` | `res_mc2_rocof.json`, `rstats.jsonl`, `res_mc10_ceemdan100.json`, `res_mc8b_pmu_exact.json`, `res_mc17_ipdft_fix_raw.json`, `res_mc9_new_scenarios.json` |
| Table 8 and Table S4 (measured background) | `mc15_semisynthetic.py`, `mc17_ipdft_fix.py`, `agg_semi.py` (see below) | `semi_R_raw.jsonl`, `semi_G_raw.jsonl`, `real_f_raw.jsonl`, `res_semi_aggregate_raw.json` |
| Section 7.3 and S4 (CEEMDAN with 200 and 500 trials, recombined IMFs) | `ap/ceemdan_ensayos.py` | `ap/ceemdan_ensayos.json` |
| Supplementary S1 (60 Hz grid) | `mc7_rocof60.py` | `res_mc7_rocof60.json` |
| Supplementary S4, Table S1 (sensitivity) | `mc14_sensitivity.py` | `sens.jsonl` |
| Supplementary S6 (numerical check of Eq. (3)) | `ap/prop1_check.py` | `ap/prop1_check.json` |
| Supplementary S8 (PMU compliance) | `mc11_pmu_compliance.py` | `res_mc11_pmu_compliance.json` |
| Supplementary S15 (VMD in Scenario N) | `ap/vmd_n.py` | `ap/vmd_n.json` |
| Table S11 (computational cost) | `timing.py` | `res_timing.json` |

Table 9 (guidance per task) and Table S12 (overall comparison) summarize the previous tables and have no script of their own.

### Analyses added in v1.1.0 (`analisis_pendientes/`)

These scripts were written during the revision of November 2026. They are run from their own folder and locate the rest of the code through the environment variable `M3_CODE` (default: the parent folder). They never overwrite the result files of `experiments/journal/`. `p_ceemdan30.py` and `ceemdan_ensayos.py` use several processes and are resumable.

### Stored results superseded by later runs

Three stored values are kept for traceability but are **not** the ones reported in the paper:

- `res_mc6_sweep.json`, column `CEEMDAN` (20 ensemble trials, 5 realizations per cell): superseded by `analisis_pendientes/p_ceemdan30.json` (30 trials, 20 realizations per cell).
- `res_mc9_new_scenarios.json`, key `U5`, row `IpDFT` (computed before the periodic Hann window was adopted): superseded by `res_mc17_ipdft_fix_raw.json`.
- `res_mc8_new_methods.json`, key `R`, row `M-class` (computed before the exact reference model of the standard): superseded by `res_mc8b_pmu_exact.json`; running `mc8_new_methods.py` with the current code reproduces the published values.

### Validation with the measured background (Tables 8 and S4)

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
             Beltran-Carbajal, Francisco},
  title   = {Revisiting Instantaneous Frequency Detection in Power Quality:
             When Does Adaptive Decomposition Actually Help?},
  journal = {Electric Power Systems Research},
  year    = {2026},
  note    = {Submitted}
}
```

## License

The code written by the authors is released under the MIT License (see `LICENSE`). Third-party libraries keep their own licenses (table above).

## Contact

Corresponding author: Irvin López García, Universidad Autónoma Metropolitana, Unidad Azcapotzalco, Mexico City, Mexico — ilg@azc.uam.mx

First author: Johinner Mauricio Sanabria Villamizar — al2212801004@azc.uam.mx
