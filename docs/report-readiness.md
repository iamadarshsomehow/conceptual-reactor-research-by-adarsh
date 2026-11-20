# Report Readiness Notes

This note reconciles the report's conclusions and run details with the models and files currently included in this repository.

## 1. Interpret the results within the model scope

The neutronics calculation reports `k-infinity = 1.05892 ± 0.00053` for a reflected single block. Applying the project's assumed 3% leakage fraction gives the arithmetic estimate `k-effective ≈ 1.027`. The 3% leakage is an assumption, not a value calculated by OpenMC. Accordingly, this repository characterizes `1.027` as a leakage-corrected screening estimate; it does not demonstrate criticality of a complete finite reactor core.

The heavy-water drain calculations show reduced reactivity in the modeled cases. They do not establish that a real passive drainage or shutdown system would operate as intended. The thermal-hydraulics calculations are one-dimensional estimates for stated assumptions, and the power map is passed one way from neutronics without temperature feedback. These results are not a validated safety analysis or proof of reactor performance.

The RBMK comparison uses simplified single-cell models. Its trends should not be presented as a full-core comparison or code-validation benchmark.

## 2. Run settings and result provenance

The report states that its single-block runs used 60 batches, 15 inactive batches, and 50,000 particles per batch. The current repository scripts document different settings:

| Current script or run mode | Batches | Inactive batches | Particles per batch |
|---|---:|---:|---:|
| `neutronics/src/model.py` baseline | 35 | 5 | 50,000 |
| `neutronics/workflows/sweep_all.py` full run | 100 | 30 | 50,000 |
| `neutronics/workflows/sweep_all.py --quick` | 60 | 15 | 20,000 |
| `neutronics/workflows/drain_runs.py` full run | 200 | 30 | 100,000 |
| `neutronics/workflows/drain_runs.py --quick` | 60 | 15 | 20,000 |

The report's stated combination of 60 batches, 15 inactive batches, and 50,000 particles does not match these current scripts. Raw statepoints are not included and no fixed random seed is set, so the repository cannot identify which exact run produced the report's quoted baseline value. The report settings should therefore be treated as reported run metadata that has not been matched to a retained statepoint; they should not be silently replaced with settings from a different script.

## 3. Reproduction limits

The repository includes STEP geometry, source scripts, compact result tables, and selected figures. It does not include the generated DAGMC `.h5m` geometry or raw OpenMC statepoints and tallies. The original material-assignment and CAD conversion commands were not fully preserved, so the exact DAGMC input cannot be rebuilt from this repository alone.

The thermal-hydraulics input data and scripts are included, but `run_cases.py` and `power_scan.py` expect input paths that differ from the repository layout and require path adjustments. The RBMK script is archived, but uses a hard-coded local geometry path; its full sweep and post-processing workflow are not included. Python package lists are unpinned, and a nuclear-data checksum is not recorded.

As currently packaged, the repository supports review of the methods and saved results, but it is not a turnkey or independently verified reproduction package.
