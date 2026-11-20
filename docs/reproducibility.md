# Reproducibility

This guide records the software environment, geometry inputs, workflow commands, and outputs for the conceptual neutronics and thermal-hydraulics study. Full neutronics reproduction requires building DAGMC geometry from the included STEP files and obtaining the ENDF/B-VII.1 HDF5 cross-section library.

## Software versions

| Component | Version or configuration |
|---|---|
| Operating system | Ubuntu 24.04 on WSL2 |
| Compiler and build tools | GCC 13.3.0, CMake 3.28.3, Git 2.43.0 |
| OpenMC | 0.16.1-dev8, commit `7ecd3a9`, built from source with DAGMC support; MPI, parallel HDF5, libMesh, and UWUW disabled |
| MOAB | 5.6.0 |
| DAGMC | Approximately 3.2; exact release identifier remains to be recorded |
| Python | 3.12.3, virtual environment named `nuclear-env` |
| Python packages | See `requirements/neutronics.txt` and `requirements/thermal-hydraulics.txt`; versions are not pinned in these files yet |
| Nuclear data | ENDF/B-VII.1 HDF5 library for OpenMC |
| CAD | Onshape, exported as STEP |
| CAD-to-DAGMC conversion | `cad_to_openmc`, STL2 backend; exact package version remains to be recorded |

## Environment setup

Activate the Python environment and point OpenMC to the cross-section index file:

```bash
source ~/nuclear-env/bin/activate
export OPENMC_CROSS_SECTIONS=/path/to/endfb-vii.1-hdf5/cross_sections.xml
```

Install the native OpenMC, DAGMC, and MOAB stack using the versions/configuration above. Install the Python packages listed in the two requirements files. The exact package lock file can be added when a final `pip freeze` from the simulation environment is captured.

## Geometry and required inputs

The repository shares the CAD geometry as STEP files. The final design STEP is `development/versions/somuv4/somuv4.step`; the RBMK-type comparison STEP is `development/versions/rbmkagain/benchmark.step`. Earlier versions and other STEP exports are preserved under `development/versions/` and `development/geometry/`.

OpenMC consumes DAGMC `.h5m` files produced from the STEP geometry through the CAD-to-DAGMC conversion workflow. The generated `somuv4.h5m` and comparison `benchmark.h5m` files are not committed; build them locally from the corresponding STEP files. The design model expects the following DAGMC material tags to match the OpenMC material names exactly:

`air`, `carrier`, `clad`, `d20`, `fuel`, `graphite`, `h20`, `pressure`

The RBMK comparison model uses its own geometry/material tags; check the corresponding comparison model source when that workflow is added to this repository.

The thermal-hydraulics input power map is `thermal_hydraulics/data/power_map.csv`. The ENDF/B-VII.1 library requires `cross_sections.xml` plus the referenced HDF5 data files. Check that the downloaded library archive is complete before extraction.

The neutronics scripts accept the design geometry path through `SOMUV4_H5M`. The default local path is `private/somuv4.h5m`, which is ignored by Git. Set `SOMUV4_H5M` to your generated `.h5m` file if it is stored elsewhere.

## Neutronics workflows

Run commands from the repository root:

| Script | Command | Full-run settings | Main outputs |
|---|---|---|---|
| Baseline transport model | `python neutronics/src/model.py` | 35 batches, 5 inactive, 50,000 particles per batch | OpenMC statepoint, `tallies.out`, XML model/settings files |
| Baseline, void and moderator-density sweep | `python neutronics/workflows/sweep_all.py` | 100 batches, 30 inactive, 50,000 particles per batch | `sweep_results.csv`, `void_results.png`, case folders and baseline mesh power map |
| Baseline, D₂O drain and LOCA cases | `python neutronics/workflows/drain_runs.py [baseline|drain|loca]` | 200 batches, 30 inactive, 100,000 particles per batch | `results_*.csv`, `results_*.md`, `results_*_table.png`, `results_*_plot.png`, `entropy_*.csv`, `entropy_*_plot.png` |
| Rebuild drain/LOCA reports from saved CSVs | `python neutronics/workflows/drain_runs.py plot` | No transport run | Tables and figures from existing CSV files |
| Void-effect chart | `python development/versions/somuv4/make_void_chart.py` | No transport run | `void_effect_chart.png` |

The two sweep workflows support `--quick` for short test runs. `sweep_all.py --quick` uses 60 batches, 15 inactive, and 20,000 particles per batch. `drain_runs.py --quick` uses the same settings. No fixed random seed is set in these scripts.

The RBMK comparison scripts `rbmk_benchmark_sweep.py` and `normalized_void.py` are not included in the current main workflow folders, so their run commands and outputs are not listed as reproducible workflows here.

## Thermal-hydraulics workflows

The thermal-hydraulics scripts are in `thermal_hydraulics/src/`; the configuration is in `thermal_hydraulics/config/config.yaml`, and the power map is in `thermal_hydraulics/data/power_map.csv`. In the current checkout, `run_cases.py` and `power_scan.py` expect `config.yaml` and `power_map.csv` alongside the scripts, while the repository stores them in separate folders. Their paths need to be connected before those two orchestration scripts can be run directly from this checkout. `channel.py` can be given the configuration path explicitly, but the power-map path in that configuration must resolve to the data file.

| Script | Purpose | Main outputs |
|---|---|---|
| `channel.py` | `python thermal_hydraulics/src/channel.py thermal_hydraulics/config/config.yaml` | `channel_results.csv`, `channel_energy_balance.png`, `summary.json` |
| `run_cases.py` | `python thermal_hydraulics/src/run_cases.py` (after resolving the input paths noted above) | `cases_summary.csv`, `cases_comparison.png`, per-case folders |
| `power_scan.py` | `python thermal_hydraulics/src/power_scan.py` (after resolving the input paths noted above) | `scan_all.csv`, `allowable_power.csv`, `power_scan.png`, per-case scan folders |

The thermal workflow uses 50 axial nodes for the base configuration; `run_cases.py` also includes a 100-node sensitivity case. The neutronics power map is passed one way into the thermal model, with no temperature feedback into OpenMC.

## Runtime

The approximate full-workflow runtime reported for a Ryzen 7 7735HS system with 16 GB LPDDR5X-6400 memory is 40 minutes to 1 hour. Runtime varies with the selected cases, build configuration, and system load.

## Regenerated outputs

OpenMC statepoints, `summary.h5`, `tallies.out`, XML runtime files, VTK files, and case/run folders are generated by the workflows and are not committed. Compact result tables, selected plots, the final STEP models, and the version history are included in the repository.
