# Reproducibility

This guide records the software environment, geometry inputs, workflow commands, and outputs for the conceptual neutronics and thermal-hydraulics study. Full neutronics reproduction requires building DAGMC geometry from the included STEP files and obtaining the ENDF/B-VII.1 HDF5 cross-section library.

## Software versions

| Component | Version or configuration |
|---|---|
| Operating system | Ubuntu 24.04 on WSL2 |
| Compiler and build tools | GCC 13.3.0, CMake 3.28.3, Git 2.43.0 |
| OpenMC | 0.16.1.dev8, commit `7ecd3a9613f06ed0dc22368c9540faf4aaacc65f`, built from source with DAGMC; Python API matches the binary |
| OpenMC build options | MPI, parallel HDF5, libMesh, and UWUW disabled |
| DAGMC | 3.2.4 |
| MOAB | 5.6.0 |
| Python for neutronics | 3.12.3, virtual environment named `nuclear-env` |
| Python for CAD conversion | Python 3.11, separate environment named `cad-to-openmc-env` |
| Python packages | See `requirements/neutronics.txt` and `requirements/thermal-hydraulics.txt`; versions are not pinned in these files yet |
| Nuclear data | ENDF/B-VII.1 HDF5 library from the [official OpenMC data page](https://openmc.org/data/); 444 files in the `neutron/` directory |
| CAD | Onshape, exported as STEP |
| CAD-to-DAGMC conversion | `CAD_to_OpenMC` 1.0.1, STL2 backend; separate Python 3.11 environment |

## Environment setup

Activate the Python environment and point OpenMC to the cross-section index file:

```bash
source ~/nuclear-env/bin/activate
export OPENMC_CROSS_SECTIONS=/path/to/endfb-vii.1-hdf5/cross_sections.xml
```

Install the native OpenMC, DAGMC, and MOAB stack using the versions/configuration above. The Python package lists are unpinned, and no package lock file or final `pip freeze` from the simulation environments is included.

The CAD conversion environment is separate from the OpenMC environment:

```bash
python3.11 -m venv cad-to-openmc-env
source cad-to-openmc-env/bin/activate
pip install CAD_to_OpenMC==1.0.1
```

## Geometry and required inputs

The repository shares the CAD geometry as STEP files. The final design STEP is `development/versions/somuv4/somuv4.step`; the RBMK-type comparison STEP is `development/versions/rbmkagain/benchmark.step`. Earlier versions and other STEP exports are preserved under `development/versions/` and `development/geometry/`.

STL2 is used as the meshing tool, then the material-tagged mesh is converted to DAGMC `.h5m` through the standard CAD_to_OpenMC workflow. OpenMC consumes those `.h5m` files. The generated `somuv4.h5m` and comparison `benchmark.h5m` files are not committed. The assembly was created with `CAD_to_OpenMC.assembly` and `Assembly(["somuv4.step"])`, but the original material-assignment and conversion setup commands were not fully preserved. As a result, this repository does not provide a complete line-by-line recipe for rebuilding the exact DAGMC geometry. The design model expects the following DAGMC material tags to match the OpenMC material names exactly:

`air`, `carrier`, `clad`, `d20`, `fuel`, `graphite`, `h20`, `pressure`

The historical RBMK comparison script is archived at `development/versions/rbmkagain/rbmkv2.py`. It uses a hard-coded local geometry path and is not wired into the current workflows. The separate `normalized_void.py` script is not included, so the RBMK comparison is not reproducible end to end from this repository.

The thermal-hydraulics input power map is `thermal_hydraulics/data/power_map.csv`. The ENDF/B-VII.1 library requires `cross_sections.xml` plus the referenced HDF5 data files. The reported library has 444 files under `neutron/`. No checksum or fingerprint for that library is recorded in the files included here.

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

The archived RBMK model is `development/versions/rbmkagain/rbmkv2.py`. It currently evaluates only 0% and 20% void and contains a hard-coded local geometry path. The `normalized_void.py` post-processing script is not included. The repository therefore does not contain a complete RBMK sweep and result-regeneration workflow.

## Thermal-hydraulics workflows

The thermal-hydraulics scripts are in `thermal_hydraulics/src/`; the configuration is in `thermal_hydraulics/config/config.yaml`, and the power map is in `thermal_hydraulics/data/power_map.csv`. In the current checkout, `run_cases.py` and `power_scan.py` expect `config.yaml` and `power_map.csv` alongside the scripts, while the repository stores them in separate folders. These two scripts require manual path adjustments to run from the checkout as organized. `channel.py` accepts the configuration path explicitly, but the power-map path in that configuration must also resolve to the data file.

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
