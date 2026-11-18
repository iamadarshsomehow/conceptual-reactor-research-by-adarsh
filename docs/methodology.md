# Methodology

## Neutronics

The final neutronics workflow is represented by `neutronics/src/model.py`, `neutronics/workflows/sweep_all.py`, and `neutronics/workflows/drain_runs.py`. It uses OpenMC with a DAGMC geometry input and evaluates a baseline plus selected void, moderator-drain, and loss-of-coolant cases. Compact CSV tables and figures are retained in `neutronics/results/`.

## Thermal-hydraulics

The final thermal-hydraulics workflow is represented by `thermal_hydraulics/src/`. It models a one-dimensional pressure channel, uses IAPWS-97 properties, accepts an axial power map, and evaluates nominal and sensitivity cases. `thermal_hydraulics/config/config.yaml` is the central record of model inputs and assumptions.

## Coupling boundary

The workflows are coupled through an axial power distribution. This repository preserves compact power-profile CSV data but does not publish raw OpenMC statepoints or geometry.
