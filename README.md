# Natural-Uranium Reactor Research

Computational research artefacts for a conceptual natural-uranium reactor study. The repository contains the two final workflows that emerged from an iterative development process:

- **Neutronics** — OpenMC/DAGMC-based analyses of a baseline configuration and selected moderator-drain, loss-of-coolant, and void cases.
- **Thermal-hydraulics** — a one-dimensional pressure-channel model using IAPWS-97 properties and an axial power distribution derived from the neutronics analysis.

## Author

Adarsh Hazra — 10th-grade student and nuclear-engineering enthusiast.

## Scope and status

This is an exploratory computational research repository, not a validated reactor design, safety case, operating manual, or engineering specification. Results depend on modelling assumptions and should not be used for design, licensing, construction, or operation.

Raw CAD, mesh, and OpenMC statepoint files are intentionally not included in this first repository version. That keeps the public-facing material focused on methodology, code, compact derived data, and figures while the geometry and detailed design information undergo a separate publication review.

## Repository layout

```text
neutronics/            OpenMC model source, workflows, compact results, and figures
thermal_hydraulics/    Channel-model source, configuration, results, and figures
docs/                  Methodology, assumptions, publication scope, and development history
requirements/          Python package lists for the two workflows
archive/               Curated record of the development process (not raw trial directories)
```

## Workflow relationship

The neutronics workflow produces axial power information. `neutronics/data/axial_power_profile.csv` and `thermal_hydraulics/data/power_map.csv` preserve the compact data interface used by the thermal-hydraulics model.

## Reproduction status

The thermal-hydraulics workflow is represented by its source, configuration, power map, compact results, and figures. The neutronics source is included for transparency, but full reruns require the withheld DAGMC geometry file and the appropriate OpenMC/DAGMC data environment. Supply the geometry path through `SOMUV4_H5M`; the fallback location is `private/somuv4.h5m`, which is ignored by Git. See [publication scope](docs/publication-scope.md) and [assumptions and limitations](docs/assumptions-and-limitations.md).

## Before publication

1. Review every assumption and reported result against the underlying research record.
2. Obtain appropriate review before making geometry, detailed design inputs, or raw statepoints public.

