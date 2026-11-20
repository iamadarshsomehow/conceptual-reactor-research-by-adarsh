# Neutronics criticality estimate

## What the calculation shows

The report gives a reflected single-block multiplication factor of `k-infinity = 1.05892 ± 0.00053`. This result is above 1 for the modeled infinite lattice. It does not by itself establish that a finite, unreflected module or full reactor is critical.

The report estimates finite-system criticality by applying an assumed leakage fraction `L`:

```text
k-effective ≈ k-infinity × (1 − L)
```

With the assumed `L = 0.03`, the estimate is `1.05892 × 0.97 ≈ 1.027`.

Under this simplified correction, the break-even leakage fraction is `1 − 1/k-infinity ≈ 0.0557`, or about 5.6%. Thus, in this estimate, 3% leakage gives a value above 1; leakage greater than 3% does not automatically make the estimate subcritical. It remains above 1 only while the assumed leakage is below roughly 5.6%.

## Benchmark status

The project's benchmark check was the 3% leakage check alone: apply the assumed 3% leakage correction to the bare, reflected single-module `k-infinity` result and check whether the resulting estimate is above 1. By this project criterion, `k-effective ≈ 1.027`, so it passes. This was the project's chosen screening check, not an external code-validation benchmark or a full-core calculation.

The 3% leakage assumption is not established by the reflected block model. The result is sensitive to geometry, material data, boundary conditions, nuclear data, and the leakage model. A finite bare-module OpenMC model with explicit leakage would test criticality without relying on this correction, but it was not part of the benchmark check used for this project.

The report’s cited textbooks and code papers provide background for reactor analysis and computational methods, but the report does not identify a specific validated benchmark dataset or derive the 3% leakage value from a source. Record such a source and a reproducible comparison case here if one is added later.
