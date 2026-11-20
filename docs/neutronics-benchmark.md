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

This is a leakage sensitivity estimate based on an infinite-lattice calculation, not a benchmark against an experimental reactor or a validated single-module transport calculation. The 3% leakage assumption is not established by the reflected block model. The result is sensitive to geometry, material data, boundary conditions, nuclear data, and the leakage model. A finite bare-module OpenMC model with explicit leakage would be needed to answer whether that module is critical without relying on the assumed correction.

The report’s cited textbooks and code papers provide background for reactor analysis and computational methods, but the report does not identify a specific validated benchmark dataset or derive the 3% leakage value from a source. Record such a source and a reproducible comparison case here if one is added later.
