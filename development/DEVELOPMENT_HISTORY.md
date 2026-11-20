# Development History

## Early learning projects

The `rbmkagain reactor test` and `rbmk rebuild` projects were early exercises in learning how reactor geometry and simulation workflows could be developed.

## Exploring the neutronics approach

During versions v1 through v4, I explored different approaches to determine which would work best for this conceptual study. I settled on simulating a single reactor module and applying an assumed 3% neutron-leakage correction to its calculated `k-infinity` to assess criticality. This became the project's chosen neutronics screening method.

## Design iterations

From v5 through v12, I made many design improvements. After analyzing those versions, I selected v10 as the strongest starting point. The later improvements involved changes to materials and operating conditions, as well as design features such as sizes, shapes, and volumes.

I then developed further versions from v10, including v10plus, somum, somux, v2, and v3. These iterations led to somuv4, which passed the neutronics benchmark criterion established for this project.

## Thermal-hydraulics work

After somuv4 passed the project's neutronics benchmark criterion, I developed its thermal-hydraulics analysis. That work was kept together in a single directory named `thermal`.
