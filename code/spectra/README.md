# Spectra of plain DGSEM and shock capturing semidiscretizations

This folder provides setups and results to compute spectra of plain DG and shock
capturing semidiscretizations in Trixi.jl and their corresponding linear CFL
factors for a range of RK methods from OrdinaryDiffEq.jl. To reproduce the
results, you need to execute

```bash
julia -e 'include("spectra.jl"); run_experiments()'
```
