# On error-based step size control for discontinuous Galerkin methods for compressible fluid dynamics

[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7078946.svg)](https://doi.org/10.5281/zenodo.7078946)

This repository contains information and code to reproduce the results presented in the
article
```bibtex
@article{ranocha2023error,
  title={On error-based step size control for discontinuous {G}alerkin methods
         for compressible fluid dynamics},
  author={Ranocha, Hendrik and Winters, Andrew R and Castro, Hugo Guillermo
          and Dalcin, Lisandro and Schlottke-Lakemper, Michael and
          Gassner, Gregor J and Parsani, Matteo},
  journal={Communications on Applied Mathematics and Computation},
  year={2023},
  month={05},
  doi={10.1007/s42967-023-00264-y},
  eprint={2209.07037},
  eprinttype={arxiv},
  eprintclass={math.NA}
}
```

If you find these results useful, please cite the article mentioned above. If you
use the implementations provided here, please **also** cite this repository as
```bibtex
@misc{ranocha2022errorRepro,
  title={Reproducibility repository for
         "{O}n error-based step size control for discontinuous {G}alerkin methods
         for compressible fluid dynamics"},
  author={Ranocha, Hendrik and Winters, Andrew R and Castro, Hugo Guillermo
          and Dalcin, Lisandro and Schlottke-Lakemper, Michael and
          Gassner, Gregor J and Parsani, Matteo},
  year={2022},
  howpublished={\url{https://github.com/trixi-framework/paper-2022-stepsize_control}},
  doi={10.5281/zenodo.7078946}
}
```


## Abstract

We study temporal step size control of explicit Runge-Kutta methods for
compressible computational fluid dynamics (CFD), including the Navier-Stokes
equations and hyperbolic systems of conservation laws such as the Euler equations.
We demonstrate that error-based approaches
are convenient in a wide range of applications and compare them to more classical
step size control based on a Courant-Friedrichs-Lewy (CFL) number. Our numerical
examples show that error-based step size control is easy to use, robust, and efficient,
e.g., for (initial) transient periods, complex geometries, nonlinear shock
capturing approaches, and schemes that use nonlinear entropy projections.
We demonstrate these properties for problems ranging from well-understood
academic test cases to industrially relevant large-scale computations with two
disjoint code bases, the open source Julia packages Trixi.jl with OrdinaryDiffEq.jl
and the C/Fortran code SSDC based on PETSc.

## Numerical experiments

The numerical experiments presented in the paper use [Trixi.jl](https://github.com/trixi-framework/Trixi.jl)
and SSDC.
To reproduce the numerical experiments using Trixi.jl, you need to install
[Julia](https://julialang.org/).

The subfolders of this repository contain `README.md` files with instructions
to reproduce the numerical experiments, including postprocessing.

The numerical experiments were carried out using Julia v1.7.


## Authors

- [Hendrik Ranocha](https://ranocha.de) (University of Hamburg, Germany)
- [Andrew Winters](https://liu.se/en/employee/andwi94) (Linköping University, Sweden)
- Hugo Guillermo Castro (KAUST, Saudi Arabia)
- Lisandro Dalcin (KAUST, Saudi Arabia)
- Michael Schlottke-Lakemper (University of Stuttgart, Germany)
- Gregor J. Gassner (University of Cologne, Germany)
- Matteo Parsani (KAUST, Saudi Arabia)


## Disclaimer

Everything is provided as is and without warranty. Use at your own risk!
