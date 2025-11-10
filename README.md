# Asymptotic-preserving and energy-conserving methods for a hyperbolic approximation of the BBM equation

[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17566794.svg)](https://doi.org/10.5281/zenodo.17566794)

This repository contains information and code to reproduce the results presented in the article

TODO: Insert ArXiv reference

If you find these results useful, please cite the article mentioned above. If you use the implementations provided here, please also cite this repository as

```bibtex
@misc{bleecke2025asymptoticRepro,
  title={Reproducibility repository for
         "{A}symptotic-preserving and energy-conserving methods for a
         hyperbolic approximation of the BBM equation"},
  author={Bleecke, Sebastian and Biswas, Abhijit and Ketcheson, David I.
          and Ranocha, Hendrik and Sch{\"u}tz, Jochen},
  year={2025},
  doi={10.5281/zenodo.17566794},
  howpublished={\url{https://github.com/sbleecke/2025_bbmh}},
}
```

## Abstract
We study the hyperbolic approximation of the Benjamin-Bona-Mahony (BBM) equation
proposed recently by Gavrilyuk and Shyue (2022).
We develop asymptotic-preserving numerical methods using implicit-explicit
(additive) Runge-Kutta methods that are implicit in the stiff linear part.
The new discretization of the hyperbolization conserves important invariants
converging to invariants of the BBM equation.
We use the entropy relaxation approach to make the fully discrete schemes
energy-preserving.
Numerical experiments demonstrate the effectiveness of these discretizations.

## Numerical experiments
To reproduce the numerical experiments presented in this article, you need to install Julia and Python.
The numerical experiments presented in this article were performed using Julia v1.10.5 and Python 3.13.7.

First, you need to download this repository, e.g., by cloning it with git or by downloading an archive via the GitHub interface. Then, you need to start Julia in the code directory of this repository and follow the instructions described in the README.md file therein.

## Authors
* Sebastian Bleecke (Johannes Gutenberg University Mainz, Germany)
* Abhijit Biswas (Indian Institute of Technology Kanpur, India)
* David I. Ketcheson (KAUST, Saudi Arabia)
* Hendrik Ranocha (Johannes Gutenberg University Mainz, Germany)
* Jochen Schütz (Hasselt University, Belgium)

## License
The code in this repository is published under the MIT license, see the `LICENSE` file.

## Disclaimer
Everything is provided as is and without warranty. Use at your own risk!

