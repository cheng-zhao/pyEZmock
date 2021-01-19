# pyEZmock

![GitHub](https://img.shields.io/github/license/cheng-zhao/pyEZmock.svg)

## Introduction

This is a python wrapper of the *Effective Zel'dovich approximation mock generator*<sup>[\[1\]](#ref1)</sup> (EZmock), as well as the clustering measurement codes [powspec](https://github.com/cheng-zhao/powspec), [FCFC](https://github.com/cheng-zhao/FCFC), and bispec<sup>[*](#note1)</sup>.

It creates bash scripts for running EZmocks and measuring clustering statistics, and submits them to the [Slurm workload manager](https://slurm.schedmd.com/documentation.html) if necessary. It is designed for the *National Energy Research Scientific Computing Center* ([NERSC](https://www.nersc.gov/)), but can be in principle used in any computing environment, provided the dependencies are given.

This package is written by Cheng Zhao (赵成), and is distributed under the [MIT license](LICENSE.txt).

If you use this tool in research that results in publications, please cite the following papers:

> Chuang et al., 2015, [MNRAS, 446, 2621](https://ui.adsabs.harvard.edu/abs/2015MNRAS.446.2621C/abstract)<br />
> Zhao et al. 2020, [arXiv:2007.08997](https://ui.adsabs.harvard.edu/abs/2020arXiv200708997Z/abstract)

<span id="note1">*: to be appeared on github soon.</span>

## Dependencies

The dependencies of this package are as follows:

-   [Python3](https://www.python.org/)  (>= 3.6)
-   [NumPy](https://numpy.org/) \[*optional*\]
-   [Matplotlib](https://matplotlib.org/) \[*optional*\]
-   EZmock<sup>[**](#note2)</sup>
-   [powspec](https://github.com/cheng-zhao/powspec) \[*optional*\]
-   [FCFC](https://github.com/cheng-zhao/FCFC) \[*optional*\]
-   bispec<sup>[***](#note2)</sup> \[*optional*\]

They are all available on NERSC, and linked to this package by default.

<span id="note2">**: can be obtained on request to Chia-Hsun Chuang (chuangch [at] <span>stanford</span>.edu).</span>

<span id="note3">***: can be obtained on request to [Cheng Zhao](mailto:zhaocheng03@gmail.com).</span>

## Usage

Please refer to the example [notebook](example_pyez.ipynb).

## Acknowledgements

This package is inspired by [Tony Zhang's work](https://github.com/neutrinonerd3333/ezmock), and I thank Dr. Chia-Hsun Chuang for his supports and suggestions.

## References

<span id="ref1">\[1\]</span> Chuang, Kitaura, Prada, Zhao, Yepes, 2015, [EZmocks: extending the Zel'dovich approximation to generate mock galaxy catalogues with accurate clustering statistics](https://doi.org/10.1093/mnras/stu2301), *MNRAS*, 446(3):2621&ndash;2628

<span id="ref2">\[2\]</span> Zhao, et al. (*eBOSS Collaboration*), 2020, [The Completed SDSS-IV extended Baryon Oscillation Spectroscopic Survey: one thousand multi-tracer mock catalogues with redshift evolution and systematics for galaxies and quasars of the final data release](https://arxiv.org/abs/2007.08997), *arXiv e-prints*, 2007.08997
