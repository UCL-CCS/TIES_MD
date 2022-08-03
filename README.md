[![GH Actions Status](https://github.com/UCL-CCS/TIES_MD/workflows/CI/badge.svg)](https://github.com/UCL-CCS/TIES_MD/actions?query=branch%3Amaster+workflow%3ACI)
[![Conda](https://anaconda.org/conda-forge/ties_md/badges/version.svg)](https://anaconda.org/conda-forge/ties_md)
[![Anaconda Cloud Badge](https://anaconda.org/conda-forge/ties_md/badges/downloads.svg)](https://anaconda.org/conda-forge/ties_md)

<img src="https://github.com/UCL-CCS/TIES_MD/blob/main/TIES_MD/doc/source/_static/images/TIES_logov2.png" width="384">

Please see our [documentation](https://UCL-CCS.github.io/TIES_MD/) for more information on `TIES MD` or our main [TIES](http://www.ties-service.org/) page for information on all TIES packages.

To install `TIES_MD` run:

`conda install -c conda-forge ties_md`

`TIES` is a collection of software packages which can be used to calculate protein ligand binding free energies with physics based alchemical methods. `TIES` stands for thermodynamic integration with enhanced sampling and this method focuses on the use of ensemble simulations to control for the alleotoric errors in molecular dynamics simulations.

Within the `TIES` collection there are two packages: `TIES 20` and `TIES MD`. `TIES 20` can be used to build all input for `TIES MD`. `TIES MD` is an implementation of the `TIES` protocol that is built on `OpenMM` and `NAMD`.