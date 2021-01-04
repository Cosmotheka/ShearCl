# Cosmic shear power spectra

This repository accompanies the paper "Cosmic shear power spectra in practice" [arXiv:2010.09717](https://arxiv.org/abs/2010.09717), which covers in detail the main steps in the process of estimating pseudo-CL power spectra for cosmic shear, including their covariance matrices, and interpreting them in a cosmological analysis.

The resources hosted here are:

## Data
The paper presents measurements of the shear power spectrum, their covariance matrix and accompanying metadata for the first public data releases from the Hyper-Suprime Cam (HSC) and Dark Energy Survey (DES) collaborations. We make these power spectra publicly available here:
1. **HSC data:**
  * Shear power spectra, Gaussian covariance, bandpower window functions and redshift distributions can be found [here](https://entangled.physics.ox.ac.uk/index.php/s/QAnxOWhL4NZRhkU/download).
  * Noise bias and non-Gaussian contributions to the shear power spectrum covariance matrix can be found [here](https://entangled.physics.ox.ac.uk/index.php/s/CELKhMclym4mBa2/download).

2. **DES data:**
  * Shear power spectra, Gaussian covariance, bandpower window functions and redshift distributions can be found [here](https://entangled.physics.ox.ac.uk/index.php/s/Sx1gzL1kMAEgPDo/download).
  * Noise bias and non-Gaussian contributions to the shear power spectrum covariance matrix can be found [here](https://entangled.physics.ox.ac.uk/index.php/s/btEU2ltEKPvpPwH/download).
  * PSF null tests can be found [here](https://entangled.physics.ox.ac.uk/index.php/s/XQi98FwzDtw4apP/download).

These data are available in [SACC](https://github.com/LSSTDESC/sacc) format. See the [documentation](https://sacc.readthedocs.io/en/latest/) in the SACC repository and the examples provided below for details on how to read and interpret the data.

## A working pipeline
The folder [pipeline_des](pipeline_des) contains all the scripts used to download and analyse the DES data as a practical demonstration of the methods described in the paper.

## Examples
The [examples](examples) folder contains two notebooks describing in some detail the contents of the [HSC](examples/ClExampleHSC.ipynb) and [DES](examples/ClExampleDES.ipynb) data files and how to use them.

## Credit and feedback
If you use these data or the methods presented in the [paper](https://arxiv.org/abs/2010.09717), we kindly ask you to cite it.

If you have questions or find any issues with the data, feel free to open an issue here or contact the authors directly:

    Andrina Nicola (@anicola)
    Carlos Garcia-Garcia (@carlosggarcia)
    David Alonso (@damonge)
