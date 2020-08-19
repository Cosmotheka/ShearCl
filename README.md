# Cosmic shear power spectra

This repository accompanies the paper "Cosmic shear power spectra in practice" [arXiv:2009.ABCDE](dead_link), which covers in details the main steps in the process of estimating pseudo-CL power spectra for cosmic shear, including their covariance matrices, and interpreting them in a cosmological analysis.

The resources hosted here are:

## Data
The paper presents measurements of the shear power spectrum, their covariance matrix and accompanying metadata for the first public data releases from the Hyper-Suprime Cam (HSC) and Dark Energy Survey (DES) collaborations. We make these power spectra publicly available here:
1. **HSC data:**
  * Shear power spectra, Gaussian covariance, bandpower window functions and redshift distributions can be found [here](dead_link).
  * Noise bias and non-Gaussian contributions to the shear power spectrum covariance matrix can be found [here](dead_link).

2. **DES data:**
  * Shear power spectra, Gaussian covariance, bandpower window functions and redshift distributions can be found [here](dead_link).
  * Noise bias and non-Gaussian contributions to the shear power spectrum covariance matrix can be found [here](dead_link).
  * PSF null tests can be found [here](dead_link).

These data are available in [SACC](https://github.com/LSSTDESC/sacc) format. See the [documentation](https://sacc.readthedocs.io/en/latest/) in the SACC repository and the examples provided below for details on how to read and interpret the data.

## A working pipeline
The folder [pipeline_des](pipeline_des) contains all the scripts used to download and analyse the DES data as a practical demonstration of the methods described in the paper.

## Examples

