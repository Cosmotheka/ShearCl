# DES shear power spectrum pipeline

This folder contains all the code needed to reproduce the shear power spectra calculated in the paper (including their covariance and null tests).

The different scripts are as follows:
- `utils.py` contains a useful printing function and the root path used by all other scripts.
- `download.py` downloads the shear data from the Y1 data release website.
- `subsample.py` joins and splits the data into separate smaller catalogs for each redshift bin, so each catalog can be read in and operated on quickly.
- `map.py` creates all necessary maps from the catalog (shear maps, rotated maps, PSF maps, weight maps, ellipticity variance maps...).
- `cls.py` computes the pseudo-Cl for a given pair of redshift bins.
- `covs.py` computes the Gaussian covariance matrix element for the power spectra of two pairs of redshift bins.
- `covs_xPSF.py` does the same as `covs.py` for power spectra involving cross-correlations with the PSF maps.
- `windows.py` computes and stores the bandpower window functions for the power spectrum between two redshift bins.
- `pipeline.sh` executes all the scripts above in the right order.
