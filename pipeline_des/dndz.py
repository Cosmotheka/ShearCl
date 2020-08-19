import numpy as np
from astropy.io import fits
import utils as ut


predir = ut.rootpath + '/'
fname_nz = predir + "DES_data/shear_catalog/y1_redshift_distributions_v1.fits"
d = fits.open(fname_nz)[1].data
for i in range(4):
    fname = predir + 'outputs/dndz_metacal_bin%d.txt' %i
    np.savetxt(fname,
               np.transpose([d['Z_LOW'], d['Z_MID'], d['Z_HIGH'],
                             d['BIN%d' % (i+1)]]))

