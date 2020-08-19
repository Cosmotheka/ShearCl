import numpy as np
import healpy as hp
from argparse import ArgumentParser
import pymaster as nmt
import utils as ut
import pyccl as ccl


parser = ArgumentParser()
parser.add_argument("--bin-number", default=0, type=int, help="Bin number")
parser.add_argument("--nside", default=4096, type=int, help="Nside")
o = parser.parse_args()


predir = ut.rootpath = '/outputs/'

npix = hp.nside2npix(o.nside)

ut.printflush("Running %d-PSF" % (o.bin_number))

ut.printflush("MCM")
predir_mcm = predir + 'cls_metacal_mcm_bins_'
fname_mcm = predir_mcm + '%d%d_ns%d.fits' % (o.bin_number,
                                             o.bin_number,
                                             o.nside)
w = nmt.NmtWorkspace()
w.read_from(fname_mcm)

ut.printflush("Theory spectra")
ls = np.arange(3*o.nside)
cl0 = np.zeros(3*o.nside)
# Signal
cosmo = ccl.Cosmology(Omega_c=0.260-0.0479,
                      Omega_b=0.0479,
                      h=0.685,
                      n_s=0.973,
                      sigma8=0.821)
fname_nz = predir + "dndz_metacal_bin%d.txt" % o.bin_number
zi, zm, zf, dndz = np.loadtxt(fname_nz, unpack=True)
tr = ccl.WeakLensingTracer(cosmo, (zi, dndz))
sl = ccl.angular_cl(cosmo, tr, tr, ls)
# Noise
nl = np.zeros(3*o.nside)
fname_nl = predir + "maps_metacal_bin%d_ns%d_nells.npz" % (o.bin_number,
                                                           o.nside)
d = np.load(fname_nl)
nl[2:] = d['nl_cov']
# G-G
prefix = predir + "maps_metacal_bin%d_ns%d" % (o.bin_number, o.nside)
msk = np.load(prefix + '_w.npz')['w']
fsky = np.sum(msk**2) / npix
cl_gg = w.couple_cell([sl, cl0, cl0, cl0])/fsky + np.array([nl, cl0, cl0, nl])
# PSF-PSF
fname_apsf = predir + "cls_metacal_cls_bins_"
fname_apsf += "%d%d_ns%d_apsf.npz" % (o.bin_number,
                                      o.bin_number,
                                      o.nside)
dpsf = np.load(fname_apsf)
cl_pp = dpsf['cls_coupled'] / fsky
# G-PSF
cl_gp = np.array([cl0, cl0, cl0, cl0])

ut.printflush("CMCM")
fname_cmcm = predir + 'cls_metacal_cmcm_bins_'
fname_cmcm += '%d%d_%d%d_ns%d.fits' % (o.bin_number, o.bin_number,
                                       o.bin_number, o.bin_number,
                                       o.nside)
cw = nmt.NmtCovarianceWorkspace()
cw.read_from(fname_cmcm)

ut.printflush("Covariance")
nbpw = w.wsp.bin.n_bands
cov = nmt.gaussian_covariance(cw, 2, 2, 2, 2,
                              cl_gg, cl_gp, cl_gp, cl_pp,
                              w, w).reshape([nbpw, 4, nbpw, 4])

ut.printflush("Writing")
fname_cov = predir + 'cls_metacal_covar_xpsf_bin'
fname_cov += '%d_ns%d.npz' % (o.bin_number, o.nside)
np.savez(fname_cov, cov=cov)
