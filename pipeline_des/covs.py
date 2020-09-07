import numpy as np
import healpy as hp
from argparse import ArgumentParser
import pymaster as nmt
import utils as ut
import os
import pyccl as ccl


parser = ArgumentParser()
parser.add_argument("--bin-a1", default=0, type=int, help="Bin number")
parser.add_argument("--bin-a2", default=0, type=int, help="Bin number")
parser.add_argument("--bin-b1", default=0, type=int, help="Bin number")
parser.add_argument("--bin-b2", default=0, type=int, help="Bin number")
parser.add_argument("--nside", default=4096, type=int, help="Nside")
parser.add_argument("--n-iter", default=0, type=int, help="n_iter")
parser.add_argument("--recompute-mcm", default=False, action='store_true',
                    help="Recompute MCM even if it exists?")
parser.add_argument("--old-nka", default=False, action='store_true',
                    help="Use old NKA")
parser.add_argument("--full-noise", default=False, action='store_true',
                    help="Full noise calculation?")
o = parser.parse_args()

predir = ut.rootpath + '/outputs/'

npix = hp.nside2npix(o.nside)
pix_area = 4*np.pi/npix

ut.printflush("Running %d-%d %d-%d" % (o.bin_a1, o.bin_a2, o.bin_b1, o.bin_b2))

ut.printflush("Theory spectra")
ls = np.arange(3*o.nside)
cl0 = np.zeros(3*o.nside)
cl1 = np.ones(3*o.nside)

cosmo = ccl.Cosmology(Omega_c=0.264,
                      Omega_b=0.0493,
                      h=0.6736,
                      n_s=0.9649,
                      sigma8=0.8111)


def get_field(bin_no, return_mask=False, mask_sigma=False):
    prefix = predir + "maps_metacal_bin%d_ns%d" % (bin_no, o.nside)
    pix = np.load(prefix + '_goodpix.npz')['pix']
    wmap = np.zeros(npix)
    wmap[pix] = np.load(prefix + '_w.npz')['w']
    if mask_sigma:
        wmap[pix] = np.sqrt(np.load(prefix + '_w2s2.npz')['w2s2'])
    else:
        wmap[pix] = np.load(prefix + '_w.npz')['w']
    if return_mask:
        return wmap
    return nmt.NmtField(wmap, [wmap, wmap], n_iter=o.n_iter)


def get_tracer(bin_no):
    fname_nz = predir + "dndz_metacal_bin%d.txt" % bin_no
    zi, zm, zf, dndz = np.loadtxt(fname_nz, unpack=True)
    tr = ccl.WeakLensingTracer(cosmo, (zi, dndz))
    return tr


tracers = {}
if o.bin_a1 not in tracers:
    tracers[o.bin_a1] = get_tracer(o.bin_a1)
if o.bin_a2 not in tracers:
    tracers[o.bin_a2] = get_tracer(o.bin_a2)
if o.bin_b1 not in tracers:
    tracers[o.bin_b1] = get_tracer(o.bin_b1)
if o.bin_b2 not in tracers:
    tracers[o.bin_b2] = get_tracer(o.bin_b2)


def get_cl(trs, b1, b2):
    nl = np.zeros(3*o.nside)
    if (b1 == b2) and (not o.full_noise):
        fname_nl = predir + "maps_metacal_bin%d_ns%d_nells.npz" % (b1, o.nside)
        d = np.load(fname_nl)
        nl[2:] = d['nl_cov']

    sl = ccl.angular_cl(cosmo, trs[b1], trs[b2], ls)
    if o.old_nka:
        return np.array([sl + nl, cl0, cl0, nl])
    else:
        w = nmt.NmtWorkspace()
        predir_mcm = predir + 'cls_metacal_mcm_bins_'
        fname_mcm = predir_mcm + '%d%d_ns%d.fits' % (b1, b2, o.nside)
        if os.path.isfile(fname_mcm):
            w.read_from(fname_mcm)
        else:
            fname_mcm = predir_mcm + '%d%d_ns%d.fits' % (b2, b1, o.nside)
            if os.path.isfile(fname_mcm):
                w.read_from(fname_mcm)
            else:
                raise ValueError("Can't find MCM " + fname_mcm)
        mskprod = get_field(b1, return_mask=True)
        if b1 == b2:
            mskprod *= mskprod
        else:
            mskprod *= get_field(b2, return_mask=True)
        fsky = np.mean(mskprod)
        return w.couple_cell([sl, cl0, cl0, cl0])/fsky + \
            np.array([nl, cl0, cl0, nl])


clt = {}
k = '%d%d' % (o.bin_a1, o.bin_b1)
if k not in clt:
    clt[k] = get_cl(tracers, o.bin_a1, o.bin_b1)
k = '%d%d' % (o.bin_a1, o.bin_b2)
if k not in clt:
    clt[k] = get_cl(tracers, o.bin_a1, o.bin_b2)
k = '%d%d' % (o.bin_a2, o.bin_b1)
if k not in clt:
    clt[k] = get_cl(tracers, o.bin_a2, o.bin_b1)
k = '%d%d' % (o.bin_a2, o.bin_b2)
if k not in clt:
    clt[k] = get_cl(tracers, o.bin_a2, o.bin_b2)


ut.printflush("CMCM")
fname_cmcm = predir + 'cls_metacal_cmcm_bins_'
fname_cmcm += '%d%d_%d%d_ns%d.fits' % (o.bin_a1, o.bin_a2, o.bin_b1,
                                       o.bin_b2, o.nside)

cw = nmt.NmtCovarianceWorkspace()
fields_s = {}
fields_n = {}
if os.path.isfile(fname_cmcm) and not o.recompute_mcm:
    ut.printflush(" - Reading")
    cw.read_from(fname_cmcm)
else:
    ut.printflush(" - Fields")
    if o.bin_a1 not in fields_s:
        fields_s[o.bin_a1] = get_field(o.bin_a1)
    if o.bin_a2 not in fields_s:
        fields_s[o.bin_a2] = get_field(o.bin_a2)
    if o.bin_b1 not in fields_s:
        fields_s[o.bin_b1] = get_field(o.bin_b1)
    if o.bin_b2 not in fields_s:
        fields_s[o.bin_b2] = get_field(o.bin_b2)
    ut.printflush(" - Computing")
    cw.compute_coupling_coefficients(fields_s[o.bin_a1], fields_s[o.bin_a2],
                                     fields_s[o.bin_b1], fields_s[o.bin_b2])
    cw.write_to(fname_cmcm)
if o.full_noise:
    ut.printflush("CMCM - SN")
    prefix_cmcm = predir + 'cls_metacal_cmcm_sn_bins_'
    cw_sn = {}

    def get_cw_sn(b_a1, b_a2, b_b1, b_b2,
                  n_a1, n_a2, n_b1, n_b2):
        name = "%d%d_%d%d" % (b_a1, b_a2, b_b1, b_b2)
        name += "_%d%d_%d%d" % (n_a1, n_a2, n_b1, n_b2)
        fields_a = [fields_s, fields_n]

        if name in cw_sn:
            pass
        else:
            cw_sn[name] = nmt.NmtCovarianceWorkspace()
            fname_cmcm = prefix_cmcm + name + '_ns%d.fits' % o.nside
            if os.path.isfile(fname_cmcm):
                ut.printflush(" - Reading")
                cw_sn[name].read_from(fname_cmcm)
            else:
                ut.printflush(" - Fields")
                if b_a1 not in fields_a[n_a1]:
                    fields_a[n_a1][b_a1] = get_field(b_a1,
                                                     mask_sigma=bool(n_a1))
                if b_a2 not in fields_a[n_a2]:
                    fields_a[n_a2][b_a2] = get_field(b_a2,
                                                     mask_sigma=bool(n_a2))
                if b_b1 not in fields_a[n_b1]:
                    fields_a[n_b1][b_b1] = get_field(b_b1,
                                                     mask_sigma=bool(n_b1))
                if b_b2 not in fields_a[n_b2]:
                    fields_a[n_b2][b_b2] = get_field(b_b2,
                                                     mask_sigma=bool(n_b2))
                ut.printflush(" - Computing " + name)
                cw_sn[name].compute_coupling_coefficients(fields_a[n_a1][b_a1],
                                                          fields_a[n_a2][b_a2],
                                                          fields_a[n_b1][b_b1],
                                                          fields_a[n_b2][b_b2])
                cw_sn[name].write_to(fname_cmcm)

    if o.bin_a1 == o.bin_b1:
        get_cw_sn(o.bin_a1, o.bin_a2, o.bin_b1, o.bin_b2, 1, 0, 1, 0)
    if o.bin_a1 == o.bin_b2:
        get_cw_sn(o.bin_a1, o.bin_a2, o.bin_b1, o.bin_b2, 1, 0, 0, 1)
    if o.bin_a2 == o.bin_b1:
        get_cw_sn(o.bin_a1, o.bin_a2, o.bin_b1, o.bin_b2, 0, 1, 1, 0)
    if o.bin_a2 == o.bin_b2:
        get_cw_sn(o.bin_a1, o.bin_a2, o.bin_b1, o.bin_b2, 0, 1, 0, 1)

    ut.printflush("CMCM - NN")
    prefix_cmcm = predir + 'cls_metacal_cmcm_nn_bins_'
    cw_nn = {}

    def get_cw_nn(b_a1, b_a2, b_b1, b_b2):
        name = "%d%d_%d%d" % (b_a1, b_a2, b_b1, b_b2)
        if name in cw_nn:
            pass
        else:
            cw_nn[name] = nmt.NmtCovarianceWorkspace()
            fname_cmcm = prefix_cmcm + name + '_ns%d.fits' % o.nside
            if os.path.isfile(fname_cmcm):
                ut.printflush(" - Reading")
                cw_nn[name].read_from(fname_cmcm)
            else:
                ut.printflush(" - Fields")
                for b in [b_a1, b_a2, b_b1, b_b2]:
                    fields_n[b] = get_field(b, mask_sigma=True)
                ut.printflush(" - Computing " + name)
                cw_nn[name].compute_coupling_coefficients(fields_n[b_a1],
                                                          fields_n[b_a2],
                                                          fields_n[b_b1],
                                                          fields_n[b_b2])
                cw_nn[name].write_to(fname_cmcm)

    if (((o.bin_a1 == o.bin_b1) and (o.bin_a2 == o.bin_b2)) or
        ((o.bin_a1 == o.bin_b2) and (o.bin_a2 == o.bin_b1))):
        get_cw_nn(o.bin_a1, o.bin_a2, o.bin_b1, o.bin_b2)


ut.printflush("MCMs")
fname_mcm_a = predir + 'cls_metacal_mcm_bins_'
fname_mcm_a += '%d%d_ns%d.fits' % (o.bin_a1, o.bin_a2, o.nside)
wa = nmt.NmtWorkspace()
wa.read_from(fname_mcm_a)
if (o.bin_a1 == o.bin_b1) and (o.bin_a2 == o.bin_b2):
    wb = wa
else:
    fname_mcm_b = predir + 'cls_metacal_mcm_bins_'
    fname_mcm_b += '%d%d_ns%d.fits' % (o.bin_b1, o.bin_b2, o.nside)
    wb = nmt.NmtWorkspace()
    wb.read_from(fname_mcm_b)
nbpw = wa.wsp.bin.n_bands


ut.printflush("Covariance")
cshape = [nbpw, 4, nbpw, 4]
cov = nmt.gaussian_covariance(cw, 2, 2, 2, 2,
                              clt['%d%d' % (o.bin_a1, o.bin_b1)],
                              clt['%d%d' % (o.bin_a1, o.bin_b2)],
                              clt['%d%d' % (o.bin_a2, o.bin_b1)],
                              clt['%d%d' % (o.bin_a2, o.bin_b2)],
                              wa, wb).reshape(cshape)
if o.full_noise:
    name = '%d%d_%d%d' % (o.bin_a1, o.bin_a2, o.bin_b1, o.bin_b2)
    cl_ones = np.array([cl1, cl0, cl0, cl1])
    cl_zeros = np.array([cl0, cl0, cl0, cl0])
    if o.bin_a1 == o.bin_b1:
        cov += nmt.gaussian_covariance(cw_sn[name + '_10_10'], 2, 2, 2, 2,
                                       cl_ones,
                                       cl_zeros,
                                       cl_zeros,
                                       clt['%d%d' % (o.bin_a2, o.bin_b2)],
                                       wa, wb).reshape(cshape)*pix_area
        if o.bin_a2 == o.bin_b2:
            cov += nmt.gaussian_covariance(cw_nn[name], 2, 2, 2, 2,
                                           cl_ones, cl_zeros,
                                           cl_zeros, cl_ones,
                                           wa, wb).reshape(cshape)*pix_area**2
    if o.bin_a1 == o.bin_b2:
        cov += nmt.gaussian_covariance(cw_sn[name + '_10_01'], 2, 2, 2, 2,
                                       cl_zeros,
                                       cl_ones,
                                       clt['%d%d' % (o.bin_a2, o.bin_b1)],
                                       cl_zeros,
                                       wa, wb).reshape(cshape)*pix_area
        if o.bin_a2 == o.bin_b1:
            cov += nmt.gaussian_covariance(cw_nn[name], 2, 2, 2, 2,
                                           cl_zeros, cl_ones,
                                           cl_ones, cl_zeros,
                                           wa, wb).reshape(cshape)*pix_area**2
    if o.bin_a2 == o.bin_b1:
        cov += nmt.gaussian_covariance(cw_sn[name + '_01_10'], 2, 2, 2, 2,
                                       cl_zeros,
                                       clt['%d%d' % (o.bin_a1, o.bin_b2)],
                                       cl_ones,
                                       cl_zeros,
                                       wa, wb).reshape(cshape)*pix_area
    if o.bin_a2 == o.bin_b2:
        cov += nmt.gaussian_covariance(cw_sn[name + '_01_01'], 2, 2, 2, 2,
                                       clt['%d%d' % (o.bin_a1, o.bin_b1)],
                                       cl_zeros,
                                       cl_zeros,
                                       cl_ones,
                                       wa, wb).reshape(cshape)*pix_area

ut.printflush("Writing")
fname_cov = predir + 'cls_metacal_covar_bins_'
if not o.old_nka:
    fname_cov += "new_nka_"
if o.full_noise:
    fname_cov += "full_noise_"
fname_cov += '%d%d_%d%d_ns%d.npz' % (o.bin_a1, o.bin_a2,
                                     o.bin_b1, o.bin_b2, o.nside)
np.savez(fname_cov, cov=cov)
