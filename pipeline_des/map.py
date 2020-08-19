import numpy as np
from astropy.io import fits
from argparse import ArgumentParser
import healpy as hp
import utils as ut


parser = ArgumentParser()
parser.add_argument("--bin-number", default=0, type=int, help="Bin number")
parser.add_argument("--nside", default=4096, type=int, help="Nside")
parser.add_argument("--nrot", default=10, type=int, help="# rotations")
o = parser.parse_args()

npix = hp.nside2npix(o.nside)
pix_area = 4 * np.pi / npix

predir = ut.rootpath + '/'
prefix_cat = predir + 'outputs/catalog_metacal_bin%d' % (o.bin_number)
prefix_map = predir + 'outputs/maps_metacal_bin%d_ns%d' % (o.bin_number,
                                                           o.nside)


ut.printflush("Computing means")
ut.printflush("0")
cat = fits.open(prefix_cat + '_zbin_mcal.fits')[1].data
ngal = len(cat)
e1_mean = np.mean(cat['e1'])
e2_mean = np.mean(cat['e2'])
e1_means = {}
e2_means = {}
for typ in ['1p', '1m', '2p', '2m']:
    ut.printflush(typ)
    f = fits.open(prefix_cat + '_zbin_mcal_' + typ + '.fits')
    e1_means[typ] = np.mean((f[1].data)['e1'])
    e2_means[typ] = np.mean((f[1].data)['e2'])
    f.close()

ut.printflush("Calibrating")
Rg = np.array([[np.mean(cat['R11']), np.mean(cat['R12'])],
               [np.mean(cat['R21']), np.mean(cat['R22'])]])
Rs = np.array([[(e1_means['1p']-e1_means['1m'])/0.02,
                (e1_means['2p']-e1_means['2m'])/0.02],
               [(e2_means['1p']-e2_means['1m'])/0.02,
                (e2_means['2p']-e2_means['2m'])/0.02]])
Rmat = Rg + Rs
one_plus_m = np.sum(np.diag(Rmat))*0.5
ut.printflush(Rg)
ut.printflush(Rs)
ut.printflush(Rmat)
ut.printflush(one_plus_m)
cat['e1'] = (cat['e1'] - e1_mean) / one_plus_m
cat['e2'] = (cat['e2'] - e2_mean) / one_plus_m

ut.printflush("Mapping")
ipix = hp.ang2pix(o.nside,
                  np.radians(90 - cat['dec']),
                  np.radians(cat['ra']))
map_w = np.bincount(ipix, minlength=npix).astype(float)
goodpix = map_w > 0
np.savez(prefix_map + "_goodpix.npz", pix=np.arange(npix)[goodpix].astype(int))
np.savez(prefix_map + "_w.npz", w=map_w[goodpix])


def get_ellip_maps(rot=False):
    if rot:
        phi = 2*np.pi*np.random.rand(ngal)
        c = np.cos(2*phi)
        s = np.sin(2*phi)
        e1 = c * cat['e1'] + s * cat['e2']
        e2 = -s * cat['e1'] + c * cat['e2']
    else:
        e1 = cat['e1']
        e2 = cat['e2']
    we1 = np.bincount(ipix, weights=e1, minlength=npix)[goodpix]
    we2 = np.bincount(ipix, weights=e2, minlength=npix)[goodpix]
    return we1, we2


map_we1, map_we2 = get_ellip_maps()
np.savez(prefix_map + "_we.npz", e1=map_we1, e2=map_we2)

map_wpsfe1 = np.bincount(ipix, weights=cat['psf_e1'], minlength=npix)[goodpix]
map_wpsfe2 = np.bincount(ipix, weights=cat['psf_e2'], minlength=npix)[goodpix]
np.savez(prefix_map + "_wpsfe.npz", e1=map_wpsfe1, e2=map_wpsfe2)

map_w2s2 = np.bincount(ipix, weights=0.5*(cat['e1']**2+cat['e2']**2),
                       minlength=npix)[goodpix]
np.savez(prefix_map + "_w2s2.npz", w2s2=map_w2s2)


ut.printflush("Rotations")
for i in range(o.nrot):
    ut.printflush("%d / %d " % (i, o.nrot))
    map_we1, map_we2 = get_ellip_maps(rot=True)
    np.savez(prefix_map + "_we_rot%d.npz" % i, e1=map_we1, e2=map_we2)

ut.printflush("N_ell quantities")
nl_bias = np.sum((cat['e1']**2+cat['e2']**2)*0.5) * pix_area / npix
nl_cov = nl_bias / np.mean(map_w**2)
np.savez(prefix_map + "_nells.npz", nl_bias=nl_bias, nl_cov=nl_cov)
