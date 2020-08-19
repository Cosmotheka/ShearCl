import numpy as np
import healpy as hp
from argparse import ArgumentParser
import pymaster as nmt
import utils as ut
import os
import sys


parser = ArgumentParser()
parser.add_argument("--bin-number", default=0, type=int, help="Bin number")
parser.add_argument("--nside", default=4096, type=int, help="Nside")
parser.add_argument("--n-iter", default=0, type=int, help="n_iter")
parser.add_argument("--bin-number-2", default=-1, type=int, help="Bin number")
parser.add_argument("--is-psf-x", default=False, action='store_true',
                    help="Compute psf cross-correlation")
parser.add_argument("--is-psf-a", default=False, action='store_true',
                    help="Compute psf auto-correlation")
parser.add_argument("--irot-0", default=-1, type=int, help="Rotation number")
parser.add_argument("--irot-f", default=-1, type=int, help="Rotation number")
parser.add_argument("--recompute-mcm", default=False, action='store_true',
                    help="Recompute MCM even if it exists?")
o = parser.parse_args()

if o.bin_number_2 == -1:
    o.bin_number_2 = o.bin_number

is_auto = o.bin_number_2 == o.bin_number

npix = hp.nside2npix(o.nside)

predir = ut.rootpath + '/'


ut.printflush("Binning")
l_edges = np.array([0, 30, 60, 90, 120, 150, 180, 210, 240, 272, 309,
                    351, 398, 452, 513, 582, 661, 750, 852, 967, 1098,
                    1247, 1416, 1608, 1826, 2073, 2354, 2673, 3035,
                    3446, 3914, 4444, 5047, 5731, 6508, 7390, 8392,
                    9529, 10821, 12288])
l_edges = l_edges[l_edges <= 3*o.nside]
if 3*o.nside not in l_edges:
    l_edges = np.append(l_edges, 3*o.nside)
b = nmt.NmtBin.from_edges(l_edges[:-1], l_edges[1:])
l_eff = b.get_effective_ells()


ut.printflush("Weights")
prefix_map1 = predir + 'outputs/maps_metacal_bin%d_ns%d' % (o.bin_number, o.nside)
pix1 = np.load(prefix_map1 + '_goodpix.npz')['pix']
w1 = np.load(prefix_map1 + '_w.npz')['w']
w1_map = np.zeros(npix)
w1_map[pix1] = w1
if is_auto:
    prefix_map2 = prefix_map1
    pix2 = pix1
    w2 = w1
    w2_map = w1_map
else:
    prefix_map2 = predir + 'outputs/maps_metacal_bin%d_ns%d' % (o.bin_number_2, o.nside)
    pix2 = np.load(prefix_map2 + '_goodpix.npz')['pix']
    w2 = np.load(prefix_map2 + '_w.npz')['w']
    w2_map = np.zeros(npix)
    w2_map[pix2] = w2


ut.printflush("Fields")
def get_field(prefix_map, pix_arr, w_map, is_psf=False, i_rot=-1):
    if is_psf:
        fname_map = prefix_map + "_wpsfe"
    else:
        fname_map = prefix_map + "_we"
    if i_rot != -1:
        fname_map += "_rot%d" % i_rot
    fname_map += '.npz'
    d = np.load(fname_map)
    e1_map = np.zeros(npix)
    e1_map[pix_arr] = d['e1'] / w_map[pix_arr]
    e2_map = np.zeros(npix)
    e2_map[pix_arr] = d['e2'] / w_map[pix_arr]
    return nmt.NmtField(w_map, [-e1_map, e2_map], n_iter=o.n_iter)
f1 = get_field(prefix_map1, pix1, w1_map, is_psf=o.is_psf_a, i_rot=-1)
if is_auto and not o.is_psf_x:
    f2 = f1
else:
    f2 = get_field(prefix_map2, pix2, w2_map, is_psf=o.is_psf_x, i_rot=-1)


ut.printflush("MCM")
fname_mcm = predir + 'outputs/cls_metacal_mcm_'
fname_mcm += 'bins_%d%d_ns%d.fits' % (o.bin_number, o.bin_number_2, o.nside)
w = nmt.NmtWorkspace()
if os.path.isfile(fname_mcm) and not o.recompute_mcm:
    ut.printflush(" - Reading")
    w.read_from(fname_mcm)
else:
    ut.printflush(" - Computing")
    w.compute_coupling_matrix(f1, f2, b)
    w.write_to(fname_mcm)


def get_fname_cl(irot):
    fname_cls = predir + 'outputs/cls_metacal_cls_'
    fname_cls += 'bins_%d%d_ns%d' % (o.bin_number, o.bin_number_2, o.nside)
    if o.is_psf_x:
        fname_cls += '_xpsf'
    if o.is_psf_a:
        fname_cls += '_apsf'
    if irot != -1:
        fname_cls += '_rot%d' % irot
    fname_cls += '.npz'
    return fname_cls

if o.irot_0 == -1:
    ut.printflush("Cl")
    cls_coupled = nmt.compute_coupled_cell(f1, f2)
    cls = w.decouple_cell(cls_coupled)
    ut.printflush("Nl")
    if is_auto and not o.is_psf_x:
        s = np.load(prefix_map1 + '_nells.npz')
        nls_coupled = np.zeros([4, 3*o.nside])
        nls_coupled[0, 2:] = s['nl_bias']
        nls_coupled[3, 2:] = s['nl_bias']
        nls = w.decouple_cell(nls_coupled)
    else:
        nls = np.zeros_like(cls)

    ut.printflush("Writing")
    np.savez(get_fname_cl(-1), ls=l_eff, cls=cls, nls=nls, cls_coupled=cls_coupled)
else:
    ut.printflush("Rotations")
    for irot in range(o.irot_0, o.irot_f):
        ut.printflush("%d (%d %d)" % (irot, o.irot_0, o.irot_f))
        f = get_field(prefix_map1, pix1, w1_map, is_psf=False, i_rot=irot)
        cls = w.decouple_cell(nmt.compute_coupled_cell(f, f))
        nls = np.zeros_like(cls)
        np.savez(get_fname_cl(irot), ls=l_eff, cls=cls, nls=nls)
