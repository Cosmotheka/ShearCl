import numpy as np
from astropy.io import fits
from argparse import ArgumentParser
import utils as ut


parser = ArgumentParser()
parser.add_argument("--bin-number", default=0, type=int, help="Bin number")
o = parser.parse_args()


nrows_per_chunk = 1000000
predir = ut.rootpath + '/'
fname_mcal = predir + 'DES_data/shear_catalog/mcal-y1a1-combined-riz-unblind-v4-matched.fits'
fname_bins = predir + 'DES_data/shear_catalog/y1_source_redshift_binning_v1.fits'
prefix_out = predir + 'outputs/catalog_metacal_bin%d' % (o.bin_number)


def get_fits_iterator(fname, colnames, hdu=1, nrows_per_chunk=None):
    import fitsio

    # Open file and select HDU
    fitf = fitsio.FITS(fname, mode='r')
    tab = fitf[hdu]

    # Count rows and number of chunks
    nrows = tab.get_nrows()
    # If None, then just one chunk
    if nrows_per_chunk is None:
        nrows_per_chunk = nrows
    # Add one chunk if there are leftovers
    nchunks = nrows // nrows_per_chunk
    if nrows_per_chunk * nchunks < nrows:
        nchunks += 1

    for i in range(nchunks):
        start = i * nrows_per_chunk
        end = min((i + 1) * nrows_per_chunk, nrows)
        data = tab.read_columns(colnames,
                                rows=range(start, end))
        yield data


def get_iterator_metacal():
    itr_cat = get_fits_iterator(fname_mcal,
                                ['coadd_objects_id', 'e1', 'e2', 'psf_e1', 'psf_e2',
                                 'ra', 'dec', 'flags_select', 'R11', 'R22', 'R12', 'R21'],
                                nrows_per_chunk=nrows_per_chunk)
    itr_bin = get_fits_iterator(fname_bins,
                                ['coadd_objects_id', 'zbin_mcal',
                                 'zbin_mcal_1p', 'zbin_mcal_1m',
                                 'zbin_mcal_2p', 'zbin_mcal_2m'],
                                nrows_per_chunk=nrows_per_chunk)
    for m, b in zip(itr_cat, itr_bin):
        dc = {'coadd_objects_id': m['coadd_objects_id'],
              'ra': m['ra'],
              'dec': m['dec'],
              'e1': m['e1'],
              'e2': m['e2'],
              'psf_e1': m['psf_e1'],
              'psf_e2': m['psf_e2'],
              'flags_select': m['flags_select'],
              'R11': m['R11'],
              'R12': m['R12'],
              'R21': m['R21'],
              'R22': m['R22'],
              'zbin_mcal': b['zbin_mcal'],
              'zbin_mcal_1p': b['zbin_mcal_1p'],
              'zbin_mcal_1m': b['zbin_mcal_1m'],
              'zbin_mcal_2p': b['zbin_mcal_2p'],
              'zbin_mcal_2m': b['zbin_mcal_2m']}
        yield dc


def create_subsamples(iterator, fnames, cols_save=None, masks=None):
    if masks is None:
        nsamples = 1
        masks = [[['all']]]
    else:
        nsamples = len(masks)

    arrs = [{} for i in range(nsamples)]
    for i_d, d in enumerate(iterator):
        print(i_d)
        for im, mm in enumerate(masks):
            mask = np.ones(len(d['ra']), dtype=bool)
            for m in mm:
                if m[0] == 'tag':
                    mask = mask & (d[m[1]] == m[2])
                elif m[0] == 'range':
                    mask = mask & (d[m[1]] < m[3]) & (d[m[1]] >= m[2])
            if i_d == 0:
                for k in d.keys():
                    arrs[im][k] = d[k][mask]
            else:
                for k in d.keys():
                    arrs[im][k] = np.concatenate((arrs[im][k], d[k][mask]))

    for im, mm in enumerate(masks):
        if cols_save is None:
            cols_save = arrs[im].keys()
        fmt_dir = {'int64': 'K',
                   'float64': 'D'}
        cols = fits.ColDefs([fits.Column(name=k,
                                         format=fmt_dir[str(arrs[im][k].dtype)],
                                         array=arrs[im][k])
                             for k in cols_save])
        hdu = fits.BinTableHDU.from_columns(cols)
        hdu.writeto(fnames[im])

binflags = ['zbin_mcal', 'zbin_mcal_1p',
            'zbin_mcal_1m', 'zbin_mcal_2p',
            'zbin_mcal_2m']
masks = [[['tag', bf, o.bin_number],
          ['tag', 'flags_select', 0],
          ['range', 'dec', -90., -35.]]
         for bf in binflags]
fnames = [prefix_out + '_' + bf + '.fits' for bf in binflags]
cols_save = ['coadd_objects_id', 'ra', 'dec', 'e1', 'e2',
             'psf_e1', 'psf_e2', 'R11', 'R12', 'R21', 'R22']

it = get_iterator_metacal()
create_subsamples(it, fnames, cols_save=cols_save, masks=masks)
