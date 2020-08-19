import healpy as hp
import numpy as np


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


def get_weighted_sums(iterator, name_weight, names_field,
                      masks=None):
    if masks is None:
        nmaps = 1
        masks = [[['all']]]
    else:
        nmaps = len(masks)

    counts = np.zeros(nmaps)
    if name_weight is not None:
        weights = np.zeros(nmaps)
    else:
        weights = None
    if names_field is not None:
        nfields = len(names_field)
        fields = np.zeros([nmaps, nfields])

    for i_d, d in enumerate(iterator):
        for im, mm in enumerate(masks):
            mask = np.ones(len(d['ra']), dtype=bool)
            print(len(mask))
            for m in mm:
                if m[0] == 'tag':
                    mask = mask & (d[m[1]] == m[2])
                elif m[0] == 'range':
                    mask = mask & (d[m[1]] < m[3]) & (d[m[1]] >= m[2])

            counts[im] += np.sum(mask)

            if name_weight is not None:
                w = d[name_weight][mask]
                weights[im] += np.sum(w)
                if names_field is not None:
                    for i_f, n_f in enumerate(names_field):
                        f = d[n_f][mask]
                        fields[im, i_f] += np.sum(w * f)
            else:
                if names_field is not None:
                    for i_f, n_f in enumerate(names_field):
                        f = d[n_f][mask]
                        fields[im, i_f] += np.sum(f)

    if nmaps == 1:
        counts = counts[0]
        if weights is not None:
            weights = weights[0]
        if fields is not None:
            fields = fields[0]
    return counts, weights, fields


def get_weighted_maps(iterator, nside, name_ra, name_dec,
                      name_weight=None, names_field=None,
                      masks=None):

    npix = hp.nside2npix(nside)

    if masks is None:
        nmaps = 1
        masks = [[['all']]]
    else:
        nmaps = len(masks)

    map_counts = np.zeros([nmaps, npix])
    if name_weight is not None:
        map_weights = np.zeros([nmaps, npix])
    else:
        map_weights = None
    if names_field is not None:
        nfields = len(names_field)
        map_field = np.zeros([nmaps, nfields, npix])
    else:
        map_field = None

    for i_d,d in enumerate(iterator):
        ipix = hp.ang2pix(nside,
                          np.radians(90 - d[name_dec]),
                          np.radians(d[name_ra]))
        print(len(ipix))
        for im, mm in enumerate(masks):
            mask = np.ones(len(ipix), dtype=bool)
            for m in mm:
                if m[0] == 'tag':
                    mask = mask & (d[m[1]] == m[2])
                elif m[0] == 'range':
                    mask = mask & (d[m[1]] < m[3]) & (d[m[1]] >= m[2])

            ip = ipix[mask]
            map_counts[im, :] += np.bincount(ip,
                                             minlength=npix)
            if name_weight is not None:
                w = d[name_weight][mask]
                map_weights[im, :] += np.bincount(ip,
                                                  minlength=npix,
                                                  weights=w)
                if names_field is not None:
                    for i_f, n_f in enumerate(names_field):
                        f = d[n_f][mask]
                        map_field[im, i_f, :] += np.bincount(ip,
                                                             minlength=npix,
                                                             weights=w * f)
            else:
                if names_field is not None:
                    for i_f, n_f in enumerate(names_field):
                        f = d[n_f][mask]
                        map_field[im, i_f, :] += np.bincount(ip,
                                                             minlength=npix,
                                                             weights=f)

    if nmaps == 1:
        map_counts = map_counts.flatten()
        if map_weights is not None:
            map_weights = map_weights.flatten()
        if map_field is not None:
            map_field = map_field.reshape([nfields, npix])

    return map_counts, map_weights, map_field


def rotate_alm_g_c(alm_in, c2g=False):
    if c2g:
        coord=['C','G']
    else:
        coord=['G','C']

    r=hp.Rotator(coord=coord)
    return r.rotate_alm(alm_in)


def rotate_map_g_c(map_in, c2g=False):
    ns = hp.npix2nside(len(map_in))
    alm_in = hp.map2alm(map_in)
    alm_out = rotate_alm_g_c(alm_in, c2g=False)
    return hp.alm2map(alm_out, ns, verbose=False)
