import wget
import os
import utils as ut


def dwl_file(fname, url, call=None,
             verbose=True):
    if not os.path.isfile(fname):
        if verbose:
            print(fname)
        wget.download(url)
        print("\n")
        if call is not None:
            call()


def mkdir(dr):
    if not os.path.isdir(dr):
        os.makedirs(dr)


def unzip(fname):
    os.system('unzip '+fname)
    os.remove(fname)


predir = os.path.abspath(ut.rootpath)

# Create data directory
mkdir(predir+"/DES_data")
os.chdir(predir+"/DES_data")

# Shear sample
pre_url = "http://desdr-server.ncsa.illinois.edu/despublic/y1a1_files/"
mkdir("shear_catalog")
os.chdir("shear_catalog")
dwl_file("mcal-y1a1-combined-riz-unblind-v4-matched.fits",
         pre_url + "shear_catalogs/" +
         "mcal-y1a1-combined-riz-unblind-v4-matched.fits")
dwl_file("y1a1-im3shape_v5_unblind_v2_matched_v4.fits",
         pre_url + "shear_catalogs/" +
         "y1a1-im3shape_v5_unblind_v2_matched_v4.fits")
dwl_file("y1_source_redshift_binning_v1.fits",
         pre_url + "redshift_bins/y1_source_redshift_binning_v1.fits")
dwl_file("y1_redshift_distributions_v1.fits",
         pre_url + "redshift_bins/y1_redshift_distributions_v1.fits")
os.chdir("../")

# Go back to root
os.chdir(predir)
