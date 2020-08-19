import numpy as np
from argparse import ArgumentParser
import pymaster as nmt
import utils as ut


parser = ArgumentParser()
parser.add_argument("--bin-number", default=0, type=int, help="Bin number")
parser.add_argument("--nside", default=4096, type=int, help="Nside")
parser.add_argument("--bin-number-2", default=-1, type=int, help="Bin number")
o = parser.parse_args()

if o.bin_number_2 == -1:
    o.bin_number_2 = o.bin_number

predir = ut.rootpath + '/'

ut.printflush("MCM")
fname_mcm = predir + 'outputs/cls_metacal_mcm_'
fname_mcm += 'bins_%d%d_ns%d.fits' % (o.bin_number, o.bin_number_2, o.nside)
w = nmt.NmtWorkspace()
w.read_from(fname_mcm)

ut.printflush("Win")
win = w.get_bandpower_windows()

ut.printflush("Writing")
fname_win = predir + 'outputs/cls_metacal_win_'
fname_win += 'bins_%d%d_ns%d.npz' % (o.bin_number, o.bin_number_2, o.nside)
np.savez(fname_win, win=win)
