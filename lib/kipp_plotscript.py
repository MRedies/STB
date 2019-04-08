#plots data from output file
import post_proc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib import cm,colors
from matplotlib.backends.backend_pdf import PdfPages
from kipp_plotbib import bandplotter
prefix = "/Users/kipp/STB/output/"
save_prefix = "/Users/kipp/STB/weyl_figures/Weyl_figures/"
filename = "path_rel_G-K-Kprime_anticol_anticol"

bandplotter(fname = filename,files = 10,bands = [0,1,2,3],combine = False,hall = True)
bandplotter(fname = filename,files = 10,bands = [0,1,2,3],combine = True,hall = True)
