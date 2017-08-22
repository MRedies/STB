import numpy as np
from scipy.optimize import minimize
from subprocess import call
from termcolor import colored
import glob

prefix = "./inis/"
cnt = 0

def make_cfg(k):
    with open(prefix + "example.cfg", "rt") as fin:
        with open(prefix + "run.cfg", "wt") as fout:
            for line in fin:
                out_line = line.replace("%kx%", str(k[0]))
                out_line = out_line.replace("%ky%", str(k[1]))
                out_line = out_line.replace("%fphi%", str(k[2]))
                out_line = out_line.replace("%ftheta%", str(k[3]))
                fout.write(out_line)

def get_gap():
    E = np.load("output/opti/band_E.npy")
    # gap = np.max(E[1,:])
    # print colored('max val = {}'.format(gap), 'green') 
    # return -gap
    gap = np.min(E[2,:])
    print colored('min val = {}'.format(gap), 'green') 
    return gap

def f(x):
    make_cfg(x)
    call(["./stb.x", "inis/run.cfg"])
    print "x = {}".format(x)
    return get_gap()



bnd =((None, None), #kx
     (None, None), #ky
     (-np.pi, np.pi),
     (0.0, np.pi)) 
res = minimize(f, [0.1, 0.1, 0.3*np.pi, 0.2*np.pi], options={'disp': True}, bounds=bnd)

print "result: {}".format(res.x)

