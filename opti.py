import numpy as np
from scipy.optimize import minimize
from subprocess import call
from termcolor import colored

prefix = "./inis/"
cnt = 0

def make_cfg(k):
    global cnt
    with open(prefix + "example.cfg", "rt") as fin:
        with open(prefix + "run.cfg", "wt") as fout:
            for line in fin:
                out_line = line.replace("%kx%", str(k[0]))
                out_line = out_line.replace("%ky%", str(k[1]))
                out_line = out_line.replace("%atan%", str(abs(k[2])))
                out_line = out_line.replace("%d%", str(abs(k[3])))
                fout.write(out_line)

def get_gap():
    global cnt
    folder = "output/dbg/"
    E = np.load(folder + "band_E.npy")
    lower = E.shape[0]/2 -1
    upper = E.shape[0]/2 
    
    print "lower = {}; upper = {}".format(E[lower,0], E[upper,0])
    print colored("Gap = {}".format(E[upper,0] - E[lower,0]), "green")
    return np.abs(E[upper,0] - E[lower,0])

def f(x):
    global cnt
    cnt += 1
    make_cfg(x)
    call(["./stb.x", "inis/run.cfg"])
    print "k = {}".format(x)
    return get_gap()



bnd =((None, None), (None, None), #k
        (1.0, None), #a
        (0.05, 0.95)) #t_so
res = minimize(f, [0.0, 0.0, 13.2898, 0.8], options={'disp': True}, bounds=bnd)

print "result: {}".format(res.x)

