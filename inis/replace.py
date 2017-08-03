import numpy as np



atan = np.linspace(0.0,0.4, 10)
#atan = np.logspace(-2, 1.0,12)
#rad = np.arange(21,60,3)
prefix = "./inis/"

cnt = 0
for i in atan:
    with open(prefix + "example.cfg", "rt") as fin:
        with open(prefix + "rand_scan_{}.cfg".format(cnt), "wt") as fout:
            for line in fin:
                #fout.write(line.replace('%rad%', "{:04d}".format(i)))
                fout.write(line.replace('%rand%', "{:011.6f}".format(i)))
    cnt +=1
