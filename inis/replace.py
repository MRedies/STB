import numpy as np
I = np.linspace(0,1,5)

prefix = "/home/matthias/MasterCode/STB/inis/"
for i in I:
    with open(prefix+"example.cfg", "rt") as fin:
        with open(prefix + "stoner={:010.8f}.cfg".format(i), "wt") as fout:
            for line in fin:
                fout.write(line.replace('%stoner%', "{:010.8f}".format(i)))
