import numpy as np
kpd = np.arange(200,1200,200)
print(kpd)

prefix = "/home/matthias/STB/inis/"

for i in kpd:
    with open(prefix+"example.cfg", "rt") as fin:
        with open(prefix + "k_scan={}.cfg".format(i), "wt") as fout:
            for line in fin:
                fout.write(line.replace('%kpd%', "{}".format(i)))
