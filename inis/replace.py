import numpy as np



#atan = np.linspace(0.01,10, 200)
#atan = np.linspace(0.01,25,150)
kpd = np.arange(2,100)
prefix = "./inis/"

cnt = 0
for i in kpd:
    with open(prefix+"example.cfg", "rt") as fin:
        with open(prefix + "kpd_{}.cfg".format(cnt), "wt") as fout:
            for line in fin:
                fout.write(line.replace('%kpd%', "{:04d}".format(i)))
                #fout.write(line.replace('%kpd%', "{:011.6f}".format(i)))
    cnt +=1
