import numpy as np
#Ef = range(10,401,10)
Ef = np.linspace(-40.0, 40.0, 50)

prefix = "./inis/"

cnt = 0
for i in Ef:
    with open(prefix+"example.cfg", "rt") as fin:
        with open(prefix + "J2000_scan{}.cfg".format(cnt), "wt") as fout:
            for line in fin:
                fout.write(line.replace('%Ef%', "{0:010.4f}".format(i)))
    cnt +=1
