import numpy as np
#Ef = np.linspace(0.0, 0.5*np.pi, 10)

apd = np.arange(1,26, 5)

prefix = "./inis/"

cnt = 0
for i in apd:
    with open(prefix+"example.cfg", "rt") as fin:
        with open(prefix + "skyrm_size{}.cfg".format(cnt), "wt") as fout:
            for line in fin:
                #fout.write(line.replace('%apd%', "{:011.8f}".format(i)))
                fout.write(line.replace('%apd%', "{}".format(i)))
    cnt +=1
