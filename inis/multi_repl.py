import itertools
import numpy as np

prefix = "./inis/"

def replace_line(line, tags, values):
    for i in range(len(tags)):
        line = line.replace(tags[i], "{:08.4f}".format(values[i]))
    return line

    
tags = ["%ferro_theta%", "%ferro_phi%"]

fer_th = np.linspace(0.0, 0.5*np.pi, 50)
fer_phi = np.linspace(-np.pi, np.pi, 50)

cnt = 0
itera = itertools.product(fer_th, fer_phi)
print "Number of files: {}".format(sum(1 for _ in itera))
itera = itertools.product(fer_th, fer_phi)

for i in itera:
    val = np.array([i[0], i[1]])
    with open(prefix + "example.cfg", "rt") as fin:
        with open(prefix + "pt1_scan_{}.cfg".format(cnt), "wt") as fout:
            for line in fin:
                fout.write(replace_line(line, tags, val))
    cnt +=1
                
#fout.write(line.replace('%rad%', "{:04d}".format(i)))
#fout.write(line.replace('%rand%', "{:011.6f}".format(i)))


