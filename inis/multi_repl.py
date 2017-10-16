import itertools
import numpy as np

prefix = "/home/matthias/STB/inis/"

def replace_line(line, tags, values):
    for i in range(len(tags)):
        try:
            if(isinstance(values[i], np.int)):
                line = line.replace(tags[i], "{}".format(values[i]))
            else:
                line = line.replace(tags[i], "{:08.4f}".format(values[i]))
        except:
            line = line.replace(tags[i], "{}".format(values))

    return line

    
tags = ["%middle%"]

middle = np.linspace(0.5, 1.0, 10)
#atan = np.logspace(-1, 2, 10)
#middle = np.linspace(0.0, 1.0, 15)
#print apd
#temp = np.linspace(1,600, 5)

cnt = 0
itera = itertools.product(middle)
print "Number of files: {}".format(sum(1 for _ in itera))
itera = itertools.product(middle)

for i in itera:
    with open(prefix + "skyrm.cfg", "rt") as fin:
        with open(prefix + "middle_dos_{}.cfg".format(cnt), "wt") as fout:
            for line in fin:
                fout.write(replace_line(line, tags, i))
    cnt +=1
                
#fout.write(line.replace('%rad%', "{:04d}".format(i)))
#fout.write(line.replace('%rand%', "{:011.6f}".format(i)))


