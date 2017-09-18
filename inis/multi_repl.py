import itertools
import numpy as np

prefix = "/home/matthias/STB/inis/"

def replace_line(line, tags, values):
    for i in range(len(tags)):
        #line = line.replace(tags[i], "{:08.4f}".format(values[i]))
        line = line.replace(tags[i], "{}".format(values[i]))
    return line

    
tags = ["%temp%", "%apd%"]

temp = np.linspace(0, 600.0, 15)
apd  = np.arange(5, 31,5)

cnt = 0
itera = itertools.product(temp, apd)
print "Number of files: {}".format(sum(1 for _ in itera))
itera = itertools.product(temp, apd)

for i in itera:
    val = np.array([i[0], i[1]])
    with open(prefix + "example.cfg", "rt") as fin:
        with open(prefix + "frank_{}.cfg".format(cnt), "wt") as fout:
            for line in fin:
                fout.write(replace_line(line, tags, val))
    cnt +=1
                
#fout.write(line.replace('%rad%', "{:04d}".format(i)))
#fout.write(line.replace('%rand%', "{:011.6f}".format(i)))


