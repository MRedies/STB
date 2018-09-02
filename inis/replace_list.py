import itertools
import numpy as np
from glob import glob 

def strip_name(filename):
    filename = filename.split("/")[-1]
    filename = filename.split(".")[0]
    return filename

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


def make_range_to_string(iterator):
    out_str = ""
    for i in iterator:
        out_str += "{} ,".format(i)  
    return out_str[:-1]


def all_but(exception):
    all = list(range(30))
    for e in exception:
        all.remove(e)
    return all
    
lower = list(range(21,28))
print(len(lower))

lower.extend(list(range(0,21, 5)))
print(len(lower))


for l in lower:
    u = l + 3
    exce = range(l,u)
    drop_layers = make_range_to_string(all_but(exce))

    with open("bobber_drop_stencil.cfg", "rt") as fin:
        with open("bobber_only_three_both_{:02d}_{:02d}.cfg".format(l,u), "wt") as fout:
            for line in fin:
                line = line.replace("%lower%", "{}".format(l))
                line = line.replace("%upper%", "{}".format(u-1))
                line = line.replace("%layerlist%", drop_layers)
                fout.write(line)