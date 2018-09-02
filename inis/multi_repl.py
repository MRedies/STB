import itertools
import numpy as np
from glob import glob 
prefix = ""#"/home/matthias/STB/inis/"

def strip_name(filename):
    filename = filename.split("/")[-1]
    filename = filename.split(".")[0]
    return filename

def replace_line(line, tags, values):
    for i in range(len(tags)):
        try:
            line = line.replace(tags[i], "{}".format(values[i]))
            # if(isinstance(values[i], np.int)):
            #     line = line.replace(tags[i], "{}".format(values[i]))
            # else:
            #     line = line.replace(tags[i], "{:08.4f}".format(values[i]))
        except:
            line = line.replace(tags[i], "{}".format(values))

    return line

    
tags = ["%filename%"]
# soc = np.arange(20,29, dtype=np.int)
# soc = np.concatenate( (soc, np.arange(0,20,4, dtype=np.int)) )
# print(soc,)

soc = glob("bobber_only*.cfg")
print(len(soc))
#for f in soc:
#    print(f)
#exit()

cnt = 0
itera = itertools.product(soc)
print("Number of files: {}".format(sum(1 for _ in itera)))
itera = itertools.product(soc)

for i in itera:
    with open(prefix + "jure.j", "rt") as fin:
        with open(prefix + "jrun_{}.j".format(i[0]), "wt") as fout: 
            for line in fin:
                fout.write(replace_line(line, tags, i))
    cnt +=1
