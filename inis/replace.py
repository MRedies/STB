#import numpy as np
num_k = range(10,401,10)
print(num_k)

prefix = "./inis/"

for i in num_k:
    with open(prefix+"example.cfg", "rt") as fin:
        with open(prefix + "k_scan{}.cfg".format(i), "wt") as fout:
            for line in fin:
                fout.write(line.replace('%num_k%', "{}".format(i)))
