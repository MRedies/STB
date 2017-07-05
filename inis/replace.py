import numpy as np
atan = np.linspace(14,17, 50)


prefix = "./inis/"

cnt = 0
for i in atan:
    with open(prefix+"example.cfg", "rt") as fin:
        with open(prefix + "chosen_param_atan_{}.cfg".format(cnt), "wt") as fout:
            for line in fin:
                #fout.write(line.replace('%apd%', "{:011.8f}".format(i)))
                fout.write(line.replace('%atan%', "{:011.6f}".format(i)))
    cnt +=1
