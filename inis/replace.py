prefix = "/home/matthias/MasterCode/STB/inis/"
for i in range(1,10):
    with open(prefix+"example.cfg", "rt") as fin:
        with open(prefix+"square=%d.cfg"%(i), "wt") as fout:
            for line in fin:
                fout.write(line.replace('%n%', str(i)))
