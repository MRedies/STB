for i in range(1,10):
    with open("example.cfg", "rt") as fin:
        with open("path=%d.cfg"%(i), "wt") as fout:
            for line in fin:
                fout.write(line.replace('%hs%', str(i)))
