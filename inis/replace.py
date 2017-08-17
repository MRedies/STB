import numpy as np



#atan = np.linspace(1.95, 2.05, 10)
#atan = np.arange(10,360,5)
#atan = np.logspace(-2, 1.0,12)
#rad = np.arange(21,60,3)

# theta = np.load("m_theta.npy")
# phi   = np.load("m_phi.npy")

r = np.linspace(0.1, 0.9, 7)


prefix = "/home/matthias/STB/inis/"

cnt = 0
for mid in r:
    with open(prefix + "example.cfg", "rt") as fin:
        with open(prefix + "shift_middle_{}.cfg".format(cnt), "wt") as fout:
            for line in fin:
                #line = line.replace('%m_phi%',   "{:011.6f}".format(phi[i]))
                #line = line.replace('%m_theta%', "{:011.6f}".format(theta[i]))
                line = line.replace('%middle%', "{:011.6f}".format(mid))
                fout.write(line)
    cnt +=1
