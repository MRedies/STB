#!/home/matthias/anaconda2/bin/ipython

import sys
from subprocess import call

call(["ipython", "inis/replace.py"])

for i in range(1,len(sys.argv)):
    print("Run: " + sys.argv[i])
    arg = str(sys.argv[i])
    command = "mpirun -np 2 ./stb.x " + arg

    call(command.split(), shell=False)

