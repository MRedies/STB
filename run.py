#!/home/matthias/anaconda3/bin/python

import sys
from subprocess import call

call(["ipython", "inis/replace.py"])

for i in range(1,len(sys.argv)):
    print("Run: " + sys.argv[i])
    arg = str(sys.argv[i])
    command = "./stb.x " + arg

    call(command.split(), shell=False)

