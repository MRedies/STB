#!/usr/bin/python

import sys
import os
import shutil
import time

def get_folder(line):
    return line.split("=", 1)[1].strip().replace("'", "")

def mkdir(folder):
    if(not("%" in folder)):
        if(not(os.path.isdir(folder))):
            os.makedirs(folder)
            print "created " + folder
    else:
        print "didn't create " + folder + "because of %"

for filename in sys.argv[1:]:
    with open(filename) as f:
        lines = f.readlines()
        for l in lines:
            if "band_prefix" in l:
                mkdir(get_folder(l))
