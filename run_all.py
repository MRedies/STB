from glob import glob
from subprocess import call

get_key = lambda str: int(str[5:-2])

files = sorted(glob("jrun_bobber_only*"))

if(len(files) == 12):
  for i,f in enumerate(files):
    print("{} run ".format(i) + f)
    call(["sbatch", f])