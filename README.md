When cloning the project you need to use
```bash
git clone --recursive
```
so that the submodules are cloned aswell.

# Install

## Install fortran_stdlib

See here for install guide: https://github.com/fortran-lang/stdlib .
Remember to select the correct compiler via the ```FC``` variable and to set ```-DCMAKE_INSTALL_PREFIX``` to a path where you have necessary permissions. ```./cmake/Modules/stdlib/``` is a good choice.

## Install third-party

You need to install the package in ./thirdparty/triangulation by calling ```cmake``` in that directory.

## Install STB code

Export the path where you installed the ```stdlib``` package, for example ```export CMAKE_PREFIX_PATH=$PROJECT/kipp1/STB/cmake/Modules/stdlib:$CMAKE_PREFIX_PATH```.
Then run cmake in a build folder to avoid clutter, ```cmake ../ -DPLATFORM=GNU -DDEBUG=off```, selecting the correct compiler for the platform via the ```-DPLATFORM``` flag and choosing the compiler optimization level via the ```-DDEBUG``` flag (see ```./CMakeLists.txt```).
Then run ```make```.

# Example batch file

```
#!/bin/bash -x
#SBATCH--account=jiff40
#SBATCH--nodes=32
#SBATCH--ntasks-per-node=128
#SBATCH--cpus-per-task=1
#SBATCH--output=mpi-out.%j
#SBATCH--error=mpi-err.%j
#SBATCH--time=24:00:00
#SBATCH--partition=dc-cpu
#SBATCH--mail-user=j.kipp@fz-juelich.de
#SBATCH--mail-type=FAIL
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
module use $OTHERSTAGES
module load Stages/2022
ml GCC CMake Python ParaStationMPI imkl

module list
#ldd /p/home/jusers/kipp1/jureca/test/STB/stb.x

srun /p/project/cjiff40/kipp1/STB/build/stb.x /p/project/cjiff40/kipp1/ML_Samples/STB_Samples/random_transport/random.cfg
```