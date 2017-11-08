#!/usr/bin/env zsh
 
### Job name
#BSUB -J cluster_files/clx-test
   
### File / path where STDOUT & STDERR will be written
###    %J is the job ID, %I is the array ID
#BSUB -o cluster_files/clx-test.%J.%I
      
### Request the time you need for execution in minutes
### The format for the parameter is: [hour:]minute,
### that means for 80 minutes you could also use this: 1:20
#BSUB -W 5:00
          
### Request memory you need for your job per PROCESS in MB
#BSUB -M 126000 
               
### Hybrid Job with <N> MPI Processes in groups to <M> processes per node
# #BSUB -n <N>
# #BSUB -R "span[ptile=<M>]"
#BSUB -n 1
##BSUB -R "span[ptile=1]"
                     
### Request a certaion node type
#BSUB -m c24m128
#BSUB -P jara0062
                            
### Use nodes exclusive
#BSUB -x
                                    
### Each MPI process with T Threads
export OMP_NUM_THREADS=24
                                             
#BSUB -a intelmpi
module switch openmpi intelmpi
                                                        
### Change to the work directory
cd /home/pb321611/STB

#BSUB -u m.redies@fz-juelich.de
#BSUB -N
#BSUB -B

### Execute your application
echo "MPI FLAG"
echo $MPIEXEC $FLAGS_MPI_BATCH
$MPIEXEC $FLAGS_MPI_BATCH ./stb.x inis/example.cfg
