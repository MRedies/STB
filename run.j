#!/usr/bin/env zsh
 
### Job name
#BSUB -J cluster_files/mpi-s-test
   
### File / path where STDOUT & STDERR will be written
###    %J is the job ID, %I is the array ID
#BSUB -o cluster_files/mpi-s-test.%J.%I
      
### Request the time you need for execution in minutes
### The format for the parameter is: [hour:]minute,
### that means for 80 minutes you could also use this: 1:20
#BSUB -W 5:00
          
### Request memory you need for your job per PROCESS in MB
#BSUB -M 21000 
               
### Hybrid Job with <N> MPI Processes in groups to <M> processes per node
# #BSUB -n <N>
# #BSUB -R "span[ptile=<M>]"
#BSUB -n 5
##BSUB -R "span[ptile=1]"
                     
### Request a certaion node type
#BSUB -m mpi-s
                            
### Use nodes exclusive
#BSUB -x
                                    
### Each MPI process with T Threads
                                    export OMP_NUM_THREADS=12
                                             
#BSUB -a intelmpi
                                                        
### Change to the work directory
                                                        cd /home/mr071525/STB

#BSUB -u m.redies@fz-juelich.de
#BSUB -N
#BSUB -B

### Execute your application
$MPIEXEC $FLAGS_MPI_BATCH /home/mr071525/STB/stb.x /home/mr071525/STB/inis/example.cfg

