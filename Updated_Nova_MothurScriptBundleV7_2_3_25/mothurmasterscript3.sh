#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=24:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=8   # 8 processor core(s) per node
#SBATCH --mem=setme   # maximum memory per node
#SBATCH --mail-user=????@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --output="sb%j-mothurpt3"

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
## slurm break
module load mothur/1.48.0
mothur mscrip4.1.sh ##will run from cluster split to count groups, no rename files step, will also cover get.oturep,
## fin