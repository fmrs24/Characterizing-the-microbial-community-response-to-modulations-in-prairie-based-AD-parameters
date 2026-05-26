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
#SBATCH --output="sb%j-mothurpt2"

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module use /opt/rit/spack-modules/lmod/linux-rhel7-x86_64/Core
module load mothur/1.43.0-py3-7cpdet4
mothur mscrip3.1.sh ##will take you through chimera.vsearch, stop just before cluster split
