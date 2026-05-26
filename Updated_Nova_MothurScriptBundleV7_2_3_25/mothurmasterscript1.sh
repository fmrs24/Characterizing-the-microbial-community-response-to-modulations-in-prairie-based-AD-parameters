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
#SBATCH --output="sb%j-mothurpt1"

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load mothur/1.48.0
gzip *.logfile ### gets any existing mothur logs out of the way for greps later
mothur mscrip1.1.sh ## will go up to 1st summary step and exit
MinLength=$(grep -n -i '2.5%' *.logfile | awk '{print $3}' ) ## Defines the MinLength variable based on the length of the 2.5% row from the last mothur log file (the third column)
MaxLength=$(grep -n -i '97.5' *.logfile | awk '{print $3}' ) ## Defines the MaxLength variable based on the length of the 97.5% row from the last mothur log file. for miseq data which is highly regular, this should be appropriate
## ^ and will help remove sequences which did not overlap properly during the make.contigs command
## now lets write the next script using the min length variable we just created. This script will finish with the alignment step.
echo $"screen.seqs(fasta=stability.trim.contigs.fasta, count=stability.contigs.count_table, maxambig=0, minlength=$MinLength, maxlength=$MaxLength, maxhomop=8, processors=8)
summary.seqs(fasta=stability.trim.contigs.good.fasta, processors=8)
count.groups(count=stability.contigs.good.count_table)
unique.seqs(fasta=stability.trim.contigs.good.fasta, count=stability.contigs.good.count_table)
summary.seqs(fasta=stability.trim.contigs.good.unique.fasta, count=stability.trim.contigs.good.count_table)
align.seqs(fasta=stability.trim.contigs.good.unique.fasta, reference=silva.nr_v138.align, processors=8)" > mscrip1.5.sh
##now lets run that script
mothur mscrip1.5.sh ## goes up to 3rd summary, break to gzip logs out of the way
gzip *.logfile ## zip zip
mothur mscrip2.1.sh ## goes to 4th summary which we need to use to set start and end points for screening the aligned files.
##lets grab those start and end points from the log file
Start=$(grep -n -i '97.5%' *.logfile | awk '{print $2}') ##pulls start from 97.5% - this is preferable due to way start is defined in mothur
End=$(grep -n -i '2.5%' *.logfile | awk '{print $3}') ##pulls end from 2.5% - this is preferable due to way start is defined in mothur
## now lets write the last set of commands before the chimera check
echo "screen.seqs(fasta=stability.trim.contigs.good.unique.align, count=stability.trim.contigs.good.count_table, start=$Start, end=$End, processors=8)
summary.seqs(fasta=stability.trim.contigs.good.unique.good.align, count=stability.trim.contigs.good.good.count_table, processors=8)
count.groups(count=stability.trim.contigs.good.good.count_table)
filter.seqs(fasta=stability.trim.contigs.good.unique.good.align, vertical=T, trump=., processors=8)
pre.cluster(fasta=stability.trim.contigs.good.unique.good.filter.fasta, count=stability.trim.contigs.good.good.count_table, diffs=2, processors=8)
summary.seqs(fasta=stability.trim.contigs.good.unique.good.filter.precluster.fasta, count=stability.trim.contigs.good.unique.good.filter.precluster.count_table, processors=8)
count.groups(count=stability.trim.contigs.good.unique.good.filter.precluster.count_table)" > mscrip2.5.sh
mothur mscrip2.5.sh ##goes up to chimera check