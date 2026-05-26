Workflow for processing 16S rRNA gene amplicon sequencing data targeting the V4 region:
A majority of this pipeline was developed by Chiron J. Anderson, PhD, and Lucas R. Koester, PhD

Phase 1 - Mothur based pipeline for generating Operational Taxonomic Units (OTUs) from Illumina sequencing data
Phase 2 - Data analysis and visualization using Phyloseq in R
Phase 3 - Statistical analysis in SAS

#######################################################

Phase 1 - Mothur based pipeline for generating Operational Taxonomic Units (OTUs) from Illumina sequencing data

Note: This Mothur based pipeline is typically ran on a HPC cluster that uses Slurm to manage jobs on the computer nodes. Adjust accordingly if another job manager is used.
Additionally, ensure that Mothur is available on the HPC.(Check with the following command: module spider mothur) This pipeline has been confirmed to work with mothur/1.48.0. 

Step 1 - Initial quality check - Run: Step1_RunFastqc.sh
a. at a minimum check for:
	i. adapters (typically removed if sequencing performed at ISU DNA facility
	ii. # of raw reads per sample (for reporting)
	iii. sequence length
	iv. sequence quality

Step 2 - Preparing files to run Mothur pipeline:
a. Upload the following into the same directory as your raw .fastq files:
	i. vsearch-2.13.6-linux-x86_64
	ii. Contents of: Updated_Nova_MothurScriptBundleV7_2_3_25

b. Ensure that all files are unzipped: command: gunzip *.gz

The following pertains to files within Updated_Nova_MothurScriptBundleV7_2_3_25

c. open mscrip3.sh and update the path to the vsearch-2.13.6-linux-x86_64 folder in this command:
chimera.vsearch(vsearch=/work/PATH/vsearch-2.13.6-linux-x86_64/bin/vsearch, fasta=stability.trim.contigs.good.unique.good.filter.precluster.fasta, dereplicate=t, reference=silva.gold.fasta, processors=12)

d. open Slurmqueuer.sh and do the following:
	i. update email (occurs in two places in eh script)
	ii. adjust ntasks-per-node (#CPUs)
	iii. adjust mem = (RAM - occurs in two places in the script)
	
Step 3 - Run Mothur script (command: sbatch Slurmqueuer.sh)
a. the script will automatically submit multiple jobs and notify when each job is completed
b. it will also generate multiple log files summarizing information from each major script (mothurpt1, mothurpt2, mothurpt3)

step 4 - Confirming the Mothur run was successful:
a. Check initial and final summary.seq
	i. informs as to the length, quality, and number of sequences

b. Check the following in mothurpt1:
	i. length of filtered alignment: ###
		-acceptable values from 450-1400
		-if value is smaller that means too much was removed

c. Check the following in morthurpt2:
	i. check what % of your sequences were chimeras 
		- I've seen as high as 25%

Step 4 - Download and prepare outputs to use in R locally
a. download the file stability.trims.contigs.good.unique.good.filter.precluster.pick.opti_mcc0.01.cons.taxonomy 
	i. rename to: stability.cons.taxonomy
	ii. this file provides the SILVA taxonomy for the consensus sequence of each OTU

b. download the file stability.trims.contigs.good.unique.good.filter.precluster.pick.opti_mcc.shared 
	i. rename to: stability.opti_mcc.shared
	ii. this is your "count table" - identifies the number of reads from each OTU in each sample

c. download the file stability.trims.contigs.good.unique.good.filter.precluster.pick.opti_mcc0.01.rep.fasta
	i. this file provides the consensus sequence of each OTU
	ii. useful if you want to use other databases to assign taxonomy to a given OTU (e.g., NCBI)

#######################################################

Phase 2 - Data analysis and visualization using Phyloseq in R

a. Prepare local project directory and make three subdirectories named the following and move the following files into each directory:
	1. design
	- file: create metadata.csv file that contains the names of the samples/sequencing files and associated data in a long format
	2. shared
	-file: stability.opti_mcc.shared
	3. taxonomy
	-file: stability.cons.taxonomy
	
	- within the same project directory, download and place the following script: working_with_phyloseq_Decontam_4_12_23.R

b. Proceed with through the working_with_phyloseq_Decontam_4_12_23.R

#######################################################

Phase 3 - Statistical analysis in SAS

- Within the phyloseq pipeline, files will be generated that are intended to be uploaded and analyzed with SAS statistical software. 
- SAS is used to compare: 1. alpha diversity metrics, 2. changes in relative abundance on a genus and OTU level

a. When analyzing alpha diversity metrics, download and use the code: alpha_diversity_code_4-21-21.sas
	-Note that the figure associated with the alpha diversity metrics is generated from the primary phyloseq R script

b. When analyzing changes in relative abundance on an OTU level, download and use the code: OTU_comparison_5-3-22.sas
	-This SAS script generates a lot of outputs
	-Many of these outputs are used in the script responsible for the figures associated with this analysis: working_with_phyloseq_otu_comparison_7_19_21.R

#######################################################

