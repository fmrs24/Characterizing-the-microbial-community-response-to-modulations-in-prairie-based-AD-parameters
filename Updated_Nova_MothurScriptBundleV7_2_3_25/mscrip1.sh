make.file(inputdir=./, type=fastq)
make.contigs(file=stability.files, processors=8)
summary.seqs(fasta=stability.trim.contigs.fasta, processors=8)