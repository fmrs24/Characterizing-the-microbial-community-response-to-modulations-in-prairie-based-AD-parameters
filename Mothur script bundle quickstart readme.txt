Mothur script bundle quickstart readme
Steps:
1) Copy the relevant script bundle somewhere locally on your computer (so you don't mess with the public copy)
2) Unzip the folder
3) Upload the full contents of the folder to the directory with your reads
4) Use gunzip *.gz on prontodtn to unzip the reference databases
5) Open the Slurmqueuer.sh file on the server (with winSCP or nano etc)
6) Change the email address to your own on lines 11 (in the slurm block) and 16 (variables).
7) Optionally, change any of the other parameters on lines 17-21. These set the ram, processors, and max wall times of the 3 segments of the mothur script.
	For most projects the defaults should be fine as is, but for particularly large datasets you may find you need to increase the memory
8) Save your changes, then start your job with sbatch Slurmqueuer.sh (on pronto, not the dtn)
9) You will get updates as the different parts of the job finish (if you set the email properly), but you can monitor progress with
	squeue | grep "<your netID>"
10) Don't forget to check the logs and make sure everything ran properly before moving on to your analysis!