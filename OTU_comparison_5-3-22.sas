**SAS Pipeline Code for analyzing count data (i.e.Microbiome)using a negative binomal distrubution, Edited by Amy Petry**	
	Steps needed prior to using this code:
		1. You need the following files in your folder pipeline to compelte this code 
			A. Shared CSV with the counts of the top OTUS you want to analysis. Use the R Code Pipeline to obtain this 
			B. Taxonomic information from the R code pipeline for the top # of OTUS you are intrested in. 
			C. Metadata file with the count data included and log transformed. All count data must be log transformed**;
**Run the following a Step at a time completing all of the edit sections;

**Step 1: Setting up your sas log to stretch out Pvalues to 18 decimal places, EDITS: none;

proc template; 			
define column common.pvalue;
notes "default p-value column";
just = r;
format = pvalue20.18;
ods listing;
ods graphics off;
end;
run;
** The below code turns off the template step above, this is optinal, but in case you need to change the settings, Just need to remove the '*' in front of the lines to do so; 
*proc template;    
*delete common.pvalue;
*run;

**Step 1: Importing the 3 files you need to run the analysis. Complete the edit steps before running in each section;

FILENAME REFFILE 'C:\Users\PATH\shared_sas.csv';** edit this pathway to your specfic pathway for the sharedfile;
**Inputing the refile from above;
PROC IMPORT DATAFILE=REFFILE 
    DBMS=CSV
    OUT=Work.Shared; **Outputing this info into your working directory; 
    GETNAMES=YES;
RUN;

FILENAME REFFILE 'C:\Users\PATH\taxonomy_100_sas.csv'; ** edit this pathway to your specfic pathway for the taxonomy file;

PROC IMPORT DATAFILE=REFFILE
    DBMS=CSV
    OUT=Work.Taxonomy; **Outputing this info into your working directory; 
    GETNAMES=YES;
RUN;

*before reading in metadata file, make sure that the count data was been log transformed (ln function in excel). Otherwise, this code will NOT run!*;

FILENAME REFFILE 'C:\Users\PATH\metadata_sas.csv'; ** edit this pathway to your specfic pathway for the metadatafile;

PROC IMPORT DATAFILE=REFFILE
    DBMS=CSV
    OUT=Work.Metadata;**Outputing this info into your working directory; 
    GETNAMES=YES;
RUN;
**Step 2:Checking the inputs abouve by printing the data into an output**;
Proc Print data=Shared; Run; 
Proc Print data=Taxonomy; Run;
Proc Print data=Metadata; Run;Quit;
**Step 3: condesing the files into one working directory**;
data work.Data;
    set work.Shared; *Create a shared file for each location with ONLY the top 200 OTUs. Otherwise, it will not read in OTU counts.*;
    array x(i) Otu000001-Otu000100;  *EDIT:depending on the number of OTUs that were present, the number of zeros in these numbers change. Make sure they match the column headings for shared file*;
do over x;
        y = x;
        output;
    end;
    keep Group i y;
    rename i=OTU;
 run;
**Step 4 Creating a table to analyze by OTU and including the metadata analysis;
PROC SQL;
Create table work.DataRun As
	Select 
	B.SampleID, B.INPUT, B.INPUT, B.INPUT, B.INPUT, B.Count, A.OTU, A.y
From work.Data as A Inner join work.Metadata as B
On A.Group=B.SampleID;
quit;
**Step 5 analysing the data using a negative binomal, you need to edit the class statement to the appropriate model for your study**;
Proc Sort data=Datarun; by OTU; Run; 
title "INPUT_OTU_comparison";
proc glimmix data=Datarun initglm maxopt=100 pconv=1e7; 
where OTU in (1:100); 
class INPUT INPUT INPUT; *include all variables to be included in the model;
model y = INPUT INPUT INPUT*INPUT / dist=negbinomial ddfm=KR offset=Count alpha=0.05; *include all the variables you would like to include, including interactions *;
*ABOVE: you can change K(ENWARD)R(OGERS) to SATTERTHWAITE or K(ENWARD)R(OGERS) depending on your experimental set-up;
random _residual_ / subject=INPUT type=CS; *for the repeated variable (ex. sheep_ID). You can select different types for type and check the BCC value for each. Lowest is best;
*random  / subject=INPUT type=CS; *this is for including simple random varible effects (ie not repeated);
output out=work.pred pred=p resid=r;
lsmeans INPUT INPUT INPUT*INPUT/ pdiff=all; *include the same variables in the model line here;
by OTU;
nloptions maxiter=100; 
ods listing exclude diffs lsmeans;
ods output tests3=work.pvalues; 
ods output FitStatistics=work.fitstats;
ods output ParameterEstimates=work.parms;
ods output LSMeans=work.LSM; *use the estimates in this file to get reportable values. In excel, use exp function and then multiply by 100.*;
ods output diffs=work.pdfiff; *Use this to differentiate superscripts*;
run; 
**Step 6 obtaining the CI to calcualte the SEM**;
proc univariate data = work.LSM cibasic;
var Estimate;
by OTU;
ods output BasicIntervals = work.CI;  *this is the file that will contain confidence intervals.*; 
run;
**Step 7 Exporting work files into one excel sheet**; **EDIT: the outfile statement to the location you want the file and the name**; 
Proc Export data=work.pvalues
	Outfile="C:\Users\PATH\pvalues.csv"
	DBMS=CSV;
	putnames=yes;
	run;
	Proc Export data=work.fitstats
	Outfile="C:\Users\PATH\fitstats.csv"
	DBMS=CSV;
	putnames=yes;
	run;
	Proc Export data=work.parms
	Outfile="C:\Users\PATH\parms.csv"
	DBMS=CSV;
	putnames=yes;
	run;
	Proc Export data=work.LSM
	Outfile="C:\Users\PATH\LSM.csv"
	DBMS=CSV;
	putnames=yes;
	run;
	Proc Export data=work.pdfiff
	Outfile="C:\Users\PATH\pdiff.csv"
	DBMS=CSV;
	putnames=yes;
	run;
	Proc Export data=work.CI
	Outfile="C:\Users\PATH\CI.csv"
	DBMS=CSV;
	putnames=yes;
	run;
* Use exp function and multiply by 100 in excel. Then convert to standard error = (upper limit-lower limit)/3.92*;
**Step 8: Obtaininig Q values using false discovery rates**;
	Proc export data=work.pvalues
	Outfile="C:\Users\PATH\pv.csv"
	DBMS=CSV;
	putnames=yes;
	run;
data qvalues;
	infile "C:\Users\PATH\pv.csv" dlm = "," firstobs = 2;
	input OTU Effect$ NumDF DenDF FValue RAW_P;  run;
	Proc print data=qvalues; run; quit;
	proc multtest inpvalues=qvalues hom hoc fdr;
	ods output pValues=QV; 
	run; 
	Proc Export data=work.QV
	Outfile="C:\Users\PATH\qvalues.csv"
	DBMS=CSV;
	putnames=yes;
	run; quit; 
