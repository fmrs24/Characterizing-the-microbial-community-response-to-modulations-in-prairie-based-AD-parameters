* Alpha diversity - READING IN THE DATA (generated in the working with phyloseq code);
data data_alpha;
infile "PATH\alpha_diversity_sas.csv" dlm = "," firstobs = 2;
input Sheep_ID$ pregnancy$ sampling$ bred$ variable$ value;
run;
** This is simply vieweing the data. You will need to close out of the viewer before moving on.;
proc print data = data_alpha;
run;

proc sort data = data_alpha;
by variable;
run;
** similar to the OTU comparison SAS code, you will need to include all variables of interest in the class statement (INCLUDING "variable");
proc mixed data= data_alpha;
by variable;
class INPUT INPUT INPUT INPUT variable;
model value = INPUT INPUT INPUT*INPUT / DDFM=OPTION outp=res; *same here, some options for ddfm= include SATTERTHWAITE, KR and others. Do what is best for your group structure;
*repeated /type = CS subject = INPUT; **repeated variable is set-up slightly different in this code, but you can still just swap in the variable that repeats (see additional repeated structures below);
lsmeans INPUT INPUT INPUT*INPUT /pdiff ADJUST=tukey; **account for multiple testing with a correction best suited for the comparison;
*random Block Rep;
ods output
run; quit;


*repeated /type = CS subject = Sheep_ID;
*repeated /type = UN subject = Sheep_ID;
*repeated /type = AR(1) subject = Sheep_ID;
*repeated /type = TOEP subject = Sheep_ID;
