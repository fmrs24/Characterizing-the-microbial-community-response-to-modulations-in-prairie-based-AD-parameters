##########################################################
############ Working with Phyloseq Master ################
###### Ordinations,phylum and genus level bar graphs #####
###### alpha diversity and OTU comparison bar graphs #####
###### including taxonomy and shared table exports #######
################### LRK 9/3/2021##########################
#                      R 4.0.5                           #
##########################################################

#PhyloSeq is a tool used to analyze amplicon sequence data that integrates
#several R functions and packages. It is great for managing 16S amplicon data,
#and has wonderful support and tutorials available to the user. Basically, 
#anything you want to do with your data, there are forums that have already 
#answered your question.

#In this tutorial, I will be using some tools to quickly summarize your data set.
#We will be generating some numbers and figures that will help you describe the
#data in quick, meaningful ways.

########################################################################### GLOSSARY ##########################################################################
    #[INPUT]-------------- some organization variable that you would like to work with. This variable must be a column title, and represent a column in the metadata file.
                                # often, we would just cluster, remove, color by etc. by the fixed effect of interest. Sometimes, however, we may color points by a "combo" color,
                                # which can help to visualize levels within the fixed effect. This is useful when the data is affected by interaction
    #[LEVEL]-------------- This term specifies a specific level within a fixed effect. This term requires a specific character string or value within a column in the metadata file
    #[FIXEDVAR------------ This stands for any fixed variable that you would like to compare. It is not proper to use the "Combo" input here.
    #[RANDOMVAR]---------- This stands for any random variable that you would like to include.
    #[ == ]--------------- This combination of characters is considered an "operator" and stands for "does equal".
    #[ != ] -------------- This combination of characters is considered an "operator" and stands for "does not equal". 
                                # These operators are throughout the script, and should not need to be scrutinized. The lines they are in (should) all have some form of annotation.
    #[LOWEST_SAMPLE_SIZE]- This will be a numerical value representing the number of sequence reads in the sample with the fewest reads. The idea here is that you will rarefy by
                                # the lowest possible value without excluding that last sample. Often, we sort out samples that are below 5,000 reads, 
                                # then use the next highest value for this input. You can find this information in the Sample_sum_df object created later.


#let me know if you find additional undefined terms!

##########################################################
#################### Installation  #######################
##########################################################

#INSTALLING PACKAGES
#most features in R come from external packages. R is an open source software 
#that anyone can use AND IMPROVE! People make packages for others, but they 
#will first need to be installed. I have included the packages that we need for
#this analysis below.

#IMPORTANT: once installed in the R folder, they will not 
#need to be installed again

#make sure you have installed these programs. If you have run the basic_info
#module, they should already be installed.

#if not, delete the "#" from the lines below, before #LOADING PACKAGES, highlight the lines, and click run

#IMPORTANT: you can watch these commands run in the console window below. IT WILL
#STOP AND ASK TO UPDATE. When prompted, click into the console window and type
# "a" and enter. This will update all packages to the latest edition.

#if (!requireNamespace("BiocManager"))
#install.packages("BiocManager")
#BiocManager::install("phyloseq")
#BiocManager::install("decontam")
#install.packages("ggplot2")
#install.packages("reshape2")
#install.packages("vegan")
#install.packages("dplyr")
#install.packages("splitstackshape")
#install.packages("viridis") #for nice color palettes
#install.packages("remotes")
#install.packages("devtools")
#library(devtools)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
#install.packages("lsmeans")
#install.packages("car")
#install.packages("multcompView")
#install.packages("cld")

#FMS 12.19.22 
#Newer versions of R will have likely lost the command "csplit"
#This package offers an alternative command "concat.split"
#The command has been added to this version of the code as well.
install.packages("splitstackshape")

#hopefully I have included all of the packages that are necessary for this script. Let me know if you run into any errors here.
    #POTENTIAL ERROR: Error in library("UNINSTALLED PACKAGE") : there is no package called 'UNINSTALLED PACKAGE'
##########################################################
################## Loading packages  #####################
##########################################################

#LOADING PACKAGES
#now that you have installed them, you will need to load them to this session. 

#IMPORTANT: unlike installing, packages will need to be loaded every session.

#make sure you have loaded these programs
#to load, highlight all of the library commands and click run

library(phyloseq)
library(ggplot2)
library(vegan)
library(dplyr)
library(scales)
library(grid)
library(reshape2)
library(data.table)
library(stringr)
library(tibble)
library(splitstackshape)
library(viridis)
library(remotes)
library(pairwiseAdonis)
library(nlme)
library(lsmeans)
library(car)
library(multcompView)
library(decontam)
library(splitstackshape)
##########################################################
############# setting up the environment  ################
##########################################################

#the way I structure my working environment
  #[PROJECT FOLDER]
    #|-->[taxonomy]
    #|      |---------> stability.cons.taxonomy
    #|-->[shared]
    #|      |---------> stability.opti_mcc.shared
    #|-->[design]
    #|      |---------> metadata.csv
    #|-->[sas]
    #|      |----------> working_with_phyloseq_OTU_comparison_7_19_21.R
    #|--> working_with_phyloseq_subsetting_7_19_21.R

#It is no problem if you do not structure it like this. You will just need to remember to change the path when trying to read files in.

#if you dont want things in scientific notation, run this (will save changes so will not need to run again)
options(scipen=999)

#Finally, ensure you are reading from the correct directory
# click "session" (above) -> "set working directory" -> "to source file location"
# this will make sure R is reading from the location of the notebook you opened

##########################################################
####### creating a taxonomy object for later use  ########
##########################################################
#reding in the cons.taxonomy file produced by mothur
tax <- read.table("taxonomy/stability.cons.taxonomy", header = TRUE)
#making a list of characters representing the bootstrap values produced in the cons.taxonomy file
bootstraps <- c("(100)", "(99)", "(98)", "(97)", "(96)", "(95)", "(94)", "(93)", "(92)", "(91)", "(90)", "(89)", "(88)", "(87)", "(86)", "(85)", "(84)", "(85)", "(84)", "(83)", "(82)", "(81)", "(80)")
#creating a function that will allow a grep command to escape and important characters.
#we will want to exscape the "()" in the bootstrap values so they will be considered a character and deleted.
regex.escape <- function(string) {
  gsub("([][{}()+*^${|\\\\?])", "\\\\\\1", string)
}
#making a new list including the escapes and a pipe character between each value
bootstrap_input <- paste0("\\b", paste0(regex.escape(bootstraps), collapse="|"), "\\b")
#using gsub to find and replace all values in the input and replace them with nothing
tax$Taxonomy <- gsub(bootstrap_input, "", tax$Taxonomy)
#spliting the file apart into columns by the semicolon
tax <- concat.split(data.frame(tax), "Taxonomy", sep=";", drop=TRUE)
#renaming the columns
names(tax)[3:8] <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
#writing the new dataframe (removal of bootstrap values, splitting by semicolon and adding colum headers)
# to a .csv file in your source directory
#write.csv(tax, "taxonomy.csv", row.names = FALSE)

##########################################################
###### exporting the tax object to an excel table ########
##########################################################
# duplicating the object you created above for modification
tax_w_rel_abund <- tax
#make a list containing the relative abundance of each OTU (the below code is kind of round-a-bout, but it works)
Relative_Abundance <- as.numeric(t((t(tax_w_rel_abund$Size)/rowSums(t(tax_w_rel_abund$Size)))*100))
#append this to the current tax variable to make a new one
tax_w_rel_abund$Relative_Abundance <- Relative_Abundance
#truncate it to the top 100 OTUs
tax_100 <- tax_w_rel_abund[1:100,]
#rename the rows
names(tax_100)[3:8] <- c("Kingdom", "Phylum", "Class", 
                         "Order", "Family", "OTU")

tax_100$Relative_Abundance <- as.numeric(tax_100$Relative_Abundance)

#write.csv(tax_100, file = "relative_abundance_OTU.csv", row.names = FALSE, quote = FALSE)
#you can open .csv files with Excel

##########################################################
##### importing mothur output into a physeq object  ######
##########################################################
## use the tab autocomplete if you do not have shortened file names
data <- import_mothur(mothur_shared_file = "shared/stability.opti_mcc.shared", mothur_constaxonomy_file = "taxonomy/stability.cons.taxonomy")

#read your meta data in as well
map <- read.csv("design/metadata.csv")
########## NOTE
## Consider making your metadata variables into factors at this point with as.factor()
# map$EXAMPLE <- as.factor(map$EXAMPLE)
## particularly, variables which are tracked as numbers, but are not meant to be interpreted as a continuous scale (IE Animal ID number or something)

#convert the metadata into a phyloseq object
map <- sample_data(map)

#naming the rows to correspond with the SampleID column within the metadatafile.
rownames(map) <- map$SampleID

#merging the shared-taxonomy phyloseq object with the metadata file (map) and viewing it
data_merge <- merge_phyloseq(data, map)
data_merge
#assigning taxonomic rankings
colnames(tax_table(data_merge)) <- c("Kingdom", "Phylum", "Class", 
                                     "Order", "Family", "Genus")


#########################################################
# using Decontam to remove potential contaminant sequences
# This code is based on the vignette presented here https://benjjneb.github.io/decontam/vignettes/decontam_intro.html
# two screening methods are used, one which requires DNA concentrations post PCR, and one that requires negative/process controls
# NOTE: to use this package you will want the following 2 metadata columns- DNA concentration (numerical), and Control (y/n)

##########################################################
#####Prior to running decontam,make a pcoa################
###########showing the differences between ###############
############your samples and controls#####################
##########################################################
###change "INPUT" to the title of the column in your
###metadata that distinguishes between controls and samples
set.seed(941996)

control_pcoa <- ordinate(
  physeq = data_merge, 
  method = "PCoA", 
  distance = "bray"
)
#plotting the newly created PCoA
plot_ordination(
  physeq = data_merge,
  ordination = control_pcoa,
  color = "INPUT",
  title = 
) + 
  geom_point(aes(color = INPUT), size = 4) +
  scale_fill_viridis_d(option = "B") 
##########################################################
#DECONTAM PART 1
# FREQUENCY BASED SCREENING - requires DNA concentrations
# Making sequencing depth distribution graph
freq_df <- as.data.frame(sample_data(data_merge)) # Put sample_data into a ggplot-friendly data.frame
freq_df$LibrarySize <- sample_sums(data_merge) # create column with read depth per sample
freq_df <- freq_df[order(freq_df$LibrarySize),] # reorder the dataframe according to read depth
freq_df$Index <- seq(nrow(freq_df)) # create new variable to enforce sample order in the following graph
## In the following graphic replace INPUT with whatever column of your metadata will differentiate between your real and negative/process control samples
## you can also use it to highlight other variables which could be of interest
ggplot(data=freq_df, aes(x=Index, y=LibrarySize, color=INPUT)) + geom_point()
### frequency based screening
## The following command will create a dataframe that indicates whether a given OTU is an contaminant or not
## The version presented here just uses the default settings but you can adjust the threshold used with threshold = _  (p value threshold, default being 0.1)
#
contamdf.freq <- isContaminant(data_merge, method="frequency", conc="INPUT")  #Replace INPUT with column containing concentrations. This line is the one that actually runs the calculations, can take a while
head(contamdf.freq) # displays the first few lines of the table, worth looking at the whole thing
#### making the graphs
table(contamdf.freq$contaminant) # makes a table that displays the number of contaminant OTUs (value returned under TRUE)
head(which(contamdf.freq$contaminant)) #displays the number of the first few OTUs that were flagged as contaminants

# The following command makes graphs which display the frequency of an OTU appearing in samples
# vs the DNA concentration of the sample. Contaminants from kits/PCR should be more common in
# low DNA concentration samples than high (following redline) while real OTUs should be flat (black line)- see vignette for examples
###  IN THE FOLLOWING COMMAND, replace the vector (c(_,_,...)) with whichever OTUs you want to check the graphs for. Likely this will be any OTUs within the top few hundred that were flagged, plus a "normal" OTU for comparison
plot_frequency(data_merge, taxa_names(data_merge)[c(INPUTOTUNUMber,INPUTOTUNUMBER,INPUTOTUNUMBER)], conc="INPUT") +  #Replace INPUT with column containing concentrations,
  xlab("DNA Concentration (PicoGreen fluorescent intensity)")

### The following code can be used if you would like to do random spot checks instead of picking OTUs to graph and observe- not always necessary, but sometimes informative

#set.seed(941996) #using same seed as below for consistancy/reproducability
#plot_frequency(data_merge, taxa_names(data_merge)[sample(which(contamdf.freq$contaminant),3)], conc="INPUT") + # this is pulling 3 random OTUs from the dataframe to plot-
# xlab("DNA Concentration (Insert Measurement Technique)") # note that it may pull rare OTUs that only have 10 reads and may not appear in most of your samples

#

######### Removing based on frequency- THIS IS THE PART THAT ACTUALLY REMOVES FLAGGED OTUS FROM YOUR DATA. #####
### Note that if you want to keep any OTUs that were flagged (based on observation of the graph, or knowledge of the system (IE. you would expect to observe that taxa in that environment))
### you will have to manually change the TRUE to FALSE  for that OTU (In a later update I may provide code to do so but for now you can do it in R, or export and re-import a csv)
# write.csv(contamdf.freq, file = "contamfreqdf.csv", row.names = TRUE, quote = FALSE)
data.noncontamF <- prune_taxa(!contamdf.freq$contaminant, data_merge) # here we are creating a new phyloseq object to carry forward, preserving the original data_merge
data.noncontamF # double check that it looks right


##############################
# DECONTAM Part 2:
# Prevalence based screening - requires negative controls
##############################
### preparing to check prevalence by taking the phyloseq object from part 1 and reducing OTUs abundance to binary
### if you did not run Part 1 for whatever reason you can run the following line -
# data.noncontamF <- data_merge # DO NOT RUN THIS IF YOU DID THE FIRST PART it will overwrite and undo it
### the basic principle of prevalence based screening is checking the ratio of how often an OTU is present in your controls vs your real samples
sample_data(data.noncontamF)$is.neg <- sample_data(data.noncontamF)$INPUT == "Control" ## replace INPUT and Control with the column and entry denoting negative controls
contamdf.prev <- isContaminant(data.noncontamF, method="prevalence", neg="is.neg") # this is the line that actually does the math, may take a while
table(contamdf.prev$contaminant) # again, displays number of flagged contaminant OTUs under "TRUE"
head(which(contamdf.prev$contaminant))

## ^^^ Following the above calculation I again recommend looking at the contamdf.prev sheet yourself and checking the ID of any OTUs in the first few hundred to see if they are things you would expect in your data.
#
### Now, even though we already ran Decontam and have its results, we will make sheets to enable us to look at graphs and make sense of it ourselves

ps.pa <- transform_sample_counts(data.noncontamF, function(abund) 1*(abund>0)) ## creating a sheet with binary 1/0 presence/absence of OTUs in samples
ps.pa.neg <- prune_samples(sample_data(ps.pa)$INPUT == "Control", ps.pa) # sheet for negative controls only 
ps.pa.pos <- prune_samples(sample_data(ps.pa)$INPUT != "Control", ps.pa) # sheet for real samples only
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant) # Creates a dataframe that contains the number of times an OTU appears in real samples and in negative controls respectively.
## Each point on the following graph represents an OTU and shows how often it appears in real samples (Y axis) vs negative controls (X axis)
## it will also color them based on whether or not Decontam flagged them as a contaminant based on prevalence. However if you want to see what OTUs were flagged, check contamdf.prev

ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() + 
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

##### again, if based on your own expertise/the specifics of your data you would like to keep some of the OTUs it flagged as potential contaminants, you will need to change their entry in the Contaminant column
# write.csv(contamdf.prev, file = "contamprevdf.csv", row.names = TRUE, quote = FALSE)
#

#### Removal based on prevalence - THIS IS THE PART THAT ACTUALLY REMOVES OTUs FROM YOUR DATA for prevalence ####
####### REPLACE CONTAMDF.PREV with contamprevdf below if you exported and reimported the csv as above to edit any of the flagged OTUs
data.noncontamPF <- prune_taxa(!contamdf.prev$contaminant, data.noncontamF) ##removing the OTUs and making the new phyloseq object
data.noncontamPF ## double check that it looks correct

#########################################################
#################### POST DECONTAM ######################

### if you ran Decontam you will need to choose one of the following options to continue with the code
### if you did not run it, ignore this

### Option 1- edit the first two commands in the next section to call on the phyloseq object you just made instead of data_merge
# taxa_sum_df <- data.frame(sum = taxa_sums(data.noncontamPF))
# data_sub <- data.noncontamPF %>%
# prune_taxa(taxa_sums(.) > 9, .)
### Pro- maintains original data_merge object in case you want it later, Con- takes up more space

### Option 2- overwrite data_merge by running the following
 #data_merge <- data.noncontamPF
### Pro- saves space, Con- original object is gone if you need it for some reason

### after that you may proceed with the rest of the code as normal, since your data_sub object will have already removed the contaminant OTUs


##########################################################
#### Removing OTUs that represent less than 10 reads  ####
##########################################################

#checking the number of OTUs before pruning
taxa_sum_df <- data.frame(sum = taxa_sums(data_merge))

#pruning the taxa so only OTUs with more than 10 sequences per read remain
data_sub <- data_merge %>%
  prune_taxa(taxa_sums(.) > 9, .)

#checking the number of OTUs left after pruning
pruned_taxa_sum_df <- data.frame(sum = taxa_sums(data_sub))

#finding the sum of sequences in filtered dataset
total_trimmed_seqs <- sum(pruned_taxa_sum_df)

#changing and checking the taxonomy names given
colnames(tax_table(data_sub)) <- c("Kingdom", "Phylum", "Class", 
                                   "Order", "Family", "Genus")
colnames(tax_table(data_sub))

##########################################################
############# basic info about each sample  ##############
##########################################################

#listing the sequence number per sample
sample_sum_df <- data.frame(sum = sample_sums(data_sub))
#Click the sample_sum_df table in your environment to view it, and sort by sum to ID the lowest values- 
# this will let you discern whether to remove samples with low read depth, and tell you what number to subsample to in later sections

#get the average sequence length and standard deviation
average_sequencing_depth <- mean(sample_sum_df$sum)
standard_deviation_of_sequencing_depth <- sd(sample_sum_df$sum)

#plotting the numbers generated above- consider changing bin width if the resolution is too low for your data.
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())

#write.csv(sample_sum_df, file = "Sample_sum_df.csv", row.names = FALSE, quote = FALSE)

##########################################################
####### Removal of samples based on metadata type ########
##########################################################

#I am removing samples based on a level from INPUT variable (column header).
# example: remove the LEVEL "treatment1" from INPUT Trt_Type

#removing samples from the working phyloseq object
  #removal of all samples within a specific level of a variable
data_sub_type <- subset_samples(data_sub, INPUT != "LEVEL")
  #removal of a specific sample from the SampleID column (removal of a single sample)
data_sub_type <- subset_samples(data_sub_type, SampleID != "LEVEL")

#### Removal of samples based on a set sequencing depth threshold
sub_df <- as.data.frame(sample_data(data_sub)) # substitute data_sub for data_sub_type if you already removed samples
sub_df$LibrarySize <- sample_sums(data_sub) # create column with read depth per sample
sub_df <- sub_df[order(sub_df$LibrarySize),] ## sorting by library size for each of use so you can look at the table and make decisions. sample_sum_df would be very similar but without the metadata
#### The line below sets the threshold for removal
sub_df$Keep <-  1*(sub_df$LibrarySize>5000) ### Replace the 5000 with whatever threshold is appropriate for your data! Our lab considers 5000 to be the absolute minimum
### We often use 10000 or even higher as a minimum however, depending on the distribution of the data in question.
sub_df <- sub_df[order(sub_df$SampleID)] # sort by SampleID
data_sub_type <- subset_samples(data_sub, sub_df$Keep != 0) ## this is the line that actually removes the samples- AGAIN, REPLACE WITH data_sub_type IF YOU ALREADY REMOVED SAMPLES
data_sub_type ## checking that smaple number matches expected number

#removing the same values from the metadata file
map_sub_type <- map[map$INPUT != "LEVEL",]
map_sub_type <- map_sub_type[map_sub_type$SampleID != "LEVEL",]
## removing by depth threshold again- NOTE: SUB map_sub_type in for map if you have already removed other samples.
map <- map[order(map$SampleID)] ## sorting by SampleID to make sure it matches sub_df's order- only works if any other removals match
map_sub_type <- subset_samples(map, sub_df$Keep != 0) ## removes samples by referencing sub_df- therefore it is critical that they are in the same order and have the same samples already removed

#making and combo colors to color by later (optional)
#map_df$combo <- with(map_df, paste0(INPUT2, "_", INPUT2))

#making a map dataframe for merging downstream
map_df <- read.csv("design/metadata.csv")
map_df
#removing the same values from the metadata dataframe
map_sub_type_df <- map_df[map_df$INPUT != "LEVEL",]
map_sub_type_df <- map_sub_type_df[map_sub_type_df$SampleID != "LEVEL",]
### or doing it by depth threshold- again sub map_sub_type_df for map_df in if you already removed by level
sub_dfm <- as.data.frame(sub_df$SampleID)
colnames(sub_dfm)[1]="SampleID"
sub_dfm$Keep <- sub_df$Keep
map_dfm <- merge(map_df, sub_dfm, by = "SampleID")
#map_dfm <- map_dfm[order(map_dfm$SampleID)] ## just to make sure the order matches between dataframes
map_sub_type_df <- map_dfm[map_dfm$Keep != 0,]
#making and combo colors to color by later (optional)
#map_df$combo <- with(map_df, paste0(INPUT2, "_", INPUT2))

#after removal of samples, it is likely that some OTUs will drop below 10 sequences
#in order to get an accurate count of OTUs left after subsetting per group, you can run the following commands

#Similar to above, removing OTUs represented by less than 10 sequences
data_sub_type_prune <- data_sub_type %>%
  prune_taxa(taxa_sums(.) > 9, .)
#creating a new taxonomy table with sequence counts
pruned_taxa_type_sum_df <- data.frame(sum = taxa_sums(data_sub_type_prune))
#calculating the sum of sequences per sample
sample_sum_df_prune <- data.frame(sum = sample_sums(data_sub_type_prune))
#average number of sequences per sample after pruning
average_sequencing_depth_prune <- mean(sample_sum_df_prune$sum)
#std.of sequences per sample after pruning
standard_deviation_of_sequencing_depth_prune <- sd(sample_sum_df_prune$sum)
#calculating the sum of sequences left after pruning
total_trimmed_seqs_pruned <- sum(pruned_taxa_type_sum_df)
taxasum <- sum(pruned_taxa_type_sum_df$sum)
#calculating new relative abundance of taxa for specific group
pruned_taxa_type_sum_df$relabund <- ((pruned_taxa_type_sum_df$sum)/taxasum)*100
#Changing the object to a data table to maintain rownames
pruned_taxa_type_sum_dt <- setDT(pruned_taxa_type_sum_df, keep.rownames = TRUE)[]
#renaming the first column "OTU" for merging
names(pruned_taxa_type_sum_dt)[1]<-paste("OTU")
#merging the taxonomy information with the taxonomy sum information
pruned_taxa_type_sum_dt <- merge(x = pruned_taxa_type_sum_dt, y = tax[ , c("OTU", "Genus")], by=c("OTU"))
#ordering by abundance - sometimes when subsetting, OTU order changes
pruned_taxa_type_sum_dt_ord <- pruned_taxa_type_sum_dt[order(-relabund),]
##writing this subset tax_w_relabund file to a csv (change INPUT to whatever specific group you are analyzing)
write.csv(pruned_taxa_type_sum_dt_ord, file = "relative_abundance_OTU_S2only.csv", row.names = FALSE, quote = FALSE)


##########################################################
############## Alpha diversity measures  #################
##########################################################

#begin by setting a repeatable seed number for random subsetting of the data.
set.seed(941996)


#subsample the data to lowest sample sequence depth (See sample_sum_df for LOWEST_SAMPLE_SIZE)
data_scaled = rarefy_even_depth(data_sub_type, sample.size = LOWEST_SAMPLE_SIZE)

#Alternatively, you can choose not to subsample (there is argument within the field as to which is right)
#data_scaled = data_sub_type

#making objects that are character strings for each test (can be input below)
alpha_meas = c("Observed", "Chao1", "Shannon", "Simpson", "InvSimpson")
#chao1 is a species richness estimator. tries to estimate total number of species independent of sampling depth
chao1 <- "Chao1"
#shannon is a diversity estimator which provides a value that is based on a combination of both species richness and evenness
shannon <- "Shannon"
#sobs is total species observed
sobs <- "Observed"
#simpson is a community evenness calculation. the higher the number, the more even the communities (comparison of relative abundance for all diff species)
simpson <-"Simpson"
#adding in the subset metadata file if needed

#using the estimate richness command to get measures for each test above (included in alpha_meas object)
data_alpha <- estimate_richness(data_scaled, split = TRUE, measures = alpha_meas)
#making a SampleID column to merge with (retains row names when setting as a datatable as column 1)
data_alpha <- setDT(data_alpha, keep.rownames = TRUE)[]
#renaming column 1
names(data_alpha)[1]<-paste("SampleID")
#if your samples are integers, this code automatically adds a capital X to the front. This code will remove the first character to allow accurate merging (DO NOT RUN MULTIPLE TIMES)
#data_alpha$SampleID <- gsub("^.", '', data_alpha$SampleID)

#combining the metadata to the alpha diversity measures output object
data_alpha <- merge(x = data_alpha, y = map_sub_type_df, by=c("SampleID"))

#Write this object to a csv file. Use the entries here for comparing means/ making statistical comparisons
write.csv(data_alpha, file = "sample_diversity_estimators.csv", row.names = FALSE, quote = FALSE)

#Making boxplots for alpha diversity measures
#estimating richness for all samples NOTE: we will not be aggregating this time.
data_alpha_full <- estimate_richness(data_scaled, split = TRUE, measures = alpha_meas)
#OPTIONAL: removing some measurments based on what you want to see in your graphic
data_alpha_full$se.chao1 <- NULL
data_alpha_full$InvSimpson <- NULL

#setting it to a data table and retaining the rownames as the first column
data_alpha_full <- setDT(data_alpha_full, keep.rownames = TRUE)[]
#renaming the new first column to SampleID
names(data_alpha_full)[1]<-paste("SampleID")
#if your samples are integers, this code automatically adds a capital X to the front. This code will remove the first character to allow accurate merging (DO NOT RUN MULTIPLE TIMES).
#data_alpha_full$SampleID <- gsub("^.", '', data_alpha_full$SampleID)

#Merging the INPUT column with your alpha diversity measurements
  #NOTE: The order should match the model you plan to use in SAS for all of the following steps.
data_alpha_full <- merge(x = data_alpha_full, y = map_sub_type_df[ , c("SampleID", "INPUT", "INPUT", "INPUT")], by=c("SampleID"))  
#melting the data based on both the SampleID and INPUT variable
molten_alpha_meas_full <- reshape2::melt(data = data_alpha_full, id.vars = c("SampleID", "INPUT", "INPUT", "INPUT"))

#exporting a csv file to use with the added alpha diversity SAS code
#merging the alpha diversity measurements with 
data_alpha_merged <- merge(x = data_alpha[ , c("SampleID", "Observed", "Chao1", "Shannon", "Simpson")], 
                           y = map_sub_type_df[ , c("SampleID", "INPUT", "INPUT", "INPUT", "INPUT", "INPUT")], by=c("SampleID"))
#removing the sampleID
data_alpha_merged$SampleID <- NULL
#melting based on important variables
data_alpha_merged_molten <- reshape2::melt(data = data_alpha_merged, id.vars = c("INPUT", "INPUT", "INPUT", "INPUT"))
#writing the df to a file (can be directly imported into SAS using the code)
write.csv(data_alpha_merged_molten, file = "alpha_diversity_sas.csv", row.names = FALSE, quote = FALSE)

#This code will create a histogram for each of your measurements and facet by variables of interest
  #NOTE: the averages presented in this graphic are not representative of the LSMeans created in SAS (but this can still be used to display info with significance.)
ggplot(data = molten_alpha_meas_full, aes(y=value, x=INPUT, color=INPUT)) +
  geom_boxplot() +
  # for facet wrapping, you can wrap by multiple variables. If you would like to simply wrap by alpha diversity measures remove the + INPUT
  facet_wrap((~variable), scales="free") +
  #scale_fill_viridis_d(option = "B") + #this is an alternative to the above scale fill
  #viridis is really nice for generating colorblind friendly figures with nice separation of color
  #to add this option, remove the # before "scale_fill_viridis" and a # before "scale_fill_manual"
  scale_color_manual(values= c("goldenrod", "firebrick", "aquamarine3", "darkorchid3")) +
  geom_jitter(shape=16, position=position_jitter(0.2))

#Alternatively, you can run the following code and never leave R!
#Some cool code from Laura Tibbs-Cortes: statistical testing for alpha diversity measures!!!

#Begin by building a model
alphadiv.lme.model <- lme(data=molten_alpha_meas_full %>% 
  #filter by the alpha diversity metric you want to test (WILL HAVE TO RUN THIS MULTIPLE TIMES)                         
  filter(variable=="VARIABLE"), 
  #including the fixed effect that you wish to test for (including interactions) 
  value ~ FIXEDVAR1 + FIXEDVAR2 + FIXEDVAR1*FIXEDVAR2,
  #including any random or repeated variables
  #option = 
  random=~1|RANDOMVAR) 
#Testing the separation in means for each variable specified above, correcting for any random variables
Anova(alphadiv.lme.model, type="III")
#getting pariwise estimates of the LSMeans
lsmeans(alphadiv.lme.model, list(pairwise ~ FIXEDVAR1|FIXEDVAR2)) # to get difference of Least Squares Means, Tukey-adjusted for multiple comparisons
#to get Least Squares Means of fixed effects (should only use if no interactions are significant)
  #lsmeans(alphadiv.lme.model, list(~FIXEDVAR1, ~FIXEDVAR2))
alpha_div_input <- lsmeans(alphadiv.lme.model, list(pairwise ~ FIXEDVAR1|FIXEDVAR2))

#writing the results to an object then extracting them as a data frame
alpha_div_results <- lsmeans(alphadiv.lme.model, list(pairwise ~ FIXEDVAR1|FIXEDVAR2))
alpha_div_input <- as.data.frame(summary(alpha_div_results$`lsmeans of FIXEDVAR1 | FIXEDVAR2`))
#adding a column to specify the levels
alpha_div_input$combo <- with(alpha_div_input, paste0(FIXEDVAR1, "_", FIXEDVAR2))

#Plotting the output
ggplot(alpha_div_input, aes(x = combo, y = lsmean, fill = combo)) + 
  #this specifies that you would like to use a bar graph that has black outlines
  geom_bar(stat = "identity", colour = "black") +
  #this option includes the sampleIDs along the x-Axis
  theme(axis.text.x = element_text(vjust=0.5, angle = 90)) +
  scale_fill_manual(values= c("goldenrod", "aquamarine3", "firebrick", "darkorchid3")) +
  #using the upper and lower confidence intervals produced by the lsmeans command as error bars
  geom_errorbar(aes(ymin=alpha_div_input$lower.CL, 
                    ymax=alpha_div_input$upper.CL), width=.2) +
  #add a title that matches the variable of interest
  xlab("VARIABLE")

#see example below for more details

##########################################################
#################### Ordinations  ########################
##########################################################

#now I'm going to create the PCoA and CCA plots, visualizing those community clusters

#here we are randomly subsampling the sequences based on a given seed. By providing the seed, we ensure that the same sequences are subsampled each time. This helps reproducibility. 
set.seed(941996)
#change the sample.size input number to the lowest sequences per sample value. (can be found in sample_sum_df)
data_scaled = rarefy_even_depth(data_sub_type, sample.size = LOWEST_SAMPLE_SIZE)

#Alternatively, you can choose not to subsample (there is argument in the field as to which is right)
#data_scaled = data_sub_type


##################### UNCONSTRAINED ORDINATIONS ########################
#This is used to visualize the largest amount of variation across all variables
#this tutorial uses PCoA, but can also use NMDS by changing method = NMDS

#Here we create the coordinated matrix for a PCoA plot using in built commands from phyloseq
data_pcoa <- ordinate(
  physeq = data_scaled, 
  method = "PCoA", 
  distance = "bray"
)

#if needed
  #data_scaled$Sample <- as.character(data_scaled$INPUT)

#plotting the newly created PCoA
plot_ordination(
  physeq = data_scaled,
  ordination = data_pcoa,
  color = "INPUT",
  # shape = "INPUT2",
  title = 
) + 
  geom_point(aes(color = INPUT), size = 4) +
  scale_color_manual(values= c("goldenrod", "aquamarine3", "firebrick", "darkorchid3"))#+
  #stat_ellipse(aes(group=INPUT)) #+
  #scale_shape_manual(values=c(0,1,2,5)) #Add shape = INPUT2 to geom_point above-Use this if you want to dictate the shapes used, check this list for values http://www.sthda.com/english/wiki/ggplot2-point-shapes

#NOTE the colors in "scale_color_manual" can be replaced by NA if you want to keep the structure but only display a subset of points
#ex: #scale_color_manual(values= c("NA", "NA", "NA", "NA", "black", "salmon", "tan", "Blue")) +

##################### CONSTRAINED ORDINATIONS ########################
#This is used to visualize the largest amount of variation for specific variables of interest
  #this module has a built-in portion to test correlations between numerical values, helping you eliminate redundant variables
  #additionally, it is also able to test the significance of all variables, but I still prefer to use the PERMANOVA module below for this purpose

#################### REFERENCE TEXT ####################
#https://esajournals.onlinelibrary.wiley.com/doi/10.1890/0012-9658%282003%29084%5B0511%3ACAOPCA%5D2.0.CO%3B2
#http://deneflab.github.io/MicrobeMiseq/demos/mothur_2_phyloseq.html
#https://www.gdc-docs.ethz.ch/MDA/handouts/MDA20_PhyloseqFormation_Mahendra_Mariadassou.pdf

#################### making a data frame based on the data above ######################## 
data_sub_df <- as(sample_data(data_sub_type), "data.frame")

#################### identify correlated variables #################################
#identify correlated variables and remove one or both from the analysis - reduce noise
####################################################################################

#this ony works with numerical values. I think it is still OK to run categorical variables below

# reading in the metadata again
map_cor <- read.csv("design/metadata.csv")
map_cor_sub_type_df <- map_cor[map_cor$INPUT != "LEVEL",]
#map_cor_sub_type_df <- map_cor_sub_type_df[map_cor_sub_type_df$INPUT != "LEVEL",]

#setting rownames to sample IDs
rownames(map_cor_sub_type_df) <- map_cor_sub_type_df$SampleID
#removing sample IDs column
map_cor_sub_type_df = subset(map_cor_sub_type_df, select = -c(SampleID) )
#checking the structure of each column (treatment and dog are character strings)
str(map_cor_sub_type_df)
#making a list of columns that are numeric or not (including integer)
num_columns <- unlist(lapply(map_cor_sub_type_df, is.numeric))
#removing character columns (treatment and dog)
map_cor_numeric <- map_cor_sub_type_df[ , num_columns]
#checking the correlations between numeric variables (cor only works on numerical variables, possibly use dummy variables?)
correlated_variables <- as.data.frame(cor(map_cor_numeric, method="pearson"))
#same as above but giving a different name
correlated_variables_edit <- as.data.frame(cor(map_cor_numeric, method="pearson"))
#creating a data frame with only the correlations about 0.5 (absolute value)
correlated_variables_edit[ abs(correlated_variables_edit)<0.5 ] <- ""

#################### calculate the distance matrix ########################
#calculate distance matrix between samples based on bray-curtis dissimilarity
###########################################################################

#calculating the bray-curtis dissimilarity between each sample
dist.bc <- distance(data_sub_type, "bray")

#################### calculating the eigen value and inertia for the cap variables ########################
#Here we used some variables to constrain by that are able to detect differences in animal treatment groups
###########################################################################################################

#I selected only a few variables based on correlation and relevance to the question (fat digestibility). 
#So, If the variable was not correlated to anything, I included it in the model.
#If there was correlation between variables, I selected the most relevant variable to represent the rest.

#I created the cap scaling object here including the selected variables and the dataframe used above.
cap <- capscale(dist.bc ~ INPUT + INPUT + INPUT*INPUT, data = data_sub_df)
#example
#cap <- capscale(dist.bc ~ Period + Room + Treatment + Fat_Digestibility + Fecal_DM + Feed_Intake_Asfed + RDW, data = data_sub_df)
#This displays the cap loadings that will be used, as well as the inertia (measure of effect size?)
cap
#and running an anova to detect significance of the variables influencing the data.
#in this case, it seems Period and Fat_Digestibility does influence the data.
anova(cap, by = "term")

#################### creating an ordination based on these variables ########################
#Here we create the coordinated matrix for a CAP plot using in built commands from phyloseq
#############################################################################################

#EXAMPLE: PLEASE ADJUST THE FORMULA TO MATCH YOUR QUESTION

#creating a repeatable color scheme
Treatment_colors <- c("darkgreen", "darkorange2", "blue", "deeppink2")

#Using the phyloseq package, we can create ordination plot loadings based on CAP calculations
#I used the same variables as above.
data_cap <- ordinate(
  physeq = data_sub_type, 
  method = "CAP", 
  distance = "bray",
  formula = ~ INPUT + INPUT + INPUT*INPUT
  #example
  #formula = ~ Period + Room + Treatment + Fat_Digestibility + Fecal_DM + Feed_Intake_Asfed + RDW
)

#this code is for defining the arrow maps that are over-layed on the bi-plot
arrowmat <- vegan::scores(data_cap, display = "bp")
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)
arrow_map <- aes(xend = CAP1, 
                 yend = CAP2, 
                 x = 0, 
                 y = 0, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

label_map <- aes(x = 1.3 * CAP1, 
                 y = 1.3 * CAP2, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

arrowhead = arrow(length = unit(0.02, "npc"))

#plotting the newly created CAP
cap_plot <- plot_ordination(
  physeq = data_sub,
  ordination = data_cap,
  color = "INPUT",
  title = , 
  axes = c(1,2)
) + 
  geom_point(aes(colour = INPUT), size = 4) + 
  geom_point(size = 1.5) + 
  scale_color_manual(values = c("goldenrod", "aquamarine3", "firebrick", "darkorchid3"))

cap_plot

#this section add arrows that demonstrate the effect size of the variable of interest
  #NOTE: they are a bit difficult to explain in a manuscript, so I usually leave them off
cap_plot + 
  geom_segment(
    mapping = arrow_map, 
    size = .5, 
    data = arrowdf, 
    color = "Black", 
    arrow = arrowhead
  ) + 
  geom_text(
    mapping = label_map, 
    size = 4,  
    data = arrowdf, 
    show.legend = FALSE
  )

##########################################################
################ Phylum level information  ###############
##########################################################

#begin by grouping all the data into phylum level
data_phylum <- data_sub_type %>%
  tax_glom(taxrank = "Phylum") %>%                     
  # agglomerate at phylum level 
  transform_sample_counts(function(x) {(x/sum(x))} ) %>% 
  # Transform to rel. abundance
  psmelt() %>%                                         
  # Melt to long format
  arrange(Phylum)                                      
  # Sort data frame alphabetically by phylum

#shows you how many different phyla are as well as the number of samples that contain it
number_of_phyla <- as.data.frame(data_phylum %>% group_by(Phylum) %>% summarize(count=n()))

#I have to subset the the data set one more time
data_phylum_only <- data_sub_type %>%
  tax_glom(taxrank = "Phylum")

#I also want to see the number of of sequences within each phylum and their % abundance
#begin by calculating the sum of sequences for each phylum 
pruned_phylum_sum_df <- data.frame(sum = taxa_sums(data_phylum_only))
#make a list containing only the unique phylum level entries
phyla_list <- as.list(get_taxa_unique(data_phylum_only, taxonomic.rank=rank_names(data_phylum_only)[2], errorIfNULL=TRUE))
#make those phyla the row names
row.names(pruned_phylum_sum_df) <- phyla_list
#adding a row with the relative abundances of each phyla
pruned_phylum_sum_df$relative_abundance <- (pruned_phylum_sum_df$sum/sum(pruned_phylum_sum_df$sum)*100)
#ordering the phyla by relative abundance
pruned_phylum_sum_df <- pruned_phylum_sum_df[order(-pruned_phylum_sum_df$relative_abundance),]
#setting the dataframe to a data table for exporting via write.xlsx
pruned_phylum_sum_df <- setDT(pruned_phylum_sum_df, keep.rownames = TRUE)[]

write.csv(pruned_phylum_sum_df, file = "relative_abundance_phylum.csv", row.names = TRUE, quote = FALSE)

#changing the name of the first column to phyla
names(pruned_phylum_sum_df)[1] <- "phyla"
#removing any unclassified bacteria
pruned_phylum_sum_df2 <- pruned_phylum_sum_df[pruned_phylum_sum_df$phyla!= "Bacteria_unclassified"]
pruned_phylum_sum_df2 <- pruned_phylum_sum_df2[pruned_phylum_sum_df2$phyla!= "unknown_unclassified"]
#making a list of the 10 most abundant phyla names
top_10_phylum <- as.list(pruned_phylum_sum_df2$phyla[1:10])

#making a matrix from the phyloseq object (maybe could have went straight with data.frame here)
phylum_shared <- as(otu_table(data_phylum_only), "matrix")
#converting it from a matrix to a data frame
phylum_shared = as.data.frame(phylum_shared)
#making the row names of this object the phyla listed above
row.names(phylum_shared) <- phyla_list

#creating a phylum level shared file with the added phyla names
write.csv(phylum_shared, file = "phylum_shared.csv", row.names = TRUE, quote = FALSE)

##########################################################
############## Phylum level visualizations  ##############
##########################################################

######################################################################
### This set of visualizations show the composition of each sample ###
######################################################################

#making that plot!
#first I need to convert the animal column to a categorical variable using as.charater()
#maybe not needed (use only for categorical variables with levels that are numbers (animal tag numbers))
#data_phylum$INPUT <- as.character(data_phylum_subset$INPUT)
###check SampleID
#now for the plot
ggplot(data_phylum, aes(x = Sample, y = Abundance, fill = Phylum)) + 
  #this specifies that you would like to use a bar graph that has black outlines
  geom_bar(stat = "identity", colour = "black") +
  
  # Remove x axis title and add a text angle to the x ticks. Additionally, I added a general text change to sans type fonts
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90)) + 
  
  scale_fill_viridis_d(option = "B") + #this is an alternative color palette tool
  #viridis is really nice for generating colorblind friendly figures with nice separation of color
  #to add this option, remove the # before "scale_fill_viridis"
  
  guides(fill = guide_legend(keywidth = 1, keyheight = 1)) +
  
  ylab("Relative Abundance (%) \n")


###########################################################
######### making that plot faceted by a variable ##########
###########################################################

#truncating the phyla represented in the graph
data_phylum_short <- data_phylum[data_phylum$Phylum %in% top_10_phylum,]

#if needed, yu can change the animal number to a character value so it is not considered continuous
#data_phylum_short$Animal <- as.character(data_phylum_short$SampleID)

ggplot(data_phylum_short, aes(x = SampleID, y = Abundance, fill = Phylum)) + 
  
  #change the input here to which variable you want to facet by
  #facet_grid(INPUT~.) +
  
  #to remove blank space, use the facet_wrap option instead of the facet_grid
  facet_wrap(INPUT~., scales = 'free', nrow = 1) +
  
  #this specifies that you would like to use a bar graph that has black outlines
  geom_bar(stat = "identity", colour = "black") +
  
  # Remove x axis title and add a text angle to the x ticks. Additionally, I added a general text change to sans type fonts
 
  #this option includes the sampleIDs along the x-Axis
  theme(axis.title.x = element_blank(), axis.text.x = element_text(vjust=0.5, angle = 90)) + 
  
  #this option excludes the sampleIDS from the x-axis
  #theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x=element_blank()) +

  scale_fill_viridis_d(option = "B") + #this is an alternative color palette tool
  #viridis is really nice for generating colorblind friendly figures with nice separation of color
 
  
  guides(fill = guide_legend(keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (%) \n")

#######################################################################
#### This set of visualizations show the data merged by a variable ####
############## and subset to the most abundant phyla ##################
#######################################################################
#NOTE: list of most abundant phyla was made above

#creating a palette of color to be used in your plot. This can be changed.
#newpalette <- c("#F0E442", "#000000", "#990000", "#CC79A7", "#E69F00", "#56B4E9", "#D55E00", "#009E73", "#CCCC33", "#0072B2")
#I would like to only view bar graphs containing phylum with higher that 1% abundance
#to do this, I will try to subset the 

#The next step is for producing bar graphs comparing phyla
data_phylum_subset <- data_sub_type %>%
  # agglomerate at phylum level
  tax_glom(taxrank = "Phylum") %>%                     
  #convert to proportions to avoid weighing samples by their read depth
  transform_sample_counts(function(x) {(x/sum(x))} ) %>% 
  #The merge_samples function sums the abundances of samples with the sample SampleType
  #if you don't want to merge by a variable, you can hash the next to code lines out
  merge_samples("INPUT") %>%
  #Then, we need to convert to proportions again, since the total abundance of each SampleType will equal the number of samples that were merged
  transform_sample_counts(function(x) {(x/sum(x))*100} ) %>%
  # Melt to long format
  psmelt() %>%                                         
  #Sort data frame alphabetically by phylum
  arrange(Phylum)                                      

###We needed to filter the dataset based on the top abundance FILE MADE ABOVE###
data_phylum_subset <- data_phylum_subset[data_phylum_subset$Phylum %in% top_10_phylum,]

#making that plot!
#first I need to convert the animal column to a categorical variable using as.charater()
data_phylum_subset$Sample <- as.character(data_phylum_subset$Sample)
#then plot (X = "Sample")
#if you look at the column that you merged by (see above) you will notice that column name has been changed to Sample. not sure Why it does this.
ggplot(data_phylum_subset, aes(x = Sample, y = Abundance, fill = Phylum)) + 
  #this specifies that you would like to use a bar graph that has black outlines
  geom_bar(stat = "identity", colour = "black") +
  
  #scale_fill_manual(values=newpalette) +
  
  scale_fill_viridis_d(option = "B") + #this is an alternative to the above scale fill
  #viridis is really nice for generating colorblind friendly figures with nice separation of color
  #to add this option, remove the # before "scale_fill_viridis" and a # before "scale_fill_manual"
  
  guides(fill = guide_legend(keywidth = 1, keyheight = 1)) +
  
  ylab("Relative Abundance (%) \n") +
  
  xlab("INPUT")

####################################################################
############ Alternative way of viewing the phylum data ############
####################################################################

ggplot(data_phylum_subset, aes(x = reorder(Phylum, Abundance), y = Abundance, fill = Sample)) + 
  #here I am denoting I would like to use a bar graph, but instead of a stacked bar, I want a side by side bar (position dodge)
  geom_bar(stat = "identity", position =position_dodge2(reverse = TRUE)) +
  #Telling ggplot I would like to use these two colors to denote each cheese type
  scale_fill_manual(values= c("goldenrod", "aquamarine3", "firebrick", "darkorchid3")) +
  #labeling the Y axis
  ylab("Relative Abundance") +
  #labeling the X axis
  xlab("Phyla") +
  #flipping the axis
  coord_flip()

#################################################
############ Genus level information ############
#################################################

#IF YOU ARE HAVING DIFFICULTIES WITH YOUR GENUS LEVEL GRAPHICS DUE TO REPEATED GENUS NAMES
#SEE EXPERIMENTAL OTUID MODULE BELOW

#I have to subset the the data set one more time
data_genus_only <- data_sub_type %>%
  tax_glom(taxrank = "Genus")

#I also want to see the number of of sequences within each phylum and their % abundance
pruned_genus_sum_df <- data.frame(sum = taxa_sums(data_genus_only))

pruned_genus_sum_df <- setDT(pruned_genus_sum_df, keep.rownames = TRUE)[]
#naming the new row "OTU"
names(pruned_genus_sum_df)[1] <- "OTU"

merged_genus_list <- merge(x = pruned_genus_sum_df, y = tax[ , c("OTU", "Genus")], by = "OTU", all.x=TRUE)
#calculating relative abundance (currently of the subset data)
merged_genus_list$relative_abundance <- (merged_genus_list$sum/sum(merged_genus_list$sum)*100)
#ordering the genuses based on relative abundance
merged_genus_list <- merged_genus_list[order(-merged_genus_list$relative_abundance),]
merged_genus_list <- data.frame(lapply(merged_genus_list, as.character), stringsAsFactors=FALSE)

#write to csv
write.csv(merged_genus_list, file = "relative_abundance_genus.csv", row.names = TRUE, quote = FALSE)

#make a top 10 list
top_10_genus <- as.list(merged_genus_list$Genus[1:10])

#alternatively, make a top 20 list
top_20_genus <- as.list(merged_genus_list$Genus[1:20])

#### MAKE GENUS SHARED!!!
data_genus_SH <- data_sub_type %>%
  tax_glom(taxrank = "Genus")
genus_shared <- as(otu_table(data_genus_SH), "matrix")
#converting it from a matrix to a data frame
genus_shared = as.data.frame(genus_shared)
#move OTU labels to column
genus_shared$OTUList <- row.names(genus_shared)
#genus_shared_IDs <- merge(x = genus_shared , y = tax[, c("OTU", "Genus")], by.x = "OTUList", by.y = "OTU", all.x = TRUE, all.y = F)
genus_shared_IDs <- merge(x = tax[, c("OTU", "Genus")] , y = genus_shared, by.x = "OTU", by.y = "OTUList", all.x = F, all.y = T)

#creating a genus level shared file with the added genus names- note that this is *not* identical to the normal shared format, as the genus labels are added as an additional column.
write.csv(genus_shared_IDs, file = "genus_shared.csv", row.names = F, quote = FALSE)

##########################################################
############## Genus level visualizations  ###############
##########################################################

######################################################################
### This set of visualizations show the composition of each sample ###
######################################################################

#I have to subset the the data set one more time
data_genus <- data_sub_type %>%
  tax_glom(taxrank = "Genus") %>%
  
  transform_sample_counts(function(x) {(x/sum(x))*100} ) %>% 
  
  # Transform to rel. abundance
  psmelt() %>%                                         
  # Melt to long format
  arrange(Genus)                                      
  # Sort data frame alphabetically by phylum


###We needed to filter the dataset based on the top abundance file made above###
data_genus <- data_genus[data_genus$Genus %in% top_20_genus,]

data_genus$Sample <- as.character(data_genus$Sample)

#now for the plot
ggplot(data_genus, aes(x = Sample, y = Abundance, fill = Genus)) + 
  #this specifies that you would like to use a bar graph that has black outlines
  geom_bar(stat = "identity", colour = "black") +
  
  facet_wrap(INPUT~., scales = 'free', nrow = 1) +
  
  # Remove x axis title and add a text angle to the x ticks.
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90)) + 
  
  #scale_fill_manual(values=newpalette) +
  
  scale_fill_viridis_d(option = "B") + #this is an alternative to the above scale fill
  #viridis is really nice for generating colorblind friendly figures with nice separation of color
  #to add this option, remove the # before "scale_fill_viridis" and a # before "scale_fill_manual"

  guides(fill = guide_legend(keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (%) \n")


#######################################################################
#### This set of visualizations show the data merged by a variable ####
############## and subset to the most abundant genus ##################
#######################################################################


#I have to subset the the data set one more time
data_genus_subset <- data_sub_type %>%
  tax_glom(taxrank = "Genus") %>%
  
  transform_sample_counts(function(x) {(x/sum(x))} ) %>% 
  
  merge_samples("INPUT") %>%
  
  transform_sample_counts(function(x) {(x/sum(x))*100} ) %>%
  # Transform to rel. abundance
  psmelt() %>%                                         
  # Melt to long format
  arrange(Genus)                                      
# Sort data frame alphabetically by phylum

#once again, we will need to remove all but the 10 most abundant genera (this can be modified)
data_genus_subset <- data_genus_subset[data_genus_subset$Genus %in% top_10_genus,]

ggplot(data_genus_subset, aes(x = Sample, y = Abundance, fill = Genus)) + 
  #this specifies that you would like to use a bar graph that has black outlines
  geom_bar(stat = "identity", colour = "black") +
  
  facet_wrap(INPUT~., scales = 'free', nrow = 1) +
  
  #scale_fill_manual(values=newpalette) +
  
  scale_fill_viridis_d(option = "B") + #this is an alternative to the above scale fill
  #viridis is really nice for generating colorblind friendly figures with nice separation of color
  #to add this option, remove the # before "scale_fill_viridis" and a # before "scale_fill_manual"
  
  guides(fill = guide_legend(keywidth = 1, keyheight = 1)) +
  
  ylab("Relative Abundance (%) \n") +
  
  xlab("INPUT")#
  

##########################################################
########### making family abundance tables  ##############
##########################################################

#I have to subset the the data set one more time
data_family_only <- data_sub_type %>%
  tax_glom(taxrank = "Family")

#I also want to see the number of of sequences within each phylum and their % abundance
#I also want to see the number of of sequences within each phylum and their % abundance
#making a data frame with the genus sums (this has breaks in genus defined by OTU)
pruned_family_sum_df <- data.frame(sum = taxa_sums(data_family_only))
#Making the row names (OTUs) their own column
pruned_family_sum_df <- setDT(pruned_family_sum_df, keep.rownames = TRUE)[]
#naming the new row "OTU"
names(pruned_family_sum_df)[1] <- "OTU"
#merging the Genus column of the modified taxonomy file with the new data frame
merged_family_list <- merge(x = pruned_family_sum_df, y = tax[ , c("OTU", "Family")], by = "OTU", all.x=TRUE)
#calculating relative abundance (currently of the subset data)
merged_family_list$relative_abundance <- (merged_family_list$sum/sum(merged_family_list$sum)*100)
#ordering the genuses based on relative abundance
merged_family_list <- merged_family_list[order(-merged_family_list$relative_abundance),]

write.csv(merged_family_list, file = "relative_abundance_family.csv", row.names = TRUE, quote = FALSE)

##########################################################
########## making kingdom abundance tables  ##############
##########################################################

#this is a similar scenario to both genus and phylum calculations.
#this will tell you how much of the data set is bacterial vs. archeal
#purely descriptive. Not super useful.

#I have to subset the the data set one more time
data_kingdom_only <- data_sub_type %>%
  tax_glom(taxrank = "Kingdom")

pruned_kingdom_sum_df <- data.frame(sum = taxa_sums(data_kingdom_only))
pruned_kingdom_sum_df <- setDT(pruned_kingdom_sum_df, keep.rownames = TRUE)[]
names(pruned_kingdom_sum_df)[1] <- "OTU"
merged_kingdom_list <- merge(x = pruned_kingdom_sum_df, y = tax[ , c("OTU", "Kingdom")], by = "OTU", all.x=TRUE)
merged_kingdom_list$relative_abundance <- (merged_kingdom_list$sum/sum(merged_kingdom_list$sum)*100)
merged_kingdom_list <- merged_kingdom_list[order(-merged_kingdom_list$relative_abundance),]

write.csv(merged_kingdom_list, file = "relative_abundance_kingdom.csv", row.names = TRUE, quote = FALSE)

#################################################
###### OTU Abundance Comparison Graphics ########
#################################################

#UPDATE TO SAS OUTPUT
#THIS SECTION has been moved to its own script. (see OTU comparison)

#################################################
############### Permanova #######################
#################################################
#set the random seed (see above)
set.seed(941996)
#randomly subsample the data to a specific read depth
data_scaled = rarefy_even_depth(data_sub_type, sample.size = LOWEST_SAMPLE_SIZE)
#calculate the bray-curtis community distances for each sample (you can use other distance calculations here)
data_bray <- phyloseq::distance(data_scaled, method = "bray")
#making the phyloseq object a dataframe
data_sub_type_df <- data.frame(sample_data(data_sub_type))
#running adonis to compare microbial communities based on a variable.
#this will output results similar to a standard anova test
  #a significant (p<0.05) difference will cause pairwise tests to be run, determining which group means are different
  #update: this needs to be the adonis2() command
adonis2(data_bray ~ INPUT + INPUT + INPUT*INPUT, data = data_sub_type_df)
#betadisperser should be run with the permanova every time!
#this test will show if the differences in means are only due to differences in variation
  #permanova shows both differences in means and differences in variation
data_beta <- betadisper(data_bray, data_sub_type_df$INPUT)
#this will display the results of the betadisperser test.
#please read about this command online for interpreting results.
permutest(data_beta)
#for pairwise comparisons dispersion between levels, you can use the following co
beta.disperser <- permutest(data_beta, pairwise = TRUE)
#please read about this command online for interpreting results.
permutest(data_beta)
#making all of the pairwise comparisons between a given levels.
pairwise.adonis(data_bray, sample_data(data_sub_type)$Combo)

##########################################################
###### making a heatmap of OTU relabund per sample  ######
##########################################################

#(re)create the otu_table_new from the subsampled phyloseq object to a dataframe
otu_table_new <- as.data.frame(otu_table(data_sub_type))
sample_sum_df_sub <- as.data.frame(colSums(otu_table_new))
sample_sum_df_sub <- setDT(sample_sum_df_sub, keep.rownames = TRUE)[]
colnames(sample_sum_df_sub) <- c("SampleID","total_seqs")

#(re)create the subsampled metadata file
map_df <- read.csv("design/metadata.csv")
map_sub_type_df <- map_df[map_df$INPUT != "LEVEL",]
map_sub_type_df <- map_df[map_df$INPUT != "LEVEL",]

#transpose the otu_table to resemble the original shared file
t_otu_table_new <- as.data.frame(t(otu_table_new))
t_otu_table_new <- setDT(t_otu_table_new, keep.rownames = TRUE)[]
names(t_otu_table_new)[1]<-paste("SampleID")
#truncate to the top 20 OTUs
t_otu_table_trunc <- t_otu_table_new[,1:21]

#merging the truncated otu_table with the subsampled, truncated otu_table and the new sample_sum_df object
t_otu_table_merged <- merge(x = t_otu_table_trunc, y = sample_sum_df_sub[ , c("SampleID", "total_seqs")], by=c("SampleID"))  

#making a dataframe with the relative abundance of each OTU per sample
t_otu_table_merged_relabund <- (t_otu_table_merged[,2:21]/t_otu_table_merged$total_seqs)*100
#adding the sampleIDs back in
t_otu_table_merged_relabund$SampleID <- t_otu_table_merged$SampleID
#merging the metadata information that you specifically want to facet by (INPUT)
t_otu_table_merged_relabund <- merge(x = t_otu_table_merged_relabund, y = map_sub_type_df[ , c("SampleID", "INPUT")], by=c("SampleID"))  
#melting the data based on both SampleID and metadata variable
t_otu_table_merged_melt <- reshape2::melt(data = t_otu_table_merged_relabund, id.vars = c("SampleID", "INPUT"))
#setting the value column to numeric
t_otu_table_merged_melt$value <- as.numeric(t_otu_table_merged_melt$value)
#renaming the third and fourth column
names(t_otu_table_merged_melt)[3]<-paste("OTU")
names(t_otu_table_merged_melt)[4]<-paste("RelativeAbundance")

#(re)create the taxonomy_table from the subsampled phyloseq object to a dataframe 
taxonomy_table <- as.data.frame(tax_table(data_sub_type))
taxonomy_table <- setDT(taxonomy_table, keep.rownames = TRUE)[]
names(taxonomy_table)[1]<-paste("OTU")
#truncating the taxonomy table to the top 20 OTUs
tax_trunc <- taxonomy_table[1:20,]

#adding the OTU and genus classifications for the otu_table + metadata (merged) object
t_otu_table_merged_melt_tax <- merge(x = t_otu_table_merged_melt, y = tax_trunc[ , c("OTU", "Genus")], by=c("OTU"))
t_otu_table_merged_melt_tax$OTU <- paste(t_otu_table_merged_melt_tax$OTU, t_otu_table_merged_melt_tax$Genus, sep=" ")


#ggplot code for the heatmap
otu_heatmap <-  ggplot(t_otu_table_merged_melt_tax, aes(x = SampleID, y = reorder(OTU, desc(OTU)), fill=RelativeAbundance)) + 
  #defining internal ggplot parameters
  geom_tile(colour = "black") +
  #calling up heatmap (tile) that has defined black borders
  theme(axis.title.x = element_blank(),
        strip.text.x = element_text(margin = margin(b = 4, t = 4), size = 12, vjust = 2),
        panel.spacing = unit(0, "lines"),
        axis.text.x = element_text(angle = 75, hjust = 1)) +
  #adjusting text spacing and positioning
  #depending on how many samples you have, you may want to remove xaxis text entirely
  scale_fill_viridis(option = "D") +
  #using viridis as the color palette
  ylab("OTU") +
  #yaxis = OTUs / xaxis = samples
  facet_grid(~ INPUT, scales = "free_x", space = "free_x", switch = "x")
  #faceting by variable
otu_heatmap

##########################################################
################# OTU-Level Statistics ###################
##### prepping SAS import files from phyloseq object #####
##########################################################

#exporting the updated otu_table to an R dataframe
otu_table_new <- as.data.frame(otu_table(data_sub_type))

#calculating the new relative abundance of each OTU after subsetting / removal of sample types.
#this line calculates the row totals for the subset OTU table (removal fo sample types)
otu_counts <- as.data.frame(rowSums(otu_table_new))
#setting it to a data table and adding in a column for OTU names
otu_counts <- setDT(otu_counts, keep.rownames = TRUE)[]
#renaming the columns
colnames(otu_counts) <- c("OTU","total_seqs")
#calculating the (new) relative abundance of each OTU after sample removal and removal of low-end reads
otu_counts$relabund <- (otu_counts$total_seqs / sum(otu_counts$total_seqs))* 100

#creating a new metadata object with the same subsampling. (different from the object above to avoid conflict)
#map_df <- read.csv("design/metadata.csv")
#map_sub_type_df <- map_df[map_df$INPUT != "LEVEL",]
#map_sub_type_df <- map_sub_type_df[map_sub_type_df$INPUT != "LEVEL",]

#recalculating the sample_sum_df object after removing sample types.
sample_sum_df_sub <- as.data.frame(colSums(otu_table_new))
sample_sum_df_sub <- setDT(sample_sum_df_sub, keep.rownames = TRUE)[]
colnames(sample_sum_df_sub) <- c("SampleID","total_seqs")

#calculating the log value of the total_seqs (used as an input for SAS)
sample_sum_df_sub$Count <- (log(sample_sum_df_sub$total_seqs))

#appending the sample_sum_df_sub$Count (log value for total sequences) column to the metadata file
map_sub_type_df <- merge(x = map_sub_type_df, y = sample_sum_df_sub[ , c("SampleID", "Count")], by=c("SampleID"))

#writing the new metadata file for sas
write.csv(map_sub_type_df, file = "metadata_sas.csv", row.names = FALSE, quote = FALSE)

#making the updated shared file for SAS
otu_table_short <- as.data.frame(t(otu_table_new))
otu_table_short <- otu_table_short[,1:100]
otu_table_short <- log(otu_table_short)
otu_table_short <- setDT(otu_table_short, keep.rownames = TRUE)[]
names(otu_table_short)[1]<-paste("Group")
#Replacing all "inf" with 0. THIS MAY NOT BE THE IDEAL MOVE.
otu_table_short <- data.frame(lapply(otu_table_short, function(x) {gsub("-Inf\\s*", "0", x)}))
otu_table_sas <- otu_table_short

#writing the new shared file for sas
write.csv(otu_table_sas, file = "shared_sas.csv", row.names = FALSE, quote = FALSE)

#making the updated taxonomy file for SAS
taxonomy_table <- as.data.frame(tax_table(data_sub_type))
taxonomy_table <- setDT(taxonomy_table, keep.rownames = TRUE)[]
names(taxonomy_table)[1]<-paste("OTU")

tax_sub_100 <- taxonomy_table[1:100,]

write.csv(tax_sub_100, file = "taxonomy_100_sas.csv", row.names = FALSE, quote = FALSE)

taxonomy_table_relabund <- merge(x = taxonomy_table, y = otu_counts[ , c("OTU", "relabund")], by=c("OTU"))

write.csv(taxonomy_table_relabund, file = "subsampled_taxonomy.csv", row.names = FALSE, quote = FALSE)

##########################################################
################  Genus-Level Statistics #################
##### prepping SAS import files from phyloseq object #####
##########################################################

#we should have already converted genus_shared from a matrix to a data frame.
#But just to be sure
genus_shared_df <- as.data.frame(genus_shared)

#need to remove non-numeric column for next command
genus_shared_df <- within(genus_shared_df, rm(OTUList))

#calculating the new relative abundance of each genus 
#calculating the row sum for each genus
genus_counts <- rowSums(genus_shared_df)

#making genus_counts into a dataframe 
genus_counts <- as.data.frame(genus_counts)

#setting it to a data table and adding in a column for OTU names
genus_counts <- setDT(genus_counts, keep.rownames = TRUE)[]

#renaming the columns
colnames(genus_counts) <- c("Genus","total_seqs")

#calculating the (new) relative abundance of each OTU 
genus_counts$relabund <- (genus_counts$total_seqs / sum(genus_counts$total_seqs))* 100

#recalculating the sample_sum_df object after removing sample types.
sample_genus_sum_df <- as.data.frame(colSums(genus_shared_df))
sample_genus_sum_df <- setDT(sample_genus_sum_df, keep.rownames = TRUE)[]
colnames(sample_genus_sum_df) <- c("SampleID","total_seqs")

#calculating the log value of the total_seqs (used as an input for SAS)
sample_genus_sum_df$Count <- (log(sample_genus_sum_df$total_seqs))

#appending the sample_genus_sum_df$Count (log value for total sequences) column to the metadata file
map_genus_df <- merge(x = map_sub_type_df, y = sample_genus_sum_df[ , c("SampleID")], by=c("SampleID"))

#writing the new metadata file for sas
write.csv(map_genus_df, file = "metadata_genus_sas.csv", row.names = FALSE, quote = FALSE)

#making the updated shared file for SAS
genus_table <- as.data.frame(genus_shared)

#need to remove non-numeric column for next command
genus_table <- within(genus_table, rm(OTUList))

#flipping columns to rows 
genus_table_short <- as.data.frame(t(genus_table))

#check to see where the number of sequences for each genera drops off
#you don't want to be comparing against 0
#Because of this, like with the OTUs, the top 100 genera should be sufficient
#Your total number of genera may also be below 100. 
#Include all genera with a satisfactory number of reads
#Adjust genus_table_short[,1:100] according to the number of genera (columns) you want to keep 
genus_table_short <- genus_table_short[,1:100]
genus_table_short <- log(genus_table_short)
genus_table_short <- setDT(genus_table_short, keep.rownames = TRUE)[]
names(genus_table_short)[1]<-paste("Group")
#Replacing all "inf" with 0. THIS MAY NOT BE THE IDEAL MOVE.
genus_table_short <- data.frame(lapply(genus_table_short, function(x) {gsub("-Inf\\s*", "0", x)}))
genus_table_sas <- genus_table_short

#writing the new shared file for sas
write.csv(genus_table_sas, file = "shared_genus_sas.csv", row.names = FALSE, quote = FALSE)

#making the updated taxonomy file for SAS
genus_taxonomy_df <- as.data.frame(genus_shared)
taxonomy_table <- as.data.frame(tax_table(data_sub_type))
taxonomy_table <- setDT(taxonomy_table, keep.rownames = TRUE)[]
names(taxonomy_table)[1]<-paste("OTU")

#appending the sample_genus_sum_df$Count (log value for total sequences) column to the metadata file
genus_tax_short <- genus_taxonomy_df$OTUList
genus_tax_short <- as.data.frame(genus_tax_short)
names(genus_tax_short)[1]<-paste("OTU")
genus_taxonomy_sum_df <- merge(genus_tax_short, taxonomy_table, by.x = "OTU", by.y = "OTU")

#Adjust genus_taxonomy_sum_df[1:100,] according to the number of genera you want to keep 
genus_tax_sub_100 <- genus_taxonomy_sum_df[1:100,]

#making the excel file to be brought into R
write.csv(genus_tax_sub_100, file = "taxonomy_100_genus_sas.csv", row.names = FALSE, quote = FALSE)

##########################################################
## OTUID graphics (modified genus, experimental module) ##
##########################################################
################## NOTE: THIS SHOULD ONLY BE USED IF YOUR GENUS LEVEL GRAPHICS ARE SPECIFICALLY GIVING YOU PROBLEMS- in rare cases this can cause it's own problems
##(IE if multiple bad genera tags are removed from the same family within your top X OTUs you will get the multiple bands problem)
## first lets make an spare copy of your phyloseq object from the earlier genus
data_genus_only_OTUID <- data_genus_only
### now lets remove any problematic genus level tags (for me it was "uncultured"). add additional gsubs to get rid of your problem genus tags as NAs
data_genus_only_OTUID@tax_table =gsub("uncultured", as.character(NA), data_genus_only_OTUID@tax_table)
data_genus_only_OTUID@tax_table =gsub("unknown_unclassified", as.character(NA), data_genus_only_OTUID@tax_table)
### now lets make out OTUID list- this uses an iterative loop to generate an OTUID that is the lowest taxonomic ID for that OTU which is not an NA
OTUlist <- list()
for (i in 1:nrow(data_genus_only_OTUID@otu_table)){
  taxstrings <- as.character(data_genus_only_OTUID@tax_table[i])
  tax_level = 6
  tax_name <- NA
  while(is.na(tax_name))
  {
    tax_name <- taxstrings[tax_level]
    tax_level <- tax_level -1
  }
  OTUlist[[i]] <- tax_name
}
### this generated OTUlist which is the list of OTUIDs: now lets edit the tax_table in our new phylhoseq to match.
tempTax <- as(tax_table(data_genus_only_OTUID), "matrix")
tempTax <- as.data.frame(tempTax)
tempTax$OTUID <-OTUlist
tempTax2 <- transform(tempTax, OTUID=unlist(OTUID))
tempTax3 <- as.matrix(tempTax2, rownames=TRUE)
tax_table(data_genus_only_OTUID) <- tempTax3
#### now your phyloseq object should have the OTUID column alongside the normal taxonomic columns. check this and make sure it looks right!
### assuming everything looks good, lets move on to a modified version of the genus level code- we need to make our lists of top x to filter with
### now lets procedede with the (slightly modified) copied code from the genus level work.
#generating a dataframe with the sum of each taxa from the merge, and the info from the tax table (I really only need OTUID but I'm not sure how to pull a single column from the tax table)
pruned_genus_OTUID_sum_df <- data.frame(sum = taxa_sums(data_genus_only_OTUID), OTUID=data_genus_only_OTUID@tax_table)
#calculating relative abundance (currently of the subset data)
pruned_genus_OTUID_sum_df$relative_abundance <- (pruned_genus_OTUID_sum_df$sum/sum(pruned_genus_OTUID_sum_df$sum)*100)
#ordering the genuses/OTUIDs based on relative abundance
pruned_genus_OTUID_sum_df <- pruned_genus_OTUID_sum_df[order(-pruned_genus_OTUID_sum_df$relative_abundance),]
pruned_genus_OTUID_sum_df <- data.frame(lapply(pruned_genus_OTUID_sum_df, as.character), stringsAsFactors=FALSE)

#write to csv
write.csv(merged_genus_list, file = "relative_abundance_genus_OTUID.csv", row.names = TRUE, quote = FALSE)

#make a top 10/x list
top_10_OTUID <- as.list(pruned_genus_OTUID_sum_df$OTUID.OTUID[1:10])
top_15_OTUID <- as.list(pruned_genus_OTUID_sum_df$OTUID.OTUID[1:15])
top_20_OTUID <- as.list(pruned_genus_OTUID_sum_df$OTUID.OTUID[1:20])
top_25_OTUID <- as.list(pruned_genus_OTUID_sum_df$OTUID.OTUID[1:25])
top_30_OTUID <- as.list(pruned_genus_OTUID_sum_df$OTUID.OTUID[1:30])
top_50_OTUID <- as.list(pruned_genus_OTUID_sum_df$OTUID.OTUID[1:50])

#### OK now lets adapt the actual graphics part! First Merge by whatever relevant characteristic (in my case I was interested in Animal or Study)
data_genus_OTUID_subset <- data_genus_only_OTUID %>%
  merge_samples("INPUT") %>%
  transform_sample_counts(function(x) {(x/sum(x))*100} ) %>%
  # Transform to rel. abundance
  psmelt() %>%
  # Melt to long format
  arrange(OTUID)
# Sort data frame alphabetically by phylum
##then get rid of things outside your top x list
data_genus_OTUID_subset <- data_genus_OTUID_subset[data_genus_OTUID_subset$OTUID %in% top_50_OTUID,]
## then graph
ggplot(data_genus_OTUID_subset, aes(x = Sample, y = Abundance, fill = OTUID)) +
  geom_bar(stat = "identity", colour = "black") +
  #scale_fill_manual(values=newpalette) +
  scale_fill_viridis_d(option = "B") +
  guides(fill = guide_legend(keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance % (Genera) \n") +
  xlab("INPUT") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
### top 30
data_genus_OTUID_subset <- data_genus_OTUID_subset[data_genus_OTUID_subset$OTUID %in% top_30_OTUID,]
## then graph
ggplot(data_genus_OTUID_subset, aes(x = Sample, y = Abundance, fill = OTUID)) +
  geom_bar(stat = "identity", colour = "black") +
  #scale_fill_manual(values=newpalette) +
  scale_fill_viridis_d(option = "B") +
  guides(fill = guide_legend(keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance % (Genera) \n") +
  xlab("INPUT") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
### top 25
data_genus_OTUID_subset <- data_genus_OTUID_subset[data_genus_OTUID_subset$OTUID %in% top_25_OTUID,]
## then graph
ggplot(data_genus_OTUID_subset, aes(x = Sample, y = Abundance, fill = OTUID)) +
  geom_bar(stat = "identity", colour = "black") +
  #scale_fill_manual(values=newpalette) +
  scale_fill_viridis_d(option = "B") +
  guides(fill = guide_legend(keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance % (Genera) \n") +
  xlab("INPUT") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
## top 20
data_genus_OTUID_subset <- data_genus_OTUID_subset[data_genus_OTUID_subset$OTUID %in% top_20_OTUID,]
## then graph
ggplot(data_genus_OTUID_subset, aes(x = Sample, y = Abundance, fill = OTUID)) +
  geom_bar(stat = "identity", colour = "black") +
  #scale_fill_manual(values=newpalette) +
  scale_fill_viridis_d(option = "B") +
  guides(fill = guide_legend(keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance % (Genera) \n") +
  xlab("INPUT") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
### top 15
data_genus_OTUID_subset <- data_genus_OTUID_subset[data_genus_OTUID_subset$OTUID %in% top_15_OTUID,]
## then graph
ggplot(data_genus_OTUID_subset, aes(x = Sample, y = Abundance, fill = OTUID)) +
  geom_bar(stat = "identity", colour = "black") +
  #scale_fill_manual(values=newpalette) +
  scale_fill_viridis_d(option = "B") +
  guides(fill = guide_legend(keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance % (Genera) \n") +
  xlab("INPUT") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
## top 10
data_genus_OTUID_subset <- data_genus_OTUID_subset[data_genus_OTUID_subset$OTUID %in% top_10_OTUID,]
## then graph
ggplot(data_genus_OTUID_subset, aes(x = Sample, y = Abundance, fill = OTUID)) +
  geom_bar(stat = "identity", colour = "black") +
  #scale_fill_manual(values=newpalette) +
  scale_fill_viridis_d(option = "B") +
  guides(fill = guide_legend(keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance % (Genera) \n") +
  xlab("INPUT") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

##########################################################
##                      Examples                        ##
##########################################################

#################### Example of alpha diversity linear model w/ random variables and interaction ###################

#Begin by building a model
alphadiv.lme.model <- lme(data=molten_alpha_meas_full %>% 
                            #filter by the alpha diversity metric you want to test (WILL HAVE TO RUN THIS MULTIPLE TIMES)                         
                            filter(variable=="Shannon"), 
                          #including the fixed effect that you wish to test for (including interactions) 
                          value ~ bred + sampling + bred*sampling,
                          #including any random or repeated variables
                          random=~1|Sheep_ID) 
#Testing the separation in means for each variable specified above, correcting for any random variables
Anova(alphadiv.lme.model, type="III")
#getting pariwise estimates of the LSMeans
lsmeans(alphadiv.lme.model, list(pairwise ~ bred|sampling))
#to get Least Squares Means of fixed effects (should only use if no interactions are significant)
#lsmeans(alphadiv.lme.model, list(~bred, ~sampling))
#writing the results to an object then extracting them as a data frame
alpha_div_results <- lsmeans(alphadiv.lme.model, list(pairwise ~ bred|sampling))
alpha_div_input <- as.data.frame(summary(alpha_div_results$`lsmeans of bred | sampling`))
#adding a column to specify the levels
alpha_div_input$combo <- with(alpha_div_input, paste0(bred, "_", sampling))

#Plotting the output
ggplot(alpha_div_input, aes(x = combo, y = lsmean, fill = combo)) + 
  #this specifies that you would like to use a bar graph that has black outlines
  geom_bar(stat = "identity", colour = "black") +
  #this option includes the sampleIDs along the x-Axis
  theme(axis.text.x = element_text(vjust=0.5, angle = 90)) +
  scale_fill_manual(values= c("goldenrod", "aquamarine3", "firebrick", "darkorchid3")) +
  #using the upper and lower confidence intervals produced by the lsmeans command as error bars
  geom_errorbar(aes(ymin=alpha_div_input$lower.CL, 
                    ymax=alpha_div_input$upper.CL), width=.2) +
  #add a title that matches the variable of interest
  xlab("Shannon Diversity")

#################### Additional Reading      ###################

https://fromthebottomoftheheap.net/slides/advanced-vegan-webinar-2020/advanced-vegan#51
