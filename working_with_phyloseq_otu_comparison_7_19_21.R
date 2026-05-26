##########################################################
############ Working with Phyloseq Master ################
######### Reading in OTU comparisons from SAS ############
################### LRK 7/19/2021 ########################
#                      R 4.0.5                           #
##########################################################

#This script is for visualizing OTU comparison results from SAS.

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

#install.packages("ggplot2")
#install.packages("data.table")
#install.packages("dplyr")
#install.packages("reshape2")
#install.packages("dplyr")
#install.packages("splitstackshape")
#install.packages("viridis") #for nice color palettes
#install.packages("remotes")
#install.packages("devtools")
#library(devtools)
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
#install.packages("DescTools")
#install.packages("gtools")

##########################################################
################## Loading packages  #####################
##########################################################

#LOADING PACKAGES
#now that you have installed them, you will need to load them to this session. 
#unlike installing, packages will need to be loaded every session.

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
library(DescTools)
library(gtools)

##########################################################
############### reading in SAS outputs ###################
###### and subsetting OTUs based on significance #########
##########################################################

#the following code requires that you have run Amy's SAS code 
#and exported all of the outputs to separate .csv files

#importing
pvalues <- read.csv("pvalues.csv")
qvalues <- read.csv("qvalues.csv")
confidence_intervals <- read.csv("CI.csv")
estimates <- read.csv("estimates.csv")
LSM <- read.csv("LSM.csv")
pdiff <- read.csv("pdiff.csv")
tax <- read.csv("taxonomy_100_sas.csv")
#I will now add a "Test" column in p values to match the Test column in qvalues
#this is necessary, because some tests to not work/converge(?) when creating the q values, and thus some tests are not included.
#you may notice in the qvalues$Test some numbers are skipped
  #Update: I figure out this is due to the fact that SAS puts a "<" sign to denote very significant. This changes the vaules to a character not a numeric and does not import them properly
    #manually remove the "<" from the p / q values file
pvalues$Test <- 1:nrow(pvalues)
#now I will merge the FDR correction column from qvalues to pvalues by the Test column
pvalues_merge <- merge(x = pvalues, y = qvalues[ , c("FalseDiscoveryRate", "Test")], by=c("Test"))
#removing any rows without information (NAs)
pvalues_merge <- pvalues_merge[complete.cases(pvalues_merge), ]
#Now I will make a subset of significant tests
pvalues_merge_sub <- pvalues_merge[pvalues_merge$FalseDiscoveryRate <= 0.05,]
#merging the pvalues_merge_sub with the estimates
LSM_merged_sub <- merge(x = LSM, y = pvalues_merge_sub[ , c("OTU", "Effect")], by=c("Effect", "OTU"))

#reverse log transforming the estimates
LSM_merged_sub$reverselogestimates <- exp(LSM_merged_sub$Estimate)*100
#adding an OTU number column to the tax file
tax$OTU <- 1:nrow(tax)
#converting to a data.table
tax <- as.data.table(tax)
################################ removing entire variables (optional) #####################################

#You can make graphics for one variable, then change this to another later
#LSM_merged_sub <- LSM_merged_sub[LSM_merged_sub$Effect != "INPUT",]

################################ removal of OTUs with interaction ()
#removal of the OTUs from INPUT that are significant for interaction as well.
#This is necessary for the code to work, but it is not an outright removal of interaction
#to look specifically at the subset interaction OTUs see section below 

#identifying the lines in the data frame that contain the "*" denoting an interaction effect
LSM_merged_sub$interaction <- grepl(pattern = "\\*", LSM_merged_sub$Effect)
  #subsetting a dataframe containing ONLY the interaction
  LSM_merged_sub_interaction <- LSM_merged_sub[LSM_merged_sub$interaction == "TRUE",]
  #subsetting the dataframe containing EVERYTHING BUT the interaction
  LSM_merged_sub_noint <- LSM_merged_sub[LSM_merged_sub$interaction != "TRUE",]
  #identifying repeated OTU patterns between the two separate dataframes and removing any OTUs that match
  LSM_merged_sub_cleaned <- LSM_merged_sub_noint[!(LSM_merged_sub_noint$OTU %in% LSM_merged_sub_interaction$OTU),]

#calculating confidence intervals
  #subsetting the dataframe to only the rows that contain "mean" in the parameter column
  confidence_intervals_sub <- confidence_intervals[confidence_intervals$Parameter == "Mean",]
  #taking the reverse natural log of the upper and lower confidence interval columns
  confidence_intervals_sub[,7:8] <- exp(confidence_intervals_sub[,5:6])*100
  #calculating the SEM based on the transformed CI estimates
  confidence_intervals_sub$SEM <- ((confidence_intervals_sub$UpperCL.1-confidence_intervals_sub$LowerCL.1)/3.92)
  
  #taking the SEM to the log2 to match scales with the log2FCs later on (this is not right)
  #confidence_intervals_sub$CI <- log(confidence_intervals_sub$SEM, 2)
  
  #It's the absolute value of the Log 2 basis difference of the CI when normalized to the distribution
  confidence_intervals_sub$CI <- log((confidence_intervals_sub$UpperCL-confidence_intervals_sub$LowerCL), 2)
  
  #then adjusting it to an absolute value
  confidence_intervals_sub$CI <- abs(confidence_intervals_sub$CI)

##########################################################
# for comparing only the OTUs affected by a fixed effect #
#####           ~no interaction effect~              #####
##########################################################
  
#for OTUs with a significant fixed effect and no interaction, it isn't correct to display the estimates for
#each level.
    #We will instead switch to using the estimates of the fixed effect only
  
#begin by retaining only the fixed effect you want to test
LSM_merged_sub_cleaned_sub <- LSM_merged_sub_cleaned[LSM_merged_sub_cleaned$Effect == "INPUT",]
#convert it to a data table so you can dcast the data
LSM_merged_sub_cleaned_sub <- as.data.table(LSM_merged_sub_cleaned_sub)
#dcasting the data (going from long to short data)
log2_comparisons_wide <- dcast(LSM_merged_sub_cleaned_sub, OTU ~ INPUT, value.var="reverselogestimates")
#pasting OTU before every number in the OTU category
log2_comparisons_wide$OTUnames <- with(log2_comparisons_wide, paste0("OTU", OTU))
#calculating the log2 fold-change between levels
  #in this case, if the resulting output is positive the estimate was higher in LEVEL1 (see how to subtract log values)
log2_comparisons_wide$log2FC_INPUT <- log((log2_comparisons_wide$LEVEL1 / log2_comparisons_wide$LEVEL2), 2)
#bringing in the genus classification for the OTUs
log2_comparisons_wide <- merge(x = log2_comparisons_wide, y = tax[,c("OTU", "Genus")], by=c("OTU"))
#combining the genus classification with the OTU column for added info in the resulting graphic
log2_comparisons_wide$OTUnames <- with(log2_comparisons_wide, paste0(OTUnames, " ", Genus))
#combining the confidence interval with the log2 fold-change values
log2_comparisons_wide <- merge(x = log2_comparisons_wide, y = confidence_intervals_sub[,c("OTU", "CI")], by=c("OTU"))
#adding a "color by" column based positive or negative values
log2_comparisons_wide$comparison <- ifelse(log2_comparisons_wide$log2FC_INPUT >=0, "LEVEL1", "LEVEL2")
#Converting all log2 Fold-change values to absolute values (can differentiate using the column above)
log2_comparisons_wide$abslog2FC_INPUT <- abs(log2_comparisons_wide$log2FC_INPUT)

#Plotting the comparisons
ggplot(data = log2_comparisons_wide, aes(y=reorder(OTUnames, -OTU), x=abslog2FC_INPUT, fill= comparison)) +  
  geom_bar(stat = "identity",  position = "dodge", colour="black") + 
  
  #NOTE: set these scales to whatever is necessary for your data
  scale_x_continuous(limits = c(0, 8), breaks = c(0:8)) +
  scale_fill_manual(values=c("gray", "white")) +
  guides(fill = guide_legend(reverse = TRUE)) +
  labs(fill = "Comparisons") +
  theme(axis.title.y = element_blank())+
  
  #NOTE: when using error bars for absolute log fold-changes, only the positive error bar should be shown
  geom_errorbar(aes(xmin=abslog2FC_INPUT, 
                    xmax=log2_comparisons_wide$abslog2FC_INPUT+log2_comparisons_wide$CI), width=.2)

#This section can be copied and pasted for the other variables of interest

##########################################################
# for comparing only the OTUs affected by an interaction #
#####           ~yes interaction effect~             #####
##########################################################

#After discussing with Amy, it is best to show the entire data relative to the interaction within the same figure
  #this is important to show readers the potential effect of both important variables.

#marking any OTU that has a "*" in the effects line. This signifies there is significant interaction between variables
pdiff$interaction <- grepl(pattern = "\\*", pdiff$Effect)
##subsetting a dataframe containing ONLY the interaction
pdiff_interaction <- pdiff[pdiff$interaction == "TRUE",]
#listing the unique entries for OTUs that have a significant interaction effect
OTUs_with_int <- unique(LSM_merged_sub_interaction[c("OTU")])
#shorting the pdiffs object to only those OTUs with a significant interaction effect
pdiff_interaction_sub <- merge(x = pdiff_interaction, y = OTUs_with_int, by=c("OTU"))
#Now I will make a subset of significant tests
pdiff_interaction_cleaned <- pdiff_interaction_sub[pdiff_interaction_sub$Probt <= 0.05,]
#I am making a combination column so I can pull level comparisons out later
pdiff_interaction_cleaned$combo <- with(pdiff_interaction_cleaned, paste0(INPUT, "_", INPUT, "_", X_INPUT, "_", X_INPUT))
#making a list of unique level combinations
unique_combinations <- unique(pdiff_interaction_cleaned[c("combo")])
#making the data.frame a data.table
pdiff_interaction_cleaned <- as.data.table(pdiff_interaction_cleaned)

#for each comparison you will need to pull rows of that comparison then make a list of unique OTUs

###############################################################################################################################
################################################################# Interaction #################################################
####################################### 1. comparing two levels of one variable with a single level ##########################################
#######################################              for the interaction variable                   ##########################
###############################################################################################################################

#example: I want to compare pregnancy status (2 levels, yes and no) at a specific time point (3 levels, 1,2 and 3)
    #in this case, you can remove levels 2 and 3 from time point

#selecting the first combination you would like to compare. Check the combo list for precise input
pdiff_INPUT_1_INPUT_2 <- pdiff_interaction_cleaned[pdiff_interaction_cleaned$combo == "INPUT_1_INPUT_2"]
#merging with the subset pdiffs file
LSM_merged_sub_int_cleaned <- merge(x = LSM_merged_sub_interaction, y = pdiff_INPUT_1_INPUT_2[,c("OTU")], by=c("OTU"))
#converting to a data.table
LSM_merged_sub_int_cleaned <- as.data.table(LSM_merged_sub_int_cleaned)
#more cleaning. retaining all of the instances of one specific variable of interest. 
  #example: if I want to compare 2 levels of pregnancy (yes and no) for a single time point (1 of 3) I would put:
    # "LSM_merged_sub_int_cleaned$timepoint == "1"" to retain only info for a single time point. See below if it is a larger cross
  #please check the "LSM_merged_sub_int_cleaned" object for info. After merging, additional unwanted levels may be present
LSM_merged_sub_int_cleaned_short <- LSM_merged_sub_int_cleaned[LSM_merged_sub_int_cleaned$INPUT == "1"]
#dcasting the data by the comparison of interest
log2_comparisons_int_subset_wide <- dcast(LSM_merged_sub_int_cleaned_short, OTU ~ INPUT, value.var="reverselogestimates")
#Appending OTU before the number
log2_comparisons_int_subset_wide$OTUnames <- with(log2_comparisons_int_subset_wide, paste0("OTU", OTU))
#calculating the log 2 fold change
  #in this case, you are only comparing the levels within one variable, so you can jus put the levels of that variable for LEVEL
log2_comparisons_int_subset_wide$log2FC <- log((log2_comparisons_int_subset_wide$LEVEL / log2_comparisons_int_subset_wide$LEVEL), 2)
#bringing in the genus classification information
log2_comparisons_int_subset_wide <- merge(x = log2_comparisons_int_subset_wide, y = tax[,c("OTU", "Genus")], by=c("OTU"))
#appending the genus classification information
log2_comparisons_int_subset_wide$OTUnames <- with(log2_comparisons_int_subset_wide, paste0(OTUnames, " ", Genus))
#merging in the confidence interval information
log2_comparisons_int_subset_wide <- merge(x = log2_comparisons_int_subset_wide, y = confidence_intervals_sub[,c("OTU", "CI")], by=c("OTU"))
##adding a "color by" column based positive or negative values
log2_comparisons_int_subset_wide$comparison <- ifelse(log2_comparisons_int_subset_wide$log2FC >=0, "INPUT-1", "INPUT-2")
#calculating and adding a column for the absolute value
log2_comparisons_int_subset_wide$abslog2FC <- abs(log2_comparisons_int_subset_wide$log2FC)

#plotting similar to above
ggplot(data = log2_comparisons_int_subset_wide, aes(y=reorder(OTUnames, -OTU), x=abslog2FC, fill= comparison)) +  
  geom_bar(stat = "identity",  position = "dodge", colour="black") + 
  scale_x_continuous(limits = c(0, 12), breaks = c(0:12)) +
  scale_fill_manual(values=c("aquamarine3", "goldenrod")) +
  guides(fill = guide_legend(reverse = TRUE)) +
  labs(fill = "Comparisons") +
  theme(axis.title.y = element_blank())+
  geom_errorbar(aes(xmin=abslog2FC, 
                    xmax=log2_comparisons_int_subset_wide$abslog2FC+log2_comparisons_int_subset_wide$CI), width=.2)


###############################################################################################################################
################################################################# Interaction #################################################
######################################## 2.  comparing two levels of one variable with 2 levels of ##########################################
########################################                    another variable                       #############################
###############################################################################################################################

#example: I want to compare pregnancy status (2 levels, yes and no) across two time point (3 levels, 1,2 and 3)
  #in this case, you would need the specific combination of levels

#selecting the first combination you would like to compare. Check the combo list for precise input
pdiff_INPUT_1_INPUT_2 <- pdiff_interaction_cleaned[pdiff_interaction_cleaned$combo == "INPUT_1_INPUT_2"]
#merging with the subset pdiffs file
LSM_merged_sub_int_cleaned <- merge(x = LSM_merged_sub_interaction, y = pdiff_INPUT_1_INPUT_2[,c("OTU")], by=c("OTU"))
#converting to a data.table
LSM_merged_sub_int_cleaned <- as.data.table(LSM_merged_sub_int_cleaned)
#more cleaning. retaining all of the instances of one specific variable of interest. 
  #example: if I want to compare 2 levels of pregnancy (yes and no) across two levels of time point (1 and 2) you will need two lines:
    # "LSM_merged_sub_int_cleaned$timepoint == "1"" to retain only info for a single time point. See below if it is a larger cross
  #please check the "LSM_merged_sub_int_cleaned" object for info. After merging, additional unwanted levels may be present
LSM_merged_sub_int_cleaned_short1 <- LSM_merged_sub_int_cleaned [which(INPUT_1 == "LEVEL" & INPUT_2 == "LEVEL"),]
LSM_merged_sub_int_cleaned_short2 <- LSM_merged_sub_int_cleaned [which(INPUT_1 == "LEVEL" & INPUT_2 == "LEVEL"),]
#combining the two objects created above
LSM_merged_sub_int_cleaned_short_combined <- rbind(LSM_merged_sub_int_cleaned_short1, LSM_merged_sub_int_cleaned_short2)
#Similar to above, you can dcast with a single variable here, as you have removed all of one type
log2_comparisons_int_subset_wide <- dcast(LSM_merged_sub_int_cleaned_short_combined, OTU ~ INPUT, value.var="reverselogestimates")
#Appending OTU before the number
log2_comparisons_int_subset_wide$OTUnames <- with(log2_comparisons_int_subset_wide, paste0("OTU", OTU))
#Also similar to above, you can use only a single variable here as you have removed all of one type
log2_comparisons_int_subset_wide$log2FC <- log((log2_comparisons_int_subset_wide$INPUT / log2_comparisons_int_subset_wide$INPUT), 2)
#bringing in the genus classification information
log2_comparisons_int_subset_wide <- merge(x = log2_comparisons_int_subset_wide, y = tax[,c("OTU", "Genus")], by=c("OTU"))
#appending the genus classification information
log2_comparisons_int_subset_wide$OTUnames <- with(log2_comparisons_int_subset_wide, paste0(OTUnames, " ", Genus))
#merging in the confidence interval information
log2_comparisons_int_subset_wide <- merge(x = log2_comparisons_int_subset_wide, y = confidence_intervals_sub[,c("OTU", "CI")], by=c("OTU"))
##adding a "color by" column based positive or negative values
log2_comparisons_int_subset_wide$comparison <- ifelse(log2_comparisons_int_subset_wide$log2FC >=0, "INPUT-1", "INPUT-2")
#calculating and adding a column for the absolute value
log2_comparisons_int_subset_wide$abslog2FC <- abs(log2_comparisons_int_subset_wide$log2FC)

#plotting similar to above
ggplot(data = log2_comparisons_int_subset_wide, aes(y=reorder(OTUnames, -OTU), x=abslog2FC, fill= comparison)) +   
  geom_bar(stat = "identity",  position = "dodge", colour="black") + 
  scale_x_continuous(limits = c(0, 15), breaks = c(0:15)) +
  scale_fill_manual(values=c("purple", "goldenrod")) +
  guides(fill = guide_legend(reverse = TRUE)) +
  labs(fill = "Comparisons") +
  theme(axis.title.y = element_blank())+
  geom_errorbar(aes(xmin=abslog2FC, 
                    xmax=log2_comparisons_int_subset_wide$abslog2FC+log2_comparisons_int_subset_wide$CI), width=.2)

#####################################################  Examples ########################################################


###################################  1. U-S1 v S-S2 ##########################################

pdiff_US1_SS2 <- pdiff_interaction_cleaned[pdiff_interaction_cleaned$combo == "unsuccessful_1_successful_2"]
LSM_merged_sub_int_cleaned <- merge(x = LSM_merged_sub_interaction, y = pdiff_US1_SS2[,c("OTU")], by=c("OTU"))
LSM_merged_sub_int_cleaned <- as.data.table(LSM_merged_sub_int_cleaned)
LSM_merged_sub_int_cleaned_short1 <- LSM_merged_sub_int_cleaned [which(bred == "successful" & sampling == "2"),]
LSM_merged_sub_int_cleaned_short2 <- LSM_merged_sub_int_cleaned [which(bred == "unsuccessful" & sampling == "1"),]
LSM_merged_sub_int_cleaned_short_combined <- rbind(LSM_merged_sub_int_cleaned_short1, LSM_merged_sub_int_cleaned_short2)
log2_comparisons_int_subset_wide <- dcast(LSM_merged_sub_int_cleaned_short_combined, OTU ~ bred, value.var="reverselogestimates")
log2_comparisons_int_subset_wide$OTUnames <- with(log2_comparisons_int_subset_wide, paste0("OTU", OTU))
log2_comparisons_int_subset_wide$log2FC <- log((log2_comparisons_int_subset_wide$successful / log2_comparisons_int_subset_wide$unsuccessful), 2)
log2_comparisons_int_subset_wide <- merge(x = log2_comparisons_int_subset_wide, y = tax[,c("OTU", "Genus")], by=c("OTU"))
log2_comparisons_int_subset_wide$OTUnames <- with(log2_comparisons_int_subset_wide, paste0(OTUnames, " ", Genus))
log2_comparisons_int_subset_wide <- merge(x = log2_comparisons_int_subset_wide, y = confidence_intervals_sub[,c("OTU", "CI")], by=c("OTU"))
log2_comparisons_int_subset_wide$comparison <- ifelse(log2_comparisons_int_subset_wide$log2FC >=0, "P-S2", "NP-S1")
log2_comparisons_int_subset_wide$abslog2FC <- abs(log2_comparisons_int_subset_wide$log2FC)

ggplot(data = log2_comparisons_int_subset_wide, aes(y=reorder(OTUnames, -OTU), x=abslog2FC, fill= comparison)) +   
  geom_bar(stat = "identity",  position = "dodge", colour="black") + 
  scale_x_continuous(limits = c(0, 13), breaks = c(0:13)) +
  scale_fill_manual(values=c("aquamarine3", "firebrick")) +
  guides(fill = guide_legend(reverse = TRUE)) +
  labs(fill = "Comparisons") +
  theme(axis.title.y = element_blank())+ 
  geom_errorbar(aes(xmin=abslog2FC, 
                    xmax=log2_comparisons_int_subset_wide$abslog2FC+log2_comparisons_int_subset_wide$CI), width=.2)


######################################   2. S v S (S1 v S2 only) ##########################################

pdiff_S_S_S1_S2 <- pdiff_interaction_cleaned[pdiff_interaction_cleaned$combo == "successful_1_successful_2"]
LSM_merged_sub_int_cleaned <- merge(x = LSM_merged_sub_interaction, y = pdiff_S_S_S1_S2[,c("OTU")], by=c("OTU"))
LSM_merged_sub_int_cleaned <- as.data.table(LSM_merged_sub_int_cleaned)
LSM_merged_sub_int_cleaned_short <- LSM_merged_sub_int_cleaned[LSM_merged_sub_int_cleaned$bred == "successful"]
LSM_merged_sub_int_cleaned_short$sampling <- with(LSM_merged_sub_int_cleaned_short, paste0("S", sampling))
log2_comparisons_int_subset_wide <- dcast(LSM_merged_sub_int_cleaned_short, OTU ~ sampling, value.var="reverselogestimates")
log2_comparisons_int_subset_wide$OTUnames <- with(log2_comparisons_int_subset_wide, paste0("OTU", OTU))
log2_comparisons_int_subset_wide$log2FC <- log((log2_comparisons_int_subset_wide$S1 / log2_comparisons_int_subset_wide$S2), 2)
log2_comparisons_int_subset_wide <- merge(x = log2_comparisons_int_subset_wide, y = tax[,c("OTU", "Genus")], by=c("OTU"))
log2_comparisons_int_subset_wide$OTUnames <- with(log2_comparisons_int_subset_wide, paste0(OTUnames, " ", Genus))
log2_comparisons_int_subset_wide <- merge(x = log2_comparisons_int_subset_wide, y = confidence_intervals_sub[,c("OTU", "CI")], by=c("OTU"))
log2_comparisons_int_subset_wide$comparison <- ifelse(log2_comparisons_int_subset_wide$log2FC >=0, "P-S1", "P-S2")
log2_comparisons_int_subset_wide$abslog2FC <- abs(log2_comparisons_int_subset_wide$log2FC)

ggplot(data = log2_comparisons_int_subset_wide, aes(y=reorder(OTUnames, -OTU), x=abslog2FC, fill= comparison)) + 
  geom_bar(stat = "identity",  position = "dodge", colour="black") + 
  scale_x_continuous(limits = c(0, 13), breaks = c(0:13)) +
  scale_fill_manual(values=c("goldenrod", "firebrick")) +
  #guides(fill = guide_legend(reverse = TRUE)) +
  labs(fill = "Comparisons") +
  theme(axis.title.y = element_blank())+
  geom_errorbar(aes(xmin=abslog2FC, 
                    xmax=log2_comparisons_int_subset_wide$abslog2FC+log2_comparisons_int_subset_wide$CI), width=.2)


