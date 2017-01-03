# Script for making a pdf containing diagnostic plots of SORT-Seq data. It goes trough the following steps:
# optional: transform CS1 style data to CS2-style layout (from 4 96-cell libraries to 1 384-cell library)
# read in ALL .cout* files in one folder
# lable cells with unique name taken from the input .cout* files and a number for each cell
# optional: rename part of specific plates in the case of two different experiments in one plate
# plot one pdf per plate containing diagnostic plots
# merge specified plates into one file that can be loaded into RaceID

## questions can be addressed to m.muraro@hubrecht.eu

####install/load packages and variables####
# install.packages("oce") # install this package upon first use
source("~/AvO_lab/R/scripts/plate_diagnostics_functions2.R") # specify path to functions file
require(RColorBrewer)
require(oce)

# specify the location of your empty wells (follows primer number order)
# if you don't have empty wells just specify O21-O24 and P21-P24.  
emptywells<-c(357:360,381:384) # this corresponds to O21-O24 and P21-P24

#variables for the script
names<-list() # for the shortened names of the libraries
split_files<-list() # for  full names of the plates without .cout* extension
tc<-list() # list of .coutt files
rc<-list() # list of .coutc files
bc<-list() # list of .coutb files

#### OPTIONAL: Merge CS1-type libraries into one 384-column long CS2-type file####
# do this ONLY if you have CS1-style (4 libraries in one plate) data that you want to run this script on
# this assumes the order of A1=lib1, A2=lib2, B1=lib3, B2=lib4 when it reorders the columns
# setwd("/Users/mauro/AvO_lab/R/GK_GateID/raw_files/GK2E_GK2F/") # set directory containing the files
# CS1files <- read_files(dir = getwd()) # lists all the unique library names in the working directory
# input required is 4 unique library handles in the right order, followed by the name of the new, merged file, like so:
# reorder.cs1(c("EP10-12-1_AHVK72BGXX_S1","EP10-12-2_AHVK72BGXX_S2","EP10-12-3_AHVK72BGXX_S3","EP10-12-4_AHVK72BGXX_S4"),"EP10-CS2-live")

####read files and lable cells ####
#read data into four lists
files <- read_files(dir = getwd()) # path to files with .cout* extention. NB: ALL files in the dir will be processed

# read in files
split_files <-sub(".*\\/","",files) 
for(i in 1:length(files)){
  names[[i]] <-  sub("\\_.*","",split_files[[i]]) # split lib name to keep only name supplied by you to cuppen group 
  cat("\n",split_files[[i]],"was renamed to",names[[i]],"\n",sep = " ") # fyi how the libraries will be named
  cat("reading .cout files for plate",i, "out of", length(files),"\n",sep = " ") # reports progress
  tc[[i]] <- read.csv(paste(files[i],".coutt.csv", sep=""), header = TRUE, sep = "\t",row.names =1)
  rc[[i]] <- read.csv(paste(files[i],".coutc.csv", sep=""), header = TRUE, sep = "\t",row.names =1) 
  bc[[i]] <- read.csv(paste(files[i],".coutb.csv", sep=""), header = TRUE, sep = "\t",row.names =1)
  cat("library",names[[i]],"contains a total of",nrow(tc[[i]]),"genes") # reports number of genes found in each library
}

#lable cells: all cells in library will get a _1 to _384 extension to the library name specified in names object
for(i in 1:length(tc)){ colnames(tc[[i]])<-paste(names[[i]],c(1:384),sep="_") }

# OPTIONAL: rename part of a specified plate if you have different experiments in one plate
# do this only for one of the .coutt files at the time by choosing the correct location in the .coutt list (eg:tc[[5]])
# and manually specify the names of the experiments (eg:"alpha"/"beta") and the cell labels (eg: 1:192), in total should be 384 cells.
names(tc[[5]])<-c(paste("alpha",c(1:192),sep="_"),paste("beta",c(1:192),sep="_")) 

#### diagnostic plots####
path<-"/Users/mauro/Desktop/diagnostic plots/" # specify where you want the diagnostic plots
dir.create(path, showWarnings = F) # directory will be made if it doesn't exist
setwd(path)
for(i in 1:length(tc)){
  pdf(paste(names[[i]],"_plate_diagnostics",".pdf",sep=""))
  par(mfrow = c(3,3)) # specify grid for plots on the pdf 
  totalreads(tc[[i]],plotmethod = "hist") # plots total UMI reads/cell, can choose 4 different plot methods
  cellgenes(tc[[i]],plotmethod= "cumulative") # plot number of detected genes/cell, can choose 4 different plot methods
  overseq2(rc[[i]],bc[[i]]) # plot oversequencing per molecule
  plate.plots(tc[[i]]) # 3 plots: total reads, ERCC reads and division between the two over a plate layout
  topgenes(tc[[i]])  # 2 plots: top expressed and most variable genes
  leakygenes(tc[[i]]) # NB: this function can give errors if you don't have ERCCs or very few succesfully sequenced cells. comment out in case of errors
  # leakygenes plots (1): numer of genes and ERCC reads in the empty corner. Will give warning if a sample has more than plate average genes/5
  # (2): top expressed genes in the empty corner and % of their total reads in empty corner compared to total reads in whole plate.
  # (3): any genes that are in the top 50 genes in empty corner, but not in the top 200 genes in rest of plate (likely artifacts)
  dev.off()
} #make pdf with diagnostic plots

####merge and write multiple dataframes into one .csv####
cdata_all<-tc[[1]] # should be first position you want to start merging from
for (i in 2:length(tc)){
  cdata_all <- intersectmatrix(cdata_all,tc[[i]]) #  JC's function intersectmatrix
} # specify second position to last position you want to merge from-to
cdata_all<- cdata_all[order(rownames(cdata_all)), ] #make row names alphabetical

write.table(cdata_all,"/Users/mauro/AvO_lab/R/data files/DM23.csv", sep="\t") # this file can be loaded into RaceID


