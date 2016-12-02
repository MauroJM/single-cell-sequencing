source("~/AvO_lab/R/scripts/SCseq_functions2.R")
# install.packages("oce")

require(RColorBrewer)
require(oce)

corner<-c(357:360,381:384) # specify the location of your empty wells (following primer order)

#if needed, merge CS1-type libraries into one 384-column long CS2-type file 
#this assumes the order of A1=lib1, A2=lib2, B1=lib3, B2=lib4 when it reorders the columns
# reorder.cs1(c("EP10-12-1_AHVK72BGXX_S1","EP10-12-2_AHVK72BGXX_S2","EP10-12-3_AHVK72BGXX_S3","EP10-12-4_AHVK72BGXX_S4"),"EP10-live")

####diagnostics file####
# file to read in plates from SORT-seq and do diagnostics on them. 
# set the working directory to where you keep your .coutc/.coutb/coutt files

#read data into four lists
files <- read_files(dir = getwd()) # path to all files with .cout* extention

#variables for the script
names<-list() # this will store the "workable" names of the plates
split_files<-list() # full names of the plates without .cout* extension
tc<-list() # list of .coutt files
rc<-list() # list of .coutc files
bc<-list() # list of .coutb files
  
# read in files
split_files <-sub(".*\\/","",files)
for(i in 1:length(files)){
  names[[i]] <-  sub("\\_.*","",split_files[[i]])
  cat("\n",split_files[[i]],"was renamed to",names[[i]],"\n",sep = " ")
  cat("reading .cout files for plate",i, "out of", length(files),"\n",sep = " ")
  tc[[i]] <- read.csv(paste(files[i],".coutt.csv", sep=""), header = TRUE, sep = "\t",row.names =1)
  rc[[i]] <- read.csv(paste(files[i],".coutc.csv", sep=""), header = TRUE, sep = "\t",row.names =1) 
  bc[[i]] <- read.csv(paste(files[i],".coutb.csv", sep=""), header = TRUE, sep = "\t",row.names =1)
  cat("library",names[[i]],"contains a total of",nrow(tc[[i]]),"genes")
}

#lable cells
for(i in 1:length(tc)){
  l<-c(1:384)
   for(j in 1:384){
   l[j]<-paste(names[[i]],j,sep="_")
    }
  colnames(tc[[i]])<-l
}

####diagnostic plots
path<-"/Users/mauro/Desktop/testfiles/" # change working directory if needed
dir.create(path, showWarnings = F) 
setwd(path)
for(i in 1:length(tc)){
  pdf(paste(names[[i]],"_plate_diagnostics",".pdf",sep=""))
  par(mfrow = c(3,3))
  totalreads(tc[[i]],plotmethod = "hist")
  cellgenes(tc[[i]],plotmethod= "cumulative")
  overseq2(rc[[i]],bc[[i]])
  plate.plots(tc[[i]])
  topgenes(tc[[i]])
  leakygenes(tc[[i]])
  dev.off()
} #make pdf with diagnostic plots

#merge all dataframes of object into one and write a .csv
cdata_all<-tc[[1]]
for (i in 2:length(tc)){
  cdata_all <- intersectmatrix(cdata_all,tc[[i]]) #  JC's function intersectmatrix
}
cdata_all<- cdata_all[order(rownames(cdata_all)), ] #make row names alphabetical

write.table(cdata_all,"/Users/mauro/AvO_lab/R/data files/GK2E_GK2F.csv", sep="\t")


