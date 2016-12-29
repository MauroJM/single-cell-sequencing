#Functions needed for reading, filtering and normalizing Cel-SEQ data

#remove or keep only spike ins from given data frame
rmspike<-function(x){
  ERCCs<-grep("ERCC-",row.names(x)) # gives vector with row # of spike ins
  data<-x[-ERCCs,] # make new data frame without the specified rows
  return(data) # output new data frame
}
keepspike<-function(x){
  ERCCs<-grep("ERCC-",row.names(x)) # gives vector with row # of spike ins
  data<-x[ERCCs,] # make new data frame with only the specified rows
  return(data) # output new data frame
}

# chop of chromosome lables (Abel)
chop_chr <- function  (name, splitcharacter = "__") {   strsplit(name, splitcharacter)[[1]][1]}

#plot expression of one gene as barplot
plotgene<-function(x,n){
  barplot(as.matrix(x[grep(n,rownames(x)),]),main=n)
}

# make GENEID the rownames and remove column GENEID
mvgeneid<-function(data){
  data <- as.data.frame(data)
  rownames(data) = data[,1]
  data= data[,-1]
  return(data)
}

#reorder cells from four CS1 primer libraries into one 384 column-long library
# libraries is a vector containing the four library names, in the order A1,A2,B1,B2
reorder.cs1<-function(libraries,name){
  tc<-list()
  rc<-list()
  bc<-list()
  for(i in 1:4){
    tc[[i]]<-  read.csv(paste(libraries[i],".coutt.csv", sep=""), header = TRUE, sep = "\t",row.names =1)
    rc[[i]] <- read.csv(paste(libraries[i],".coutc.csv", sep=""), header = TRUE, sep = "\t",row.names =1) 
    bc[[i]] <- read.csv(paste(libraries[i],".coutb.csv", sep=""), header = TRUE, sep = "\t",row.names =1)
  }
  merge.tc<-intersectmatrix(tc[[1]],intersectmatrix(tc[[2]],intersectmatrix(tc[[3]],tc[[4]])))
  merge.bc<-intersectmatrix(bc[[1]],intersectmatrix(bc[[2]],intersectmatrix(bc[[3]],bc[[4]])))
  merge.rc<-intersectmatrix(rc[[1]],intersectmatrix(rc[[2]],intersectmatrix(rc[[3]],rc[[4]])))

  order<-c(matrix(c(96*0+seq(1,96), 96*1+seq(1,96)), 2, byrow = T)) 
  order2<-c(matrix(c(96*2+seq(1,96), 96*3+seq(1,96)), 2, byrow = T)) 
  all<-c()
  for(i in 0:7){
    all<-c(all,order[(1+i*24):((i+1)*24)],order2[(1+i*24):((i+1)*24)])
  }
  merge.order.tc<-merge.tc[all]
  merge.order.bc<-merge.bc[all]
  merge.order.rc<-merge.rc[all]
  
  merge.order.tc<- merge.order.tc[order(rownames( merge.order.tc)), ]
  merge.order.bc<- merge.order.bc[order(rownames( merge.order.bc)), ]
  merge.order.rc<- merge.order.rc[order(rownames( merge.order.rc)), ]
  
  write.table(merge.order.tc,paste(name,".coutt.csv",sep=""),sep="\t")
  write.table(merge.order.bc,paste(name,".coutb.csv",sep=""),sep="\t")
  write.table(merge.order.rc,paste(name,".coutc.csv",sep=""),sep="\t")
}

#JC's merge function (produces new rows on bottom of new dataframe, so reorder rows alphabetically afterwards)
intersectmatrix<-function(x,y){
  a<-setdiff(row.names(x),row.names(y))
  b<-setdiff(row.names(y),row.names(x))
  d<-matrix(data = 0,nrow = length(a),ncol = ncol(y))
  row.names(d)<-a
  colnames(d)<-colnames(y)
  c<-matrix(data = 0,nrow = length(b),ncol = ncol(x))
  row.names(c)<-b
  colnames(c)<-colnames(x)
  y<-rbind(y,d)
  x<-rbind(x,c)
  e <- match(rownames(x), rownames(y))
  f <- cbind( x, y[e,])
  return(f)
}


#overseq2, plot oversequencing per transcript
overseq2 <- function(x,y){
  main=paste("oversequencing_molecules")   # mixes string + name of choice
  xlab=bquote(log[10] ~ "read counts / barcode counts")    #  subscript in string
  rc.v<-as.vector(unlist(x))[as.vector(unlist(x>0))]  
  bc.v<-as.vector(unlist(y))[as.vector(unlist(y>0))]
  results<-rc.v/bc.v
  sub=paste("median",round(median(rc.v/bc.v),3),sep=" ")
  hist(log10(results),breaks=75, col="red", main=main,xlab=xlab,sub=sub)
}


#plot total number of reads per sample
totalreads <- function(data,plotmethod=c("barplot","hist","cumulative","combo")){
  if ( ! plotmethod %in% c("barplot","hist","cumulative","combo") ) stop("invalid method")
  if(plotmethod == "hist"){
    a<-hist(log10(colSums(data)),breaks=100,xlab="log10(counts)",ylab="frequency",main="total unique reads",col="grey",xaxt="n",col.sub="red") 
    mtext(paste("mean:",round(mean(colSums(data)))," median:",round(median(colSums(data)))),side=3,col="red",cex=0.8)
    axis(1,at=a$breaks[which(a$breaks %in% seq(0,max(a$breaks),1))],labels=a$breaks[which(a$breaks %in% c(0,1,2,3,4,5))])
    axis(1,at=a$breaks[which(a$breaks %in% seq(0,max(a$breaks),1000))],labels=a$breaks[which(a$breaks %in% seq(0,max(a$breaks),1000))])
    
    abline(v=log10(mean(colSums(data))/2),col="red")
    text(log10(mean(colSums(data))/2),max(a$counts)-2, round(mean(colSums(data))/2), srt=0.2, col = "red",pos=2)
  }
  if(plotmethod == "barplot"){
    b<-barplot(colSums(data),xaxt="n",xlab="cells",sub=paste("mean total read:",round(mean(colSums(data)))),main="total unique reads",col="black",border=NA) 
    axis(1,at=b,labels=c(1:length(data))) # 1=horizontal at = position of marks
    abline(h=mean(colSums(data)),col="red")
  }
  if(plotmethod == "cumulative"){
    plot(ecdf(colSums(data)),xlab="total reads",ylab="fraction",main="total unique reads",col="red",tck=1,pch=19,cex=0.5,cex.axis=0.8) 
    abline(v=mean(colSums(data)/2),col="red")
    mtext(paste("mean:",round(mean(colSums(data)))," median:",round(median(colSums(data)))),side=3,col="red",cex=0.8)
  }
  
  if(plotmethod == "combo"){
    a<-hist(log10(colSums(data)),breaks=100,xlab="log10(counts)",ylab="frequency",main="total unique reads",col="grey",xaxt="n",col.sub="red") 
    mtext(paste("mean:",round(mean(colSums(data)))," median:",round(median(colSums(data)))),side=3,col="red",cex=0.8)
    axis(1,at=a$breaks[which(a$breaks %in% c(0,1,2,3,4,5))],labels=a$breaks[which(a$breaks %in% c(0,1,2,3,4,5))])
    abline(v=log10(mean(colSums(data))/2),col="red")
    text(log10(mean(colSums(data))/2),max(a$counts)-2, round(mean(colSums(data))/2), srt=0.2, col = "red",pos=2)
    plotInset(log10(1),max(a$counts)/4,log10(250), max(a$counts),mar=c(1,1,1,1),
              plot(ecdf(colSums(data)),pch=".",col="red",cex=0.5,ylab=NA,xlab=NA,main=NA,cex.axis=0.8,xaxt="n",las=3,mgp=c(2,0.1,0),tck=1,bty="n"),
              debug = getOption("oceDebug"))
  }
}


#plot amount of genes detected per cell
cellgenes<-function(data,plotmethod=c("hist","cumulative","combo")){
  if ( ! plotmethod %in% c("hist","cumulative","combo") ) stop("invalid plotting method")
    genes<-apply(data,2,function(x) sum(x>=1))
  if(plotmethod == "hist"){
    a<-hist(genes,breaks=100,xlab="total genes",ylab="frequency",main="detected genes/cell",col="steelblue1",xaxt="n") 
    mtext(paste("mean:",round(mean(genes))," median:",round(median(genes))),side=3,col="red",cex=0.8)
    axis(1,at=a$breaks[which(a$breaks %in% seq(0,max(a$breaks),1000))],labels=a$breaks[which(a$breaks %in% seq(0,max(a$breaks),1000))])
  }
  if(plotmethod == "cumulative"){
    plot(ecdf(genes),pch=19,col="red",cex=0.5,ylab="frequency",xlab="detected genes/cell",main="cumulative dist genes",cex.axis=1,las=1,tck=1)
    mtext(paste("mean:",round(mean(genes))," median:",round(median(genes))),side=3,col="red",cex=0.8)
  }
  if(plotmethod == "combo"){
    a<-hist(genes,breaks=100,xlab="log10(counts)",ylab="frequency",main="detected genes/cell",col="steelblue1",xaxt="n") 
    mtext(paste("mean:",round(mean(genes))," median:",round(median(genes))),side=3,col="red",cex=0.8)
    axis(1,at=a$breaks[which(a$breaks %in% seq(0,max(a$breaks),1000))],labels=a$breaks[which(a$breaks %in% seq(0,max(a$breaks),1000))])
    plotInset(max(genes)/3,max(a$counts)/3,max(genes), max(a$counts),mar=c(1,1,1,1),
              plot(ecdf(colSums(data)),pch=19,col="red",cex=0.5,ylab=NA,xlab=NA,main=NA,cex.axis=0.6,las=3),
              debug = getOption("oceDebug"))
  }
}

#plot ERCC reads
plotspike<-function(data){
  erccs<-data[grep("ERCC-",rownames(data)),]
    b<-barplot(colSums(erccs),main="ERCC reads",ylab="total ERCC reads",xlab="cells",col="orange",xaxt="n",border=NA)
    axis(1,at=b,labels=c(1:length(data)))
}

#plot number of available transcripts vs cutoffs of median detected transcripts
testcutoff<-function(data,n,pdf=FALSE){
  main=paste("genes cutoff test",n)
  for(l in 1:15){
    z = apply(data,1,median) > l
    if(l==1){
      rc.cutoff = z
    } else {
       rc.cutoff = cbind(rc.cutoff,z)
    }
  }
  if (pdf){
    pdf(paste(getwd(),main,".pdf",sep="")) 
    plot(apply(rc.cutoff,2,sum),ylab = "number of transcripts",col="black",
    xlab = "cutoff (mean transcript no.)",main=main,type="b",lty=2,pch=19)
    dev.off()
  }
  else{
    plot(apply(rc.cutoff,2,sum),ylab = "number of transcripts",col="black",
    xlab = "cutoff (mean transcript no.)",main=main,type="b",lty=2,pch=19)
  }    
}
    

#plot number of total reads, ERCC-reads and genes/cell over a 384-well plate layout
plate.plots<-function(data){
  # genes<-apply(data,2,function(x) sum(x>=1))# calculate detected genes/cell
  spike<-colSums(keepspike(data))+0.1
  # calculate sum of spike in per cell
  total<-colSums(rmspike(data+0.1)) # sum of unique reads after removing spike ins
  palette <- colorRampPalette(rev(brewer.pal(n = 11,name = "RdYlBu")))(10) # pick which palette for plate plotting
  coordinates<-expand.grid(seq(1,24),rev(seq(1,16)))
  plot(expand.grid(x = c(1:24), y = c(1:16)),main="Unique non ERCC reads",ylab=NA,xlab=NA) #plate layout
  mtext(paste(">1500 unique reads :",round(length(which(colSums(data)>1500))/384*100),"%"),col="red",cex=0.9)
  points(coordinates,pch=19,col=palette[cut(log10(total),10)]) # plot total non-ERCC reads/cell over layout

  plot(expand.grid(x = c(1:24), y = c(1:16)),main="sum of all ERCCs",ylab=NA,xlab=NA) #plate layout
  points(coordinates,pch=19,col=palette[cut(log10(spike),10)]) #plot sum of spike ins over plate
  mtext(paste(">100 ERCCs :",round(length(which(colSums(keepspike(data))>100))/384*100),"%"),col="red",cex=0.9)
  
  plot(expand.grid(x = c(1:24), y = c(1:16)),main="sum ERCC/sum non ERCC reads",ylab=NA,xlab=NA) 
  points(coordinates,pch=19,col=palette[cut(spike/total,10)]) #plot ERCC reads/non-ERCC reads/cell
  mtext(paste(">10% spike in reads:",round(length(which(spike/total>0.05))/384*100),"%"),col="red",cex=0.9)
  
}

# plot the top 20 genes with expresion bar and then a CV plot for the same genes
topgenes<-function(data){
data<-rmspike(data)
means<-apply(data,1,mean)
vars<-apply(data,1,var)
cv<-vars/means
means<-means[order(means, decreasing = TRUE)]
cv<-cv[order(cv, decreasing = TRUE)]
names(means)<-sapply(names(means),chop_chr)
names(cv)<-sapply(names(cv),chop_chr)
barplot(log2(rev(means[1:20])),las=1,cex.names = 0.5, main="top expressed genes", xlab="log2(mean expression)",horiz=TRUE)
barplot(log2(rev(cv[1:20])),las=1,cex.names = 0.5, main="top noisy genes",xlab="log2(var/mean)",horiz=TRUE)
}


#Read files in specified directory automatically (based on Thoms script)
read_files <- function(dir = "", name = Sys.Date()){
  
  #add "/" to dir
  if(substr(dir, start = nchar(dir), stop = nchar(dir)) != "/" && dir != ""){
    dir <- paste(dir, "/", sep = "")
  }
  
  #Read files
  files <- list.files(dir, ".cout(t|b|c).csv")
  split <- strsplit(files,split = ".cout")
  file_names <- unique(as.character(data.frame(split, stringsAsFactors = FALSE)[1,]))
  
  #This check if all necessary files are in the script
  error <- ""
  for(i in 1:length(file_names)){
    
    if(file.exists(paste(dir, file_names[i],".coutb.csv", sep="")) == FALSE){
      f <- paste(file_names[i], ".coutb.csv", " is not found!", sep = "")
      error <- paste(error, "\n", f)
    }
    
    if(file.exists(paste(dir, file_names[i],".coutc.csv", sep="")) == FALSE){
      f <- paste(file_names[i], ".coutc.csv", " is not found!", sep = "")
      error <- paste(error, "\n", f)
    }
    
    if(file.exists(paste(dir,file_names[i],".coutt.csv", sep="")) == FALSE){
      f <- paste(file_names[i], ".coutt.csv", " is not found!", sep = "")
      error <- paste(error, "\n", f)
    }
  }
  
  if(error != ""){
    stop(error)
  }
  cat("the following plates will be processed:\n")
  print(file_names)

  output <- paste(dir,file_names, sep="")
  return(output)
}

# check expression in empty corner of plate and calculate "leakyness" from highly expressed genes
leakygenes<-function(data){
  corner<-data[emptywells] # subset data to 8 wells specified in diagnotics script as empty corner
  names(corner)<-c("O21","O22","O23","O24","P21","P22","P23","P24")
  genes<-apply(data,2,function(x) sum(x>=1)) # check how many genes are detected
  genes.corner<-apply(rmspike(corner),2,function(x) sum(x>=1)) # remove ERCC reads
  spike.corner<-colSums(keepspike(corner)) # keep only ERCC reads
  genespike<-data.frame(genes=genes.corner,ERCC=spike.corner)
  if(length(which(genes.corner > mean(genes/5))) != 0){
    stop(paste("Not all 8 corner samples are empty in", names[[i]],": won't be plotted"))
  } else {# check if the corner wells were actually empty, otherwise stop
        # plot genes/cell and ERCC reads/cell for corner wells
    par(mar = c(5, 4, 6, 1))
    barplot(t(genespike),main="total genes and ERCCs \n in empty corner",
          col=c("blue","red"),space=rep(c(0.7,0),8),cex.names = 0.8,las=3,beside=TRUE,
          legend=colnames(genespike),args.legend = list(x = "topright", bty = "n",horiz=TRUE,inset=c(0,-0.25)))
    }
  # determine top expressed genes in corner and compare to mean expressed genes in plate
  if( length(which(spike.corner > 75)) == 0){
    stop(paste("There are no samples with more than 75 ERCC reads in", names[[i]]))
  }  
  cornerz<-corner[which(spike.corner>75)]  # take only wells which worked (>75 ERCC reads)
  cornerz<-rmspike(cornerz) # remove ERCCs
  mean.corner<-apply(cornerz,1,sum)[order(apply(cornerz,1,sum),decreasing=TRUE)][1:50] # pick top 50 in corner
  mean.all<-apply(data,1,sum)[order(apply(data,1,sum),decreasing=TRUE)][1:200] # pick top 200 in plate
  names(mean.corner)<-sapply(names(mean.corner),chop_chr) # remove __chr* from name
  names(mean.all)<-sapply(names(mean.all),chop_chr) # remove __chr* from name
  overlap<-mean.corner[names(mean.corner) %in% names(mean.all)] # check overal between top 50 corner and 200 in plate
  non.overlap<-mean.corner[!names(mean.corner) %in% names(mean.all)]
  b<-barplot(log2(rev(overlap[1:10])),las=1,cex.names = 0.6, main="top 10 overlapping genes",sub="barcode leaking in %", xlab="log2(sum of reads in corner)",horiz=TRUE)
  text(0.5,b, round((mean.corner[names(overlap)[1:10]]/mean.all[names(overlap)[1:10]])*100,2))
  if (length(overlap)==50){
    warning(paste("there is complete overlap between corner genes and plate genes in ", names[[i]]))
  }
  else{
    barplot(log2(rev(non.overlap[1:length(non.overlap)])),las=1,cex.names = 0.6, main="top 50 empty corner genes \n not in top 200 plate genes", xlab="log2(mean expression)",horiz=TRUE)
  }
}