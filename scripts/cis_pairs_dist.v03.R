#! /software/crgadm/software/R/4.2.0-foss-2021b/bin/Rscript

### USAGE CNAG cluster
##  module purge; module load R/4.2.0-foss-2021b; 
##  Usage: cis_pairs_dist.v03.R <infile.pairs> <out_basename>
##  example: cis_pairs_dist.v03.R mapped.mq40.Capra_hircus.pairs mRupRup

# Name: cis_pairs_dist.v03.R
# date: 2025-10-16

# NOTES: Based on an R recipe provided by Dovetail Genomics' bioinformatic support


# Pass parameters via command line
args <- commandArgs(trailingOnly = TRUE) ##  one argument only -  start number 

pairs_file <- args[1]; # .pairs file contains columns with scaffold/contig ID and position for each member of the read pair
prefix <- args[2]; #often the ToLID and/or library name and/or million read pairs 

pdfName<-paste(prefix,".Cis_pairs_histogram.pdf",sep="")

##Make Pair Distance frequency plot##
       
       #import pairs file (this may take some time depending on sequencing depth)
       	    d <- read.table(pairs_file, header=F) # this matches the total No-Dup Read Pairs in Rupicapra 2,105,550
	       	 
	    #Select cis pairs
		c <- d[which(d$V2 == d$V4),] # nrow(c) is equal to No-Dup Cis Read Pairs
			   
	        #Add column for distance between reads of each pair
	            c$V9 <- (c$V5 - c$V3)

#Store the histgram class to detect the maximal frequency
hist_val<-(hist(log10(c$V9),breaks=1000, plot=FALSE)) 
	# Note that hist automatically plots inside a pdf unless disabled
	# plot = FALSE tells R to compute the histogram bins and counts without actually plotting

#summary(hist_val$counts)

pdf(file=pdfName) #  open pdf file

 par(bg="white")# set the background color as white, default is transparent
 
 #without x axis xaxt="n"
 hist(log10(c$V9), breaks=1000,main="Histogram Cis Pairs Distance hic_qc Test",xlab=("log10(Pair Distance bp)"),xaxt="n")
 axis(side=1, at=seq(0, max(log10(c$V9)), 1))
 abline(v=3,col="red",lwd=3) #10^3 is 1,000 bp = 1Kb
 text(2, max(hist_val$counts), '<1 Kb',col="red",cex=1.5)
 text(4, max(hist_val$counts), '>1 Kb',col="darkgreen",cex=1.5)
  
dev.off()



