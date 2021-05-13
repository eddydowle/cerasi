#Eddy October 2016

#source("https://bioconductor.org/biocLite.R")
#biocLite("edgeR")
#biocLite("locfit")

#this is the orginal cerasi edgeR comparisons. However there is three samples at 3M low and one 2.5 at 3M low that are really funky looking so there is also a edgeR comparison with those out.


library(edgeR)
library(stringr)
library(locfit)
library(statmod)
library(plyr)
library(ggplot2)
library(reshape)
library(heatmap3)
library(dplyr)

#setwd("/media/raglandlab/ExtraDrive2/WinterLengthRnaSeq/RawData/HO_PE100_20160112_1mismatch/mappedcounts")
setwd("~/Documents/Cerasi/Cerosi/salmonResults2018MainTranscript/")

#all cerosi work is now done with file that have had funny characters that end the sequence name removed, causing too many issues

warnings()
file_list<-list.files(pattern='*sf.cleanedNames')
file_list
table<-c()

for (file in file_list){
  sampleID<-str_match(file,"(^[a-z]*_[A-Za-z]*_[0-9]*_[A-Za-z0-9]*_[0-9]*)")[,1]
  if (is.na(sampleID)){
    sampleID<-str_match(file,"(^[a-z]*_[A-Za-z]*_[A-Z0-9]*_[0-9]*)")[,1]
  }
  print(sampleID)
  df<-read.table(file,header=T,sep="\t",row.names=NULL)
  df2<-df[,c(1,5)]
  colnames(df2)[colnames(df2)=="NumReads"] <- paste(sampleID,"_exp_count",sep="")
  if (length(table) ==0){table<-df2
  } else {
    table<-merge(table,df2,by="Name") }
}  

#cerosi.annotate<-read.table("../annotating/cerosi.annotated.pomonella.flybase.swissport",sep="\t",header=TRUE,row.names=NULL,stringsAsFactors = FALSE,quote = "",fill=TRUE)
#colnames(cerosi.annotate)[colnames(cerosi.annotate)=="X1transcript_id" ]<- "Name"

#this one is flybase,swiss and trembl e-6
cerosi.annotate<-read.table("../annotating/cerosi.annotated.pomonella.flybase.swissporte6.tophit.trembl.cleanedNames",sep="\t",header=TRUE,row.names=NULL,stringsAsFactors = FALSE,quote = "",fill=TRUE)
#count(unique(cerosi.annotate$Name))
#table1<-left_join(table,cerosi.annotate,by="Name")
table1<-right_join(cerosi.annotate,table,by="Name")

#tail(table1)
#head(table1)

#sum(is.na(table1$c_Hi_2_5M_01_exp_count))
#test<-table1[is.na(table1$c_Hi_2_5M_01_exp_count),]
#test2<-table1[!is.na(table1$c_Hi_2_5M_01_exp_count),]
#head(test$Name)
#tail(test$Name)
#tail(test)

sum(table1$c_Hi_3M_02_exp_count) #
sum(table1$c_Hi_3M_03_exp_count) #
sum(table1$c_Hi_3M_04_exp_count) #
sum(table1$c_Hi_3M_05_exp_count) #

sum(table1$c_Hi_2M_01_exp_count) #
sum(table1$c_Hi_2M_02_exp_count) #
sum(table1$c_Hi_2M_03_exp_count) #
sum(table1$c_Hi_2M_04_exp_count) #

sum(table1$c_Lo_3M_01_exp_count) #
sum(table1$c_Lo_3M_02_exp_count) #
sum(table1$c_Lo_3M_03_exp_count) #
sum(table1$c_Lo_3M_04_exp_count) #


#########greg's filtering function represented in >= 50% of samples by at least one count#######
filterMinCount<- function(x) {
  pres<-x >=1
  out=F
  if ((sum(pres)/length(pres)) >= 0.5) {out=T}
  return(out)
}



#standard
colnames(table1)
filterInd<-apply(table1[,(-1:-38)],1,filterMinCount)
#table <- table[filterInd,] #
table1 <- table1[filterInd,] #
tail(table1)
colnames(table1[,(-1:-38)])




sum(table1$c_Hi_3M_02_exp_count) #
sum(table1$c_Hi_3M_03_exp_count) #
sum(table1$c_Hi_3M_04_exp_count) #
sum(table1$c_Hi_3M_05_exp_count) #

sum(table1$c_Hi_2M_01_exp_count) #
sum(table1$c_Hi_2M_02_exp_count) #
sum(table1$c_Hi_2M_03_exp_count) #
sum(table1$c_Hi_2M_04_exp_count) #

sum(table1$c_Lo_3M_01_exp_count) #
sum(table1$c_Lo_3M_02_exp_count) #  
sum(table1$c_Lo_3M_03_exp_count) #
sum(table1$c_Lo_3M_04_exp_count) #

colnames(table1)
class(table1$c_Hi_2_5M_01_exp_count)
sum(is.na(table1$c_Hi_2_5M_01_exp_count))
grouping<-c(rep('Lo25M',4),rep('Lo2M',4),rep('Lo35M',4),rep('Lo3M',4),rep('Lo45M',4),rep('Lo4M',4),rep('Hi25M',4),rep('Hi2M',4),rep('Hi35M',4),rep('Hi3M',4),rep('Hi45M',4),rep('Hi4M',4))
grouping

colnames(table1)
table.dge<-DGEList(counts=table1[,39:86],genes=table1[,(1:38)], group=grouping)
table.dge<- calcNormFactors(table.dge)
table.dge$samples

keep<-rowSums(cpm(table.dge)>1) >=2
y2 <- table.dge[keep, , keep.lib.sizes=FALSE]
nrow(y2$counts) 

keep <- filterByExpr(table.dge)
#pre-filtering we have
nrow(table.dge$counts)
y <- table.dge[keep, , keep.lib.sizes=FALSE]
nrow(y$counts)

########create DGElist object#######

#standard

#table.dge<-DGEList(counts=table[,39:86],genes=table[,(1:38)])
#table.dge<- calcNormFactors(table.dge)
#table.dge$samples

#colnames(table.Filtered)

#table.high.dge<-DGEList(counts=table[,29:52],genes=table[,(1:28)])
#table.high.dge$samples
#table.low.dge<-DGEList(counts=table[,53:76],genes=table[,(1:28)])
#table.low.dge$samples


#table.high.dge<-calcNormFactors(table.high.dge)
#table.low.dge<-calcNormFactors(table.low.dge)

plotMDS(table.dge)
colnames(table.dge)

plotMDS(table.dge) 
plotMDS(table.dge, top=500,pch = c(rep(15,4),rep(2,4),rep(4,4),rep(19,4),rep(7,4),rep(5,4),rep(15,4),rep(2,4),rep(4,4),rep(19,4),rep(7,4),rep(5,4)),col=c(rep("green4",24), rep("purple4",24)))
#bottom right
legend(1.3,0.15,bty = "n",legend=(c("2 mo.","2.5 mo.","3 mo.","3.5 mo.","4 mo.","4.5 mo.","Hi","Low")),pch = c(2,15,19,4,5,7,18,18),col=c("black","black","black","black","black","black","green4","purple4"))

#plotMDS(table.high.dge)
#high
#plotMDS(table.high.dge, top=500,pch = c(rep(15,4),rep(2,4),rep(4,4),rep(19,4),rep(7,4),rep(5,4)),col=c(rep("blue",24)))
#bottom right
#legend(1.7,-1.1,bty = "n",legend=(c("2 mo.","2.5 mo.","3 mo.","3.5 mo.","4 mo.","4.5 mo.")),pch = c(2,15,19,4,5,7),col=c(rep("blue",6)))

#plotMDS(table.low.dge)
#low
#plotMDS(table.low.dge, top=500,pch = c(rep(15,4),rep(2,4),rep(4,4),rep(19,4),rep(7,4),rep(5,4)),col=c(rep("red",24)))
#top left
#legend(-3.6,2.0,bty = "n",legend=(c("2 mo.","2.5 mo.","3 mo.","3.5 mo.","4 mo.","4.5 mo.")),pch = c(2,15,19,4,5,7),col=c(rep("red",6)))

########cerosi design 1########


month<-factor(c(rep("2_5M",4),rep("2M",4),rep("3_5M",4),rep("3M",4),rep("4_5M",4),rep("4M",4),rep("2_5M",4),rep("2M",4),rep("3_5M",4),rep("3M",4),rep("4_5M",4),rep("4M",4)))
altitute<-factor(c(rep("High",24),rep("Low",24)))

table.sampleData<-data.frame(Sample=colnames(table.dge),month,altitute)

#model
levels(month)
month_<-factor(month, levels = c("2M", "2_5M", "3M","3_5M","4M","4_5M"))
#altitute_<-factor(altitute,levels = c ("Low","High"))
levels(month_)
levels(altitute)
design <- model.matrix(~month_*altitute)
rownames(design) <- colnames(table.dge)
design

#                        (Intercept)	month_2_5M	month_3M	month_3_5M	month_4M	month_4_5M	altituteLow	month_2_5M:altituteLow	month_3M:altituteLow	month_3_5M:altituteLow	month_4M:altituteLow	month_4_5M:altituteLow
#c_Hi_2_5M_01_exp_count	1		          1	        	0	    		0		    		0	  			0				  	0		      	0			    	    				0				    					0								  			0					  					0
#c_Hi_2_5M_02_exp_count	1		        	1		    		0		    	0		  			0					0			  		0	    			0				    						0											0						  					0											0
#c_Hi_2_5M_04_exp_count	1		        	1		    		0		    	0		  			0					0			  		0	    			0				    						0											0						  					0											0
#c_Hi_2_5M_05_exp_count	1		        	1		    		0		    	0		  			0					0			  		0	    			0				    						0											0						  					0											0
#c_Hi_2M_01_exp_count  	1		        	0		    		0		    	0		  			0					0			  		0		    		0				    						0											0						  					0											0
#c_Hi_2M_02_exp_count  	1	        		0		    		0		    	0		  			0					0			  		0	    			0				    						0											0							  				0											0
#c_Hi_2M_03_exp_count  	1		        	0		    		0	    		0		  			0					0			  		0		    		0				    						0											0						  					0											0
#c_Hi_2M_04_exp_count  	1		        	0		    		0		    	0	  				0					0			  		0	    			0				    						0											0						  					0											0
#c_Hi_3_5M_01_exp_count	1		        	0		    		0		    	1		  			0					0			  		0	    			0				    						0											0						  					0											0
#c_Hi_3_5M_02_exp_count	1		        	0		    		0	    		1	  				0					0			  		0	    			0				    						0											0						  					0											0
#c_Hi_3_5M_03_exp_count	1		        	0		    		0	    		1		  			0					0			  		0	    			0				    						0											0						  					0											0
#c_Hi_3_5M_05_exp_count	1		        	0		    		0	    		1		  			0					0			  		0	    			0				    						0											0						  					0											0
#c_Hi_3M_02_exp_count  	1		        	0		    		1	    		0		  			0					0		  			0	    			0			    							0											0						  					0											0
#c_Hi_3M_03_exp_count  	1		        	0		    		1	    		0		  			0					0			  		0		    		0			    							0											0						  					0											0
#c_Hi_3M_04_exp_count	  1		        	0		    		1	    		0		  			0					0			  		0		    		0			    							0											0						  					0											0
#c_Hi_3M_05_exp_count	  1		        	0		    		1	    		0	  				0					0			  		0		    		0			    							0											0						  					0											0
#c_Hi_4_5M_02_exp_count	1		        	0		    		0	    		0		  			0					1			  		0		    		0			    							0											0						  					0											0
#c_Hi_4_5M_03_exp_count	1	        		0	    			0	    		0	  				0					1		  			0		    		0			    							0											0						  					0											0
#c_Hi_4_5M_04_exp_count	1	        		0	    			0	    		0	  				0					1		  			0		    		0			    							0											0					  						0											0
#c_Hi_4_5M_05_exp_count	1	        		0	    			0	    		0		  			0					1			  		0		    		0			    							0											0						  					0											0
#c_Hi_4M_01_exp_count	  1	        		0	    			0	    		0		  			1					0			  		0		    		0			    							0											0						  					0											0
#c_Hi_4M_02_exp_count	  1		        	0	    			0	    		0		  			1					0			  		0	    			0			    							0											0						  					0											0
#c_Hi_4M_03_exp_count	  1	        		0	    			0	    		0		  			1					0			  		0	    			0			    							0											0						  					0											0
#c_Hi_4M_05_exp_count	  1	        		0		    		0		    	0		  			1					0		  			0		    		0			    							0											0					  						0											0
#c_Lo_2_5M_01_exp_count	1	        		1		    		0		    	0		  			0					0		  			1		    		1			    							0											0					  						0											0
#c_Lo_2_5M_02_exp_count	1	        		1		    		0		    	0		  			0					0		  			1	    			1			    							0											0					  						0											0
#c_Lo_2_5M_03_exp_count	1	        		1		    		0		    	0  		  		0					0		  			1	    			1			    							0											0					  						0											0
#c_Lo_2_5M_04_exp_count	1	        		1		    		0		    	0		  			0					0		  			1	    			1			    							0											0					  						0											0
#c_Lo_2M_01_exp_count	  1	        		0	    			0		    	0		  			0					0		  			1	    			0			    							0											0					  						0											0
#c_Lo_2M_02_exp_count	  1	        		0	    			0		    	0		  			0					0			  		1	    			0			    							0											0					  						0											0
#c_Lo_2M_03_exp_count	  1	        		0	    			0		    	0		  			0					0			  		1	    			0			    							0											0					  						0											0
#c_Lo_2M_04_exp_count	  1	        		0	    			0		    	0		  			0					0			  		1		    		0			    							0											0					  						0											0
#c_Lo_3_5M_01_exp_count	1	        		0		    		0		    	1		  			0					0			  		1	    			0				    						0											1					  						0											0
#c_Lo_3_5M_02_exp_count	1	        		0		    		0		    	1		  			0					0			  		1	    			0				    						0											1					  						0											0
#c_Lo_3_5M_03_exp_count	1	        		0		    		0		    	1		  			0					0		  			1	    			0				    						0											1					  						0											0
#c_Lo_3_5M_04_exp_count	1		        	0		    		0		    	1		  			0					0		  			1	    			0				    						0											1					  						0											0
#c_Lo_3M_01_exp_count	  1		        	0		    		1		    	0		  			0					0		  			1	    			0				    						1											0					  						0											0
#c_Lo_3M_02_exp_count	  1	        		0		    		1		    	0		  			0					0		  			1		    		0				    						1											0					  						0											0
#c_Lo_3M_03_exp_count	  1	        		0		    		1	    		0		  			0					0		  			1	    			0			    							1											0					  						0											0
#c_Lo_3M_04_exp_count	  1	        		0	    			1	    		0		  			0					0		  			1	    			0				    						1											0					  						0											0
#c_Lo_4_5M_01_exp_count	1	        		0		    		0		    	0		  			0					1		  			1	    			0				    						0											0						  					0											1
#c_Lo_4_5M_02_exp_count	1	        		0		    		0		    	0	  				0					1			  		1		    		0			    							0											0						  					0											1
#c_Lo_4_5M_03_exp_count	1	        		0	    			0	    		0	  				0					1		  			1		    		0			    							0											0						  					0											1
#c_Lo_4_5M_04_exp_count	1	        		0		    		0	    		0		  			0					1		  			1	    			0			    							0											0						  					0											1
#c_Lo_4M_01_exp_count	  1	        		0		    		0	    		0		  			1					0		  			1	    			0			    							0											0						  					1											0
#c_Lo_4M_02_exp_count	  1	        		0		    		0		    	0		  			1					0		  			1	    			0				    						0											0						  					1											0
#c_Lo_4M_03_exp_count	  1	        		0		    		0		    	0		  			1					0		  			1	    			0				    						0											0						  					1											0
#c_Lo_4M_04_exp_count	  1		        	0		    		0		    	0		  			1					0		  			1	    			0				    						0											0						  					1											0


#########estimate disperson squareroot=coefficient of variation of biological variation#######
table.dge<- estimateDisp(table.dge, design, robust=TRUE)
table.dge$common.dispersion
#0.05539425
sqrt(0.05539425)
# 0.2353598 coefficent of variation


#plot
plotBCV(table.dge)

#can look at trended dispersion that is used for QL (Quasi-Likelihood pipeline)
fit <- glmFit(table.dge, design)
colnames(design)

#fit <- glmQLFit(table.dge, design)
#plotQLDisp(fit)


##############contrasts##################

fit <- glmFit(table.dge, design)


colnames(design)
#[1] "(Intercept)"              "month_2_5M"               "month_3M"                 "month_3_5M"               "month_4M"                
#[6] "month_4_5M"               "altitute_High"            "month_2_5M:altitute_High" "month_3M:altitute_High"   "month_3_5M:altitute_High"
#[11] "month_4M:altitute_High"   "month_4_5M:altitute_High"

#2_5 to two months low
lrt2_5monthvs2low <- glmLRT(fit,contrast=c(0,1,0,0,0,0,0,1,0,0,0,0))
Low2_5monthvs2 <- data.frame(topTags(lrt2_5monthvs2low,n=nrow(table),sort="none"))
#back to 2M i.e. 2M->2_5M 2M->3M 2M->3_5M 2M->4M 2M->4_5M 
FDR.table.Lowvs2<-Low2_5monthvs2[,c(1,43)]
colnames(FDR.table.Lowvs2)[colnames(FDR.table.Lowvs2)=="FDR"] <- "month2_5vs2_low_FDR"
LogFC.table.Lowvs2<-Low2_5monthvs2[,c(1,39)]
colnames(LogFC.table.Lowvs2)[colnames(LogFC.table.Lowvs2)=="logFC"] <- "month2_5vs2_low_logFC"
#sequential i.e. 2M->2_5M->3M->3_5M->4M->4_5M
FDR.table.LowSeq<-Low2_5monthvs2[,c(1,43)]
colnames(FDR.table.LowSeq)[colnames(FDR.table.LowSeq)=="FDR"] <- "month2_5vs2_low_FDR"
LogFC.table.LowSeq<-Low2_5monthvs2[,c(1,39)]
colnames(LogFC.table.LowSeq)[colnames(LogFC.table.LowSeq)=="logFC"] <- "month2_5vs2_low_logFC"

#3 to two months low
lrt3monthvs2low <- glmLRT(fit,contrast=c(0,0,1,0,0,0,0,0,1,0,0,0))
Low3monthvs2 <- data.frame(topTags(lrt3monthvs2low,n=nrow(table),sort="none"))
FDR.table.Lowvs2$month3vs2_low_FDR<-Low3monthvs2[,43]
LogFC.table.Lowvs2$month3vs2_low_logFC<-Low3monthvs2[,39]

#3_5 to two months low
lrt3_5monthvs2low <- glmLRT(fit,contrast=c(0,0,0,1,0,0,0,0,0,1,0,0))
Low3_5monthvs2 <- data.frame(topTags(lrt3_5monthvs2low,n=nrow(table),sort="none"))
FDR.table.Lowvs2$month3_5vs2_low_FDR<-Low3_5monthvs2[,43]
LogFC.table.Lowvs2$month3_5vs2_low_logFC<-Low3_5monthvs2[,39]

#4 to two months low
lrt4monthvs2low <- glmLRT(fit,contrast=c(0,0,0,0,1,0,0,0,0,0,1,0))
Low4monthvs2 <- data.frame(topTags(lrt4monthvs2low,n=nrow(table),sort="none"))
FDR.table.Lowvs2$month4vs2_low_FDR<-Low4monthvs2[,43]
LogFC.table.Lowvs2$month4vs2_low_logFC<-Low4monthvs2[,39]

#4_5 to two months low
lrt4_5monthvs2low <- glmLRT(fit,contrast=c(0,0,0,0,0,1,0,0,0,0,0,1))
Low4_5monthvs2 <- data.frame(topTags(lrt4_5monthvs2low,n=nrow(table),sort="none"))
FDR.table.Lowvs2$month4_5vs2_low_FDR<-Low4_5monthvs2[,43]
LogFC.table.Lowvs2$month4_5vs2_low_logFC<-Low4_5monthvs2[,39]

#3 to 2_5 months low
lrt3monthvs2_5low <- glmLRT(fit,contrast=c(0,-1,1,0,0,0,0,-1,1,0,0,0))
Low3monthvs2_5 <- data.frame(topTags(lrt3monthvs2_5low,n=nrow(table),sort="none"))
FDR.table.LowSeq$month3vs2_5_low_FDR<-Low3monthvs2_5[,43]
LogFC.table.LowSeq$month3vs2_5_low_logFC<-Low3monthvs2_5[,39]

#3_5 to 2_5 months low
lrt3_5monthvs2_5low <- glmLRT(fit,contrast=c(0,-1,0,1,0,0,0,-1,0,1,0,0))

#4 to 2_5 months low
lrt4monthvs2_5low <- glmLRT(fit,contrast=c(0,-1,0,0,1,0,0,-1,0,0,1,0))

#4_5 to 2_5 months low
lrt4_5monthvs2_5low <- glmLRT(fit,contrast=c(0,-1,0,0,0,1,0,-1,0,0,0,1))

#3_5 to 3 months low
lrt3_5monthvs3low <- glmLRT(fit,contrast=c(0,0,-1,1,0,0,0,0,-1,1,0,0))
Low3_5monthvs3 <- data.frame(topTags(lrt3_5monthvs3low,n=nrow(table),sort="none"))
FDR.table.LowSeq$month3_5vs3_low_FDR<-Low3_5monthvs3[,43]
LogFC.table.LowSeq$month3_5vs3_low_logFC<-Low3_5monthvs3[,39]

#4 to 3 months low
lrt4monthvs3low <- glmLRT(fit,contrast=c(0,0,-1,0,1,0,0,0,-1,0,1,0))

#4_5 to 3 months low
lrt4_5monthvs3low <- glmLRT(fit,contrast=c(0,0,-1,0,0,1,0,0,-1,0,0,1))

#4 to 3_5 months low
lrt4monthvs3_5low <- glmLRT(fit,contrast=c(0,0,0,-1,1,0,0,0,0,-1,1,0))
Low4monthvs3_5 <- data.frame(topTags(lrt4monthvs3_5low,n=nrow(table),sort="none"))
FDR.table.LowSeq$month4vs3_5_low_FDR<-Low4monthvs3_5[,43]
LogFC.table.LowSeq$month4vs3_5_low_logFC<-Low4monthvs3_5[,39]

#4_5 to 3_5 months low
lrt4_5monthvs3_5low <- glmLRT(fit,contrast=c(0,0,0,-1,0,1,0,0,0,-1,0,1))

#4_5 to 4 months low
lrt4_5monthvs4low <- glmLRT(fit,contrast=c(0,0,0,0,-1,1,0,0,0,0,-1,1))
Low4_5monthvs4 <- data.frame(topTags(lrt4_5monthvs4low,n=nrow(table),sort="none"))
FDR.table.LowSeq$month4_5vs4_low_FDR<-Low4_5monthvs4[,43]
LogFC.table.LowSeq$month4_5vs4_low_logFC<-Low4_5monthvs4[,39]

#high to low
#two months High vs Low
lrt2monthHighvsLow <- glmLRT(fit, contrast = c(0,0,0,0,0,0,1,0,0,0,0,0))
LowvsHigh2months <- data.frame(topTags(lrt2monthHighvsLow,n=nrow(table),sort="none"))
FDR.table.LowvsHigh<-LowvsHigh2months[,c(1,43)]
colnames(FDR.table.LowvsHigh)[colnames(FDR.table.LowvsHigh)=="FDR"] <- "month2lowvs2high_FDR"
LogFC.table.LowvsHigh<-LowvsHigh2months[,c(1,39)]
colnames(LogFC.table.LowvsHigh)[colnames(LogFC.table.LowvsHigh)=="logFC"] <- "month2lowvs2high_logFC"

#2_5 months High vs Low
lrt2_5monthHighvsLow <- glmLRT(fit, contrast = c(0,0,0,0,0,0,1,1,0,0,0,0))
LowvsHigh2_5months <- data.frame(topTags(lrt2_5monthHighvsLow,n=nrow(table),sort="none"))
FDR.table.LowvsHigh$month2_5lowvs2_5high_FDR<-LowvsHigh2_5months[,43]
LogFC.table.LowvsHigh$month2_5lowvs2_5high_logFC<-LowvsHigh2_5months[,39]

#3 months High vs Low
lrt3monthHighvsLow <- glmLRT(fit, contrast = c(0,0,0,0,0,0,1,0,1,0,0,0))
LowvsHigh3months <- data.frame(topTags(lrt3monthHighvsLow,n=nrow(table),sort="none"))
FDR.table.LowvsHigh$month3lowvs3high_FDR<-LowvsHigh3months[,43]
LogFC.table.LowvsHigh$month3lowvs3high_logFC<-LowvsHigh3months[,39]

#3_5 months High vs Low
lrt3_5monthHighvsLow <- glmLRT(fit, contrast = c(0,0,0,0,0,0,1,0,0,1,0,0))
LowvsHigh3_5months <- data.frame(topTags(lrt3_5monthHighvsLow,n=nrow(table),sort="none"))
FDR.table.LowvsHigh$month3_5lowvs3_5high_FDR<-LowvsHigh3_5months[,43]
LogFC.table.LowvsHigh$month3_5lowvs3_5high_logFC<-LowvsHigh3_5months[,39]

#4 months High vs Low
lrt4monthHighvsLow <- glmLRT(fit, contrast = c(0,0,0,0,0,0,1,0,0,0,1,0))
LowvsHigh4months <- data.frame(topTags(lrt4monthHighvsLow,n=nrow(table),sort="none"))
FDR.table.LowvsHigh$month4lowvs4high_FDR<-LowvsHigh4months[,43]
LogFC.table.LowvsHigh$month4lowvs4high_logFC<-LowvsHigh4months[,39]

#4_5 months High vs Low
lrt4_5monthHighvsLow <- glmLRT(fit, contrast = c(0,0,0,0,0,0,1,0,0,0,0,1))
LowvsHigh4_5months <- data.frame(topTags(lrt4_5monthHighvsLow,n=nrow(table),sort="none"))
FDR.table.LowvsHigh$month4_5lowvs4_5high_FDR<-LowvsHigh4_5months[,43]
LogFC.table.LowvsHigh$month4_5lowvs4_5high_logFC<-LowvsHigh4_5months[,39]

######withinhigh####

#2 months to 2_5 high
lrt2_5monthvs2high <- glmLRT(fit, contrast = c(0,1,0,0,0,0,0,0,0,0,0,0))
High2_5monthvs2 <- data.frame(topTags(lrt2_5monthvs2high,n=nrow(table),sort="none"))
#back to 2M i.e. 2M->2_5M 2M->3M 2M->3_5M 2M->4M 2M->4_5M 
FDR.table.Highvs2<-High2_5monthvs2[,c(1,43)]
colnames(FDR.table.Highvs2)[colnames(FDR.table.Highvs2)=="FDR"] <- "month2_5vs2_high_FDR"
LogFC.table.Highvs2<-High2_5monthvs2[,c(1,39)]
colnames(LogFC.table.Highvs2)[colnames(LogFC.table.Highvs2)=="logFC"] <- "month2_5vs2_high_logFC"
#sequential i.e. 2M->2_5M->3M->3_5M->4M->4_5M
FDR.table.HighSeq<-High2_5monthvs2[,c(1,43)]
colnames(FDR.table.HighSeq)[colnames(FDR.table.HighSeq)=="FDR"] <- "month2_5vs2_high_FDR"
LogFC.table.HighSeq<-High2_5monthvs2[,c(1,39)]
colnames(LogFC.table.HighSeq)[colnames(LogFC.table.HighSeq)=="logFC"] <- "month2_5vs2_high_logFC"

#2 months to 3 high
lrt3monthvs2high <- glmLRT(fit, contrast = c(0,0,1,0,0,0,0,0,0,0,0,0))
High3monthvs2 <- data.frame(topTags(lrt3monthvs2high,n=nrow(table),sort="none"))
FDR.table.Highvs2$month3vs2_high_FDR<-High3monthvs2[,43]
LogFC.table.Highvs2$month3vs2_high_logFC<-High3monthvs2[,39]

#2 months to 3_5 high
lrt3_5monthvs2high <- glmLRT(fit, contrast = c(0,0,0,1,0,0,0,0,0,0,0,0))
High3_5monthvs2 <- data.frame(topTags(lrt3_5monthvs2high,n=nrow(table),sort="none"))
FDR.table.Highvs2$month3_5vs2_high_FDR<-High3_5monthvs2[,43]
LogFC.table.Highvs2$month3_5vs2_high_logFC<-High3_5monthvs2[,39]

#2 months to 4 high
lrt4monthvs2high <- glmLRT(fit, contrast = c(0,0,0,0,1,0,0,0,0,0,0,0))
High4monthvs2 <- data.frame(topTags(lrt4monthvs2high,n=nrow(table),sort="none"))
FDR.table.Highvs2$month4vs2_high_FDR<-High4monthvs2[,43]
LogFC.table.Highvs2$month4vs2_high_logFC<-High4monthvs2[,39]

#2 months to 4_5 high
lrt4_5monthvs2high <- glmLRT(fit, contrast = c(0,0,0,0,0,1,0,0,0,0,0,0))
High4_5monthvs2 <- data.frame(topTags(lrt4_5monthvs2high,n=nrow(table),sort="none"))
FDR.table.Highvs2$month4_5vs2_high_FDR<-High4_5monthvs2[,43]
LogFC.table.Highvs2$month4_5vs2_high_logFC<-High4_5monthvs2[,39]

#2_5 months to 3 high
lrt3monthvs2_5high <- glmLRT(fit, contrast = c(0,-1,1,0,0,0,0,0,0,0,0,0))
High3monthvs2_5 <- data.frame(topTags(lrt3monthvs2_5high,n=nrow(table),sort="none"))
FDR.table.HighSeq$month3vs2_5_high_FDR<-High3monthvs2_5[,43]
LogFC.table.HighSeq$month3vs2_5_high_logFC<-High3monthvs2_5[,39]

#2_5 months to 3_5 high
lrt3_5monthvs2_5high <- glmLRT(fit, contrast = c(0,-1,0,1,0,0,0,0,0,0,0,0))

#2_5 months to 4 high
lrt4monthvs2_5high <- glmLRT(fit, contrast = c(0,-1,0,0,1,0,0,0,0,0,0,0))

#2_5 months to 4_5 high
lrt4_5monthvs2_5high <- glmLRT(fit, contrast = c(0,-1,0,0,0,1,0,0,0,0,0,0))

#3 months to 3_5 high
lrt3_5monthvs3high <- glmLRT(fit, contrast = c(0,0,-1,1,0,0,0,0,0,0,0,0))
High3_5monthvs3 <- data.frame(topTags(lrt3_5monthvs3high,n=nrow(table),sort="none"))
FDR.table.HighSeq$month3_5vs3_high_FDR<-High3_5monthvs3[,43]
LogFC.table.HighSeq$month3_5vs3_high_logFC<-High3_5monthvs3[,39]

#3 months to 4 high
lrt4monthvs3high <- glmLRT(fit, contrast = c(0,0,-1,0,1,0,0,0,0,0,0,0))

#3 months to 4_5 high
lrt4_5monthvs3high <- glmLRT(fit, contrast = c(0,0,-1,0,0,1,0,0,0,0,0,0))

#3_5 months to 4 high
lrt4monthvs3_5high <- glmLRT(fit, contrast = c(0,0,0,-1,1,0,0,0,0,0,0,0))
High4monthvs3_5 <- data.frame(topTags(lrt4monthvs3_5high,n=nrow(table),sort="none"))
FDR.table.HighSeq$month4vs3_5_high_FDR<-High4monthvs3_5[,43]
LogFC.table.HighSeq$month4vs3_5_high_logFC<-High4monthvs3_5[,39]

#3_5 months to 4_5 high
lrt4_5monthvs3_5high <- glmLRT(fit, contrast = c(0,0,0,-1,0,1,0,0,0,0,0,0))

#4 months to 4_5 high
lrt4_5monthvs4high <- glmLRT(fit, contrast = c(0,0,0,0,-1,1,0,0,0,0,0,0))
High4_5monthvs4 <- data.frame(topTags(lrt4_5monthvs4high,n=nrow(table),sort="none"))
FDR.table.HighSeq$month4_5vs4_high_FDR<-High4_5monthvs4[,43]
LogFC.table.HighSeq$month4_5vs4_high_logFC<-High4_5monthvs4[,39]

#low back to 2M high#

#2M high to 2M low
lrt2Mlowvs2Mhigh <- glmLRT(fit, contrast = c(0,0,0,0,0,0,1,0,0,0,0,0))
Low2MvssHigh2M <- data.frame(topTags(lrt2Mlowvs2Mhigh,n=nrow(table),sort="none"))
#back to 2M high i.e. 2Mhigh->2_5Mlow 2Mhigh->3Mlow 2Mhigh->3_5Mlow 2Mhigh->4Mlow 2Mhigh->4_5Mlow 
FDR.table.Lowvs2MHigh<-Low2MvssHigh2M[,c(1,43)]
colnames(FDR.table.Lowvs2MHigh)[colnames(FDR.table.Lowvs2MHigh)=="FDR"] <- "month2lowvs2high_FDR"
LogFC.table.Lowvs2MHigh<-Low2MvssHigh2M[,c(1,39)]
colnames(LogFC.table.Lowvs2MHigh)[colnames(LogFC.table.Lowvs2MHigh)=="logFC"] <- "month2lowvs2high_logFC"

#2M high to 2_5M low 
lrt2_5Mlowvs2Mhigh <- glmLRT(fit, contrast = c(0,1,0,0,0,0,1,1,0,0,0,0))
Low2_5MvsHigh2M <- data.frame(topTags(lrt2_5Mlowvs2Mhigh,n=nrow(table),sort="none"))
FDR.table.Lowvs2MHigh$month2_5lowvs2high_FDR<-Low2_5MvsHigh2M[,43]
LogFC.table.Lowvs2MHigh$month2_5lowvs2high_logFC<-Low2_5MvsHigh2M[,39]

#2M high to 3M low
lrt3Mlowvs2Mhigh <- glmLRT(fit, contrast = c(0,0,1,0,0,0,1,0,1,0,0,0))
Low3MvsHigh2M <- data.frame(topTags(lrt3Mlowvs2Mhigh,n=nrow(table),sort="none"))
FDR.table.Lowvs2MHigh$month3lowvs2high_FDR<-Low3MvsHigh2M[,43]
LogFC.table.Lowvs2MHigh$month3lowvs2high_logFC<-Low3MvsHigh2M[,39]

#2M high to 3_5M low
lrt3_5Mlowvs2Mhigh <- glmLRT(fit, contrast = c(0,0,0,1,0,0,1,0,0,1,0,0))
Low3_5MvsHigh2M <- data.frame(topTags(lrt3_5Mlowvs2Mhigh,n=nrow(table),sort="none"))
FDR.table.Lowvs2MHigh$month3_5lowvs2high_FDR<-Low3_5MvsHigh2M[,43]
LogFC.table.Lowvs2MHigh$month3_5lowvs2high_logFC<-Low3_5MvsHigh2M[,39]

#2M high to 4M low
lrt4Mlowvs2Mhigh <- glmLRT(fit, contrast = c(0,0,0,0,1,0,1,0,0,0,1,0))
Low4MvsHigh2M <- data.frame(topTags(lrt4Mlowvs2Mhigh,n=nrow(table),sort="none"))
FDR.table.Lowvs2MHigh$month4lowvs2high_FDR<-Low4MvsHigh2M[,43]
LogFC.table.Lowvs2MHigh$month4lowvs2high_logFC<-Low4MvsHigh2M[,39]

#2M high to 4_5M low
lrt4_5Mlowvs2Mhigh <- glmLRT(fit, contrast = c(0,0,0,0,0,1,1,0,0,0,0,1))
Low4_5MvsHigh2M <- data.frame(topTags(lrt4_5Mlowvs2Mhigh,n=nrow(table),sort="none"))
FDR.table.Lowvs2MHigh$month4_5lowvs2high_FDR<-Low4_5MvsHigh2M[,43]
LogFC.table.Lowvs2MHigh$month4_5lowvs2high_logFC<-Low4_5MvsHigh2M[,39]

#######Alt time interaction#########
colnames(design)
AltTimeinteraction<-glmLRT(fit, coef=8:12)
AltTimeinteraction.tags <- data.frame(topTags(AltTimeinteraction,n=nrow(table),sort="none"))
AltTimeinteraction.tags.FDR05 <- subset(AltTimeinteraction.tags, FDR < 0.05)
length(AltTimeinteraction.tags.FDR05$Name) #43446

##########month################

Month<-glmLRT(fit, coef=2:6)
Month.tags <- data.frame(topTags(Month,n=nrow(table),sort="none"))
Month.tags.FDR05 <- subset(Month.tags, FDR < 0.05)
length(Month.tags.FDR05$Name) #19338

########alltitude###############

Alt<-glmLRT(fit, coef=7)
Alt.tags <- data.frame(topTags(Alt,n=nrow(table),sort="none"))
Alt.tags.FDR05 <- subset(Alt.tags, FDR < 0.05)
length(Alt.tags.FDR05$Name) #2370



#some figures

#low to 2m
FDR.table.low.from2M.sigrows<-FDR.table.Lowvs2[apply(FDR.table.Lowvs2[, -1], MARGIN = 1, function(x) any(x < 0.05)), ]
sig_geneid<-FDR.table.low.from2M.sigrows$Name
LogFC.table.low.from2M.sigFDR<- LogFC.table.Lowvs2[LogFC.table.Lowvs2$Name %in% sig_geneid, ]

df <- melt(LogFC.table.low.from2M.sigFDR, "Name")
ggplot(df, aes(variable, value,group = Name)) +
  geom_line()

LogFC.table.low.from2M.sigFDR$month2vs2_haw_logFC<-rep(0,nrow(LogFC.table.low.from2M.sigFDR))
LogFC.table.low.from2M.sigFDR<-LogFC.table.low.from2M.sigFDR[,c(1,7,2,3,4,5,6)]
df <- melt(LogFC.table.low.from2M.sigFDR, "Name")

ggplot(df, aes(variable, value,group = Name)) +
  geom_line() +
  ggtitle("Time Series 2M->3M 2M->4M 2M->5M 2M->6M") +
  scale_x_discrete("Time series from 2M Low",labels=c('2M','2_5M','3M','3_5M', '4M','4_5')) +
  scale_y_continuous("Log Fold Change")

#low consequtive

FDR.table.low.seq.sigrows<-FDR.table.LowSeq[apply(FDR.table.LowSeq[, -1], MARGIN = 1, function(x) any(x < 0.05)), ]
sig_geneid<-FDR.table.low.seq.sigrows$Name
LogFC.table.low.seq.sigFDR<- LogFC.table.LowSeq[LogFC.table.LowSeq$Name %in% sig_geneid, ]

df <- melt(LogFC.table.low.seq.sigFDR, "Name")
ggplot(df, aes(variable, value,group = Name)) +
  geom_line()


#high from 2M

FDR.table.high.from2M.sigrows<-FDR.table.Highvs2[apply(FDR.table.Highvs2[, -1], MARGIN = 1, function(x) any(x < 0.05)), ]
sig_geneid<-FDR.table.high.from2M.sigrows$Name
LogFC.table.high.from2M.sigFDR<- LogFC.table.Highvs2[LogFC.table.Highvs2$Name %in% sig_geneid, ]

df <- melt(LogFC.table.high.from2M.sigFDR, "Name")
ggplot(df, aes(variable, value,group = Name)) +
  geom_line()

LogFC.table.high.from2M.sigFDR$month2vs2_haw_logFC<-rep(0,nrow(LogFC.table.high.from2M.sigFDR))
LogFC.table.high.from2M.sigFDR<-LogFC.table.high.from2M.sigFDR[,c(1,7,2,3,4,5,6)]
df <- melt(LogFC.table.high.from2M.sigFDR, "Name")

ggplot(df, aes(variable, value,group = Name)) +
  geom_line() +
  ggtitle("Time Series 2M->3M 2M->4M 2M->5M 2M->6M") +
  scale_x_discrete("Time series from 2M High",labels=c('2M','2_5M','3M','3_5M', '4M','4_5')) +
  scale_y_continuous("Log Fold Change")


#sequential high

FDR.table.high.seq.sigrows<-FDR.table.HighSeq[apply(FDR.table.HighSeq[, -1], MARGIN = 1, function(x) any(x < 0.05)), ]
sig_geneid<-FDR.table.high.seq.sigrows$Name
LogFC.table.high.seq.sigFDR<- LogFC.table.HighSeq[LogFC.table.HighSeq$Name %in% sig_geneid, ]

df <- melt(LogFC.table.high.seq.sigFDR, "Name")
ggplot(df, aes(variable, value,group = Name)) +
  geom_line()

#from 2M high

FDR.table.low.from2Mhigh.sigrows<-FDR.table.Lowvs2MHigh[apply(FDR.table.Lowvs2MHigh[, -1], MARGIN = 1, function(x) any(x < 0.05)), ]
sig_geneid<-FDR.table.low.from2Mhigh.sigrows$Name
LogFC.table.Lowvs2MHigh.sigFDR<- LogFC.table.Lowvs2MHigh[LogFC.table.Lowvs2MHigh$Name %in% sig_geneid, ]

df <- melt(LogFC.table.Lowvs2MHigh.sigFDR, "Name")
ggplot(df, aes(variable, value,group = Name)) +
  geom_line()




#making tables

#Low
annotation.cerosi<-table[,c(1,4,8,22,36,33,34)]
colnames(table)
table.low<-merge(annotation.cerosi,LogFC.table.Lowvs2,by="Name")
FDR.table.Lowvs2$Lowest_FDR_TimeSeriesLowvs2<-apply(FDR.table.Lowvs2[,2:6],1,min)
table.low<-merge(table.low,FDR.table.Lowvs2,by="Name")
table.low<-merge(table.low,LogFC.table.LowSeq,by="Name")
FDR.table.LowSeq$Lowest_FDR_TimeSeriesLowSeq<-apply(FDR.table.LowSeq[,2:6],1,min)
table.low<-merge(table.low,FDR.table.LowSeq,by="Name")
table.low<-merge(table.low,LogFC.table.LowvsHigh,by="Name")
FDR.table.LowvsHigh$Lowest_FDR_SeriesLowvsHigh<-apply(FDR.table.LowvsHigh[,2:7],1,min)
table.low<-merge(table.low,FDR.table.LowvsHigh,by="Name")
table.low<-merge(table.low,LogFC.table.Lowvs2MHigh,by="Name")
FDR.table.Lowvs2MHigh$Lowest_FDR_SeriesLowvsHigh2M<-apply(FDR.table.Lowvs2MHigh[,2:7],1,min)
table.low<-merge(table.low,FDR.table.Lowvs2MHigh,by="Name")

head(LogFC.table.Lowvs2)
head(FDR.table.Lowvs2)
head(LogFC.table.LowSeq)
head(FDR.table.LowSeq)
head(LogFC.table.LowvsHigh)
head(FDR.table.LowvsHigh)
head(LogFC.table.Lowvs2MHigh)
head(FDR.table.Lowvs2MHigh)

Month.tags.FDR<-Month.tags[,c(1,47)]
Alt.tags.FDR<-Alt.tags[,c(1,43)]
AltTimeinteraction.tags.FDR<-AltTimeinteraction.tags[,c(1,47)]
colnames(Month.tags.FDR)[colnames(Month.tags.FDR)=="FDR" ]<- "Month_FDR"
colnames(Alt.tags.FDR)[colnames(Alt.tags.FDR)=="FDR" ]<- "Alt_FDR"
colnames(AltTimeinteraction.tags.FDR)[colnames(AltTimeinteraction.tags.FDR)=="FDR" ]<- "AltTimeInteraction_FDR"


table.low<-merge(table.low,Month.tags.FDR,by="Name")
table.low<-merge(table.low,Alt.tags.FDR,by="Name")
table.low<-merge(table.low,AltTimeinteraction.tags.FDR,by="Name")
colnames(table.low)
write.table(table.low,"BackToHigh2M/Table.cerosi.low",quote=FALSE,row.names=FALSE,sep="\t")
#High

table.high<-merge(annotation.cerosi,LogFC.table.Highvs2,by="Name")
FDR.table.Highvs2$Lowest_FDR_TimeSeriesHighvs2<-apply(FDR.table.Highvs2[,2:6],1,min)
table.high<-merge(table.high,FDR.table.Highvs2,by="Name")
table.high<-merge(table.high,LogFC.table.HighSeq,by="Name")
FDR.table.HighSeq$Lowest_FDR_TimeSeriesHighSeq<-apply(FDR.table.HighSeq[,2:6],1,min)
table.high<-merge(table.high,FDR.table.HighSeq,by="Name")
table.high<-merge(table.high,LogFC.table.LowvsHigh,by="Name")
table.high<-merge(table.high,FDR.table.LowvsHigh,by="Name")

head(LogFC.table.Highvs2)
head(FDR.table.Highvs2)
head(LogFC.table.HighSeq)
head(FDR.table.HighSeq)
head(LogFC.table.LowvsHigh)
head(FDR.table.LowvsHigh)

table.high<-merge(table.high,Month.tags.FDR,by="Name")
table.high<-merge(table.high,Alt.tags.FDR,by="Name")
table.high<-merge(table.high,AltTimeinteraction.tags.FDR,by="Name")
colnames(table.high)
write.table(table.high,"BackToHigh2M/Table.cerosi.high",quote=FALSE,row.names=FALSE,sep="\t")



########Heatmaps#########
setwd("~/Documents/Cerasi/Cerosi/salmonResults/")
library(RColorBrewer)
library(heatmap3)
library(heatmap.plus)
library(gplots)
cerosi.high<-read.table("Table.cerosi.high",sep="\t",header=TRUE,row.names=NULL,stringsAsFactors = FALSE)
cerosi.low<-read.table("Table.cerosi.low",sep="\t",header=TRUE,row.names=NULL,stringsAsFactors = FALSE)
colnames(cerosi.high)

colnames(cerosi.low)
test<-merge(cerosi.high,cerosi.low,by="Name")
colnames(test)

seq<-as.data.frame(test[c(6,8:12,65:69)])
seq<-seq[!with(seq,is.na(gene_id_pom_blast_tophit.x)), ]

#aggregate based on pom gene ID
seq<-aggregate(seq[-1], seq[1], mean)
seq[seq < -2] <- -2
seq[seq > 2] <- 2

seq.mat<-as.matrix(seq[-1]) 


#heatmap3(choppedData,col=brewer.pal(10,"PiYG"),Colv=F,labRow=F)

heatmap.2(seq.mat,col=brewer.pal(10,"PiYG"),trace="none",Colv=F,labRow=F,margins=c(13,13))

pdf("~/Documents/Cerasi/Cerosi/salmonResults/Cerasi.From2M.pdf")
heatmap.2(seq.mat,col=brewer.pal(10,"PiYG"),trace="none",Colv=F,labRow=F,margins=c(13,13))
dev.off()

tiff("~/Documents/Cerasi/Cerosi/salmonResults/Cerasi.From2M.tiff", height = 6, width = 6, units = 'in', res = 300)
heatmap.2(seq.mat,col=brewer.pal(10,"PiYG"),trace="none",Colv=F,labRow=F,margins=c(13,13))
dev.off()

#both from 2M low
seq.from2Mlow<-as.data.frame(test[c(6,43:48,65:69)])
seq.from2Mlow<-seq.from2Mlow[!with(seq.from2Mlow,is.na(gene_id_pom_blast_tophit.x)), ]
seq.from2Mlow<-aggregate(seq.from2Mlow[-1], seq.from2Mlow[1], mean)

seq.from2Mlow.mat<-as.matrix(seq.from2Mlow[-1])
seq.from2Mlow.mat[seq.from2Mlow.mat < -2] <- -2
seq.from2Mlow.mat[seq.from2Mlow.mat > 2] <- 2

tiff("~/Documents/Cerasi/Cerosi/salmonResults/Cerasi.From2Mlow.tiff", height = 6, width = 6, units = 'in', res = 300)
heatmap.2(seq.from2Mlow.mat,col=brewer.pal(10,"PiYG"),trace="none",Colv=F,labRow=F,margins=c(13,13))
dev.off()


#True sequential
seq.true<-as.data.frame(test[c(6,19:23,76:80)])
seq.true<-seq.true[!with(seq.true,is.na(gene_id_pom_blast_tophit.x)), ]
seq.true<-aggregate(seq.true[-1], seq.true[1], mean)

seq.true.mat<-as.matrix(seq.true[-1])
seq.true.mat[seq.true.mat < -2] <- -2
seq.true.mat[seq.true.mat > 2] <- 2

tiff("~/Documents/Cerasi/Cerosi/salmonResults/Cerasi.sequental.tiff", height = 6, width = 6, units = 'in', res = 300)
heatmap.2(seq.true.mat,col=brewer.pal(10,"PiYG"),trace="none",Colv=F,labRow=F,margins=c(13,13))
dev.off()

#end

























#######contrasts between in high and low######

#high
monthhighlowsep<-factor(c(rep("2_5M",4),rep("2M",4),rep("3_5M",4),rep("3M",4),rep("4_5M",4),rep("4M",4)))

table.sampleDatahigh<-data.frame(Sample=colnames(table.high.dge),monthhighlowsep)

#model
levels(monthhighlowsep)
monthHigh_<-factor(monthhighlowsep, levels = c("2M", "2_5M", "3M","3_5M","4M","4_5M"))
levels(monthHigh_)
designHigh <- model.matrix(~monthHigh_)
rownames(designHigh) <- colnames(table.high.dge)
designHigh

#> designHigh
#(Intercept) monthHigh_2_5M monthHigh_3M monthHigh_3_5M monthHigh_4M monthHigh_4_5M
#c_Hi_2_5M_01_exp_count           1              1            0              0            0              0
#c_Hi_2_5M_02_exp_count           1              1            0              0            0              0
#c_Hi_2_5M_04_exp_count           1              1            0              0            0              0
#c_Hi_2_5M_05_exp_count           1              1            0              0            0              0
#c_Hi_2M_01_exp_count             1              0            0              0            0              0
#c_Hi_2M_02_exp_count             1              0            0              0            0              0
#c_Hi_2M_03_exp_count             1              0            0              0            0              0
#c_Hi_2M_04_exp_count             1              0            0              0            0              0
#c_Hi_3_5M_01_exp_count           1              0            0              1            0              0
#c_Hi_3_5M_02_exp_count           1              0            0              1            0              0
#c_Hi_3_5M_03_exp_count           1              0            0              1            0              0
#c_Hi_3_5M_05_exp_count           1              0            0              1            0              0
#c_Hi_3M_02_exp_count             1              0            1              0            0              0
#c_Hi_3M_03_exp_count             1              0            1              0            0              0
#c_Hi_3M_04_exp_count             1              0            1              0            0              0
#c_Hi_3M_05_exp_count             1              0            1              0            0              0
#c_Hi_4_5M_02_exp_count           1              0            0              0            0              1
#c_Hi_4_5M_03_exp_count           1              0            0              0            0              1
#c_Hi_4_5M_04_exp_count           1              0            0              0            0              1
#c_Hi_4_5M_05_exp_count           1              0            0              0            0              1
#c_Hi_4M_01_exp_count             1              0            0              0            1              0
#c_Hi_4M_02_exp_count             1              0            0              0            1              0
#c_Hi_4M_03_exp_count             1              0            0              0            1              0
#c_Hi_4M_05_exp_count             1              0            0              0            1              0

#########estimate disperson squareroot=coefficient of variation of biological variation#######
table.high.dge<- estimateDisp(table.high.dge, designHigh, robust=TRUE)
table.high.dge$common.dispersion
#0.1347917
sqrt(0.1347917)
#0.367 coefficent of variation


#plot
plotBCV(table.high.dge)

#can look at trended dispersion that is used for QL (Quasi-Likelihood pipeline)
fit <- glmFit(table.high.dge, designHigh)
colnames(designHigh)

fit <- glmQLFit(table.high.dge, designHigh)
plotQLDisp(fit)


##############contrasts##################

fit <- glmFit(table.high.dge, designHigh)

colnames(designHigh)
#[1] "(Intercept)"    "monthHigh_2_5M" "monthHigh_3M"   "monthHigh_3_5M" "monthHigh_4M"   "monthHigh_4_5M"

#2_5 to two months
lrtHigh2_5monthvs2 <- glmLRT(fit,coef=2)

#3 to two months
lrtHigh3monthvs2 <- glmLRT(fit,coef=3)

#3_5 to two months
lrtHigh3_5monthvs2 <- glmLRT(fit,coef=4)

#4 to two months
lrtHigh4monthvs2 <- glmLRT(fit,coef=5)

#4_5 to two months
lrtHigh4_5monthvs2 <- glmLRT(fit,coef=6)

#3 to 2_5 months
lrtHigh3monthvs2_5 <- glmLRT(fit,contrast=c(0,-1,1,0,0,0))

#3_5 to 2_5 months
lrtHigh3_5monthvs2_5 <- glmLRT(fit,contrast=c(0,-1,0,1,0,0))

#4 to 2_5 months
lrtHigh4monthvs2_5 <- glmLRT(fit,contrast=c(0,-1,0,0,1,0))

#4_5 to 2_5 months
lrtHigh4_5monthvs2_5 <- glmLRT(fit,contrast=c(0,-1,0,0,0,1))

#3_5 to 3 months
lrtHigh3_5monthvs3 <- glmLRT(fit,contrast=c(0,0,-1,1,0,0))

#4 to 3 months
lrtHigh4monthvs3 <- glmLRT(fit,contrast=c(0,0,-1,0,1,0))

#4_5 to 3 months
lrtHigh4_5monthvs3 <- glmLRT(fit,contrast=c(0,0,-1,0,0,1))

#4 to 3_5 months
lrtHigh4monthvs3_5 <- glmLRT(fit,contrast=c(0,0,0,-1,1,0))

#4_5 to 3_5 months
lrtHigh4_5monthvs3_5 <- glmLRT(fit,contrast=c(0,0,0,-1,0,1))

#4_5 to 4 months
lrtHigh4_5monthvs4 <- glmLRT(fit,contrast=c(0,0,0,0,-1,1))

#between any months
lrtanydiffinHighMonths <- glmLRT(fit,coef=2:6)


#low

monthhighlowsep<-factor(c(rep("2_5M",4),rep("2M",4),rep("3_5M",4),rep("3M",4),rep("4_5M",4),rep("4M",4)))

table.sampleDatalow<-data.frame(Sample=colnames(table.low.dge),monthhighlowsep)

#model
levels(monthhighlowsep)
monthLow_<-factor(monthhighlowsep, levels = c("2M", "2_5M", "3M","3_5M","4M","4_5M"))
levels(monthLow_)
designLow <- model.matrix(~monthLow_)
rownames(designLow) <- colnames(table.low.dge)
designLow

#> designLow
#(Intercept) monthLow_2_5M monthLow_3M monthLow_3_5M monthLow_4M monthLow_4_5M
#c_Lo_2_5M_01_exp_count           1             1           0             0           0             0
#c_Lo_2_5M_02_exp_count           1             1           0             0           0             0
#c_Lo_2_5M_03_exp_count           1             1           0             0           0             0
#c_Lo_2_5M_04_exp_count           1             1           0             0           0             0
#c_Lo_2M_01_exp_count             1             0           0             0           0             0
#c_Lo_2M_02_exp_count             1             0           0             0           0             0
#c_Lo_2M_03_exp_count             1             0           0             0           0             0
#c_Lo_2M_04_exp_count             1             0           0             0           0             0
#c_Lo_3_5M_01_exp_count           1             0           0             1           0             0
#c_Lo_3_5M_02_exp_count           1             0           0             1           0             0
#c_Lo_3_5M_03_exp_count           1             0           0             1           0             0
#c_Lo_3_5M_04_exp_count           1             0           0             1           0             0
#c_Lo_3M_01_exp_count             1             0           1             0           0             0
#c_Lo_3M_02_exp_count             1             0           1             0           0             0
#c_Lo_3M_03_exp_count             1             0           1             0           0             0
#c_Lo_3M_04_exp_count             1             0           1             0           0             0
#c_Lo_4_5M_01_exp_count           1             0           0             0           0             1
#c_Lo_4_5M_02_exp_count           1             0           0             0           0             1
#c_Lo_4_5M_03_exp_count           1             0           0             0           0             1
#c_Lo_4_5M_04_exp_count           1             0           0             0           0             1
#c_Lo_4M_01_exp_count             1             0           0             0           1             0
#c_Lo_4M_02_exp_count             1             0           0             0           1             0
#c_Lo_4M_03_exp_count             1             0           0             0           1             0
#c_Lo_4M_04_exp_count             1             0           0             0           1             0

#########estimate disperson squareroot=coefficient of variation of biological variation#######
table.low.dge<- estimateDisp(table.low.dge, designLow, robust=TRUE)
table.low.dge$common.dispersion
#0.1744747
sqrt(0.1744747)
#0.4177 coefficent of variation


#plot
plotBCV(table.low.dge)

#can look at trended dispersion that is used for QL (Quasi-Likelihood pipeline)
fit <- glmFit(table.low.dge, designLow)
colnames(designLow)

fit <- glmQLFit(table.low.dge, designLow)
plotQLDisp(fit)


##############contrasts##################

fit <- glmFit(table.low.dge, designLow)

colnames(designLow)
#[1] "(Intercept)"   "monthLow_2_5M" "monthLow_3M"   "monthLow_3_5M" "monthLow_4M"   "monthLow_4_5M"

#2_5 to two months
lrtLow2_5monthvs2 <- glmLRT(fit,coef=2)

#3 to two months
lrtLow3monthvs2 <- glmLRT(fit,coef=3)

#3_5 to two months
lrtLow3_5monthvs2 <- glmLRT(fit,coef=4)

#4 to two months
lrtLow4monthvs2 <- glmLRT(fit,coef=5)

#4_5 to two months
lrtLow4_5monthvs2 <- glmLRT(fit,coef=6)

#3 to 2_5 months
lrtLow3monthvs2_5 <- glmLRT(fit,contrast=c(0,-1,1,0,0,0))

#3_5 to 2_5 months
lrtLow3_5monthvs2_5 <- glmLRT(fit,contrast=c(0,-1,0,1,0,0))

#4 to 2_5 months
lrtLow4monthvs2_5 <- glmLRT(fit,contrast=c(0,-1,0,0,1,0))

#4_5 to 2_5 months
lrtLow4_5monthvs2_5 <- glmLRT(fit,contrast=c(0,-1,0,0,0,1))

#3_5 to 3 months
lrtLow3_5monthvs3 <- glmLRT(fit,contrast=c(0,0,-1,1,0,0))

#4 to 3 months
lrtLow4monthvs3 <- glmLRT(fit,contrast=c(0,0,-1,0,1,0))

#4_5 to 3 months
lrtLow4_5monthvs3 <- glmLRT(fit,contrast=c(0,0,-1,0,0,1))

#4 to 3_5 months
lrtLow4monthvs3_5 <- glmLRT(fit,contrast=c(0,0,0,-1,1,0))

#4_5 to 3_5 months
lrtLow4_5monthvs3_5 <- glmLRT(fit,contrast=c(0,0,0,-1,0,1))

#4_5 to 4 months
lrtLow4_5monthvs4 <- glmLRT(fit,contrast=c(0,0,0,0,-1,1))

#between any months
lrtanydiffinLowMonths <- glmLRT(fit,coef=2:6)


