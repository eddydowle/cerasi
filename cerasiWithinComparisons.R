#this code is based around the WGCNA pipeline

setwd("~/Documents/Cerasi/Cerosi/salmonResults/BackToHigh2Mremove3M/")
setwd("~/Documents/Cerasi/Cerosi/salmonResults2018MainTranscript/BackToHigh2Mremove3Mand25sample/")
High<-read.table('Table.cerosi.high.fixannot',header=T,row.names=NULL,sep='\t')
low<-read.table('Table.cerosi.low.fixannot',header=T,row.names=NULL,sep='\t')

High<-read.table('Table.cerosi.high',header=T,row.names=NULL,sep='\t')
low<-read.table('Table.cerosi.low',header=T,row.names=NULL,sep='\t')
colnames(High)
colnames(Low)

library(gplots)
library(ggplot2)
library(reshape)
library(heatmap3)
library(RColorBrewer)
library(WGCNA)

#pom.apple.from2Mhaw<-apple[c(1,12:16)]
#pom.haw.from2Mhaw<-haw[c(1,2:5)] 

#various subsets
#pom.apple.from2Mhaw.6M<-apple[c(1,12:16,21,44,2)]
#pom.haw.from2Mhaw.6M<-haw[c(1,2:5,9,32)]


#test<-subset(pom.apple.from2Mhaw.6M , Apple6monthvs2Haw_FDR < 0.05 & ApplevsHaw6Month_FDR < 0.05)

#library(reshape)
#test2<-test[1:6]
#df <- melt(test2, "gene_id")
#ggplot(df, aes(variable, value,group = gene_id)) +
  #geom_line()

#test$FlyBase_FBgn_tophit.x.x






colnames(High)
pom.high.from2Mhigh<-High[c(1,8:11,16)]
colnames(low)
pom.low.from2Mhigh<-low[c(1,37:41,47)] 

head(pom.high.from2Mhigh)
head(pom.low.from2Mhigh)
high.low.from2Mhigh<-merge(pom.high.from2Mhigh,pom.low.from2Mhigh,by="Name")
high.low.from2Mhigh.sigAny<-subset(high.low.from2Mhigh,Lowest_FDR_TimeSeriesHighvs2 < 0.05 | Lowest_FDR_SeriesLowvsHigh2M < 0.05)
high.low.from2Mhigh.sigBoth<-subset(high.low.from2Mhigh,Lowest_FDR_TimeSeriesHighvs2< 0.05 & Lowest_FDR_SeriesLowvsHigh2M < 0.05)

head(high.low.from2Mhigh.sigAny)
head(high.low.from2Mhigh.sigBoth)
high.low.from2Mhigh.sigAny$Lowest_FDR_TimeSeriesHighvs2<-NULL
high.low.from2Mhigh.sigAny$Lowest_FDR_SeriesLowvsHigh2M<-NULL
high.low.from2Mhigh.sigBoth$Lowest_FDR_TimeSeriesHighvs2<-NULL
high.low.from2Mhigh.sigBoth$Lowest_FDR_SeriesLowvsHigh2M<-NULL

high.low.from2Mhigh.sigBoth.mat<-as.matrix(high.low.from2Mhigh.sigBoth[-1])
high.low.from2Mhigh.sigBoth.mat[high.low.from2Mhigh.sigBoth.mat < -2] <- -2
high.low.from2Mhigh.sigBoth.mat[high.low.from2Mhigh.sigBoth.mat > 2] <- 2

heatmap.2(high.low.from2Mhigh.sigBoth.mat,col=brewer.pal(10,"PiYG"),trace="none",Colv=F,labRow=F,margins=c(13,13),cexCol=0.8)
high.low.from2Mhigh.sigBoth.mat<-as.matrix(high.low.from2Mhigh.sigBoth[-1])
heatmap.2(high.low.from2Mhigh.sigBoth.mat,col=brewer.pal(10,"PiYG"),trace="none",Colv=F,labRow=F,margins=c(13,13))

allowWGCNAThreads(nThreads = 4)
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(t(high.low.from2Mhigh.sigBoth.mat), powerVector = powers, verbose = 5)

sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")
#20? in this case (softpower)

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

softPower = 20; #change from above
adjacency = adjacency(t(high.low.from2Mhigh.sigBoth.mat), power = softPower);

TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM
geneTree = hclust(as.dist(dissTOM), method = "average");

sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);

minModuleSize = 50;
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)


#pdf('~/Documents/Cerasi/ComparisonsBetweenSpecies/AllSamplesInEdgeRNoHighRemovedPom/wgcnaClusterDrosPomCerDec2016.pdf')
#tiff("~/Documents/Cerasi/ComparisonsBetweenSpecies/AllSamplesInEdgeRNoHighRemovedPom/wgcnaClusterPomCerDec2016.tiff")

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
#sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
#calculate eigengenes and then cluster modules eigengenes
MEList = moduleEigengenes(t(high.low.from2Mhigh.sigBoth.mat), colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average");
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

MEDissThres = 0.25
abline(h=MEDissThres, col = "red")
merge = mergeCloseModules(t(high.low.from2Mhigh.sigBoth.mat), dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors;
mergedMEs = merge$newMEs;
mergedMEs #to see newMEs merged

out<-data.frame(gene_id=high.low.from2Mhigh.sigBoth$Name,module=mergedColors)
write.table(out,'WgcnaModules.pom.2Mhaw.significantboth',sep="\t",row.names=F,quote=F)

sizeGrWindow(12, 9)
#tiff(file = "~/Documents/Cerasi/ComparisonsBetweenSpecies/AllSamplesInEdgeRNoHighRemovedPom/wgcnaClusterPomCerDec2016.eigengenes.tiff")
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#dev.off()


#to plot what each of the eigengenes do across the time series
xMat<-matrix(rep(1:9,1),nrow=9,ncol=ncol(mergedMEs),byrow=F)
xMat<-matrix(rep(1:12,1),nrow=9,ncol=ncol(mergedMEs),byrow=F)

matplot(xMat,as.matrix(mergedMEs),type="l")

#n.
#xMat<-matrix(rep(1:19,3),nrow=18,ncol=3,byrow=F)
#xMat<-matrix(rep(1:19,3),nrow=18,ncol=3,byrow=F)

#matplot(xMat,as.matrix(mergedMEs[,12:14]),type="l")
#matplot(xMat,as.matrix(mergedMEs[,1:12]),type="l")

months<-c(2,2.5,3.5,4,4.5)
test<-as.data.frame(cbind(c(2.5,3.5,4,4.5,2,2.5,3.5,4,4.5),as.numeric(mergedMEs$MEdarkgrey),c(rep("High",4),rep("Low",5))))
#dev.off()
test$V2=as.numeric(levels(test$V2))[test$V2]

ggplot(test, aes(V1, V2,group = V3,colour =V3)) +
  geom_line() +
  scale_y_continuous("LogFC",limits = c(-2, 2)) +
  scale_x_discrete("Month") +
  scale_color_manual(values=c("green4", "purple4"),guide = guide_legend(title = "Population")) +
  theme_bw()


#which are largest ME's
head(out)
library(plyr)
count (out$module)

#> count (out$module)
#x freq
#1         black  308 #both up
#2        brown4  997 #low down then up high down then flat
#3   darkorange2   91 #low flat high down up flat
#4 darkslateblue  617 #low up then down high down then flat
#5         green 6849 #low up then down high down then flat
#6          grey   99 #low down up down high down up
#7         ivory  239 #up then low down and high up
#8 paleturquoise  210 #down then low up and high down
#9     steelblue  765 #flat then low up and high down


#do it with one that has flybase IDS in it
high.low.from2Mhigh.sigBoth.ME<-merge(out,high.low.from2Mhigh.sigBoth,by="gene_id")

head(high.low.from2Mhigh.sigAny.ME)

flybase.geneid<-apple[,1:2]

high.low.from2Mhigh.sigBoth.ME<-merge(flybase.geneid,high.low.from2Mhigh.sigBoth.ME,by="gene_id")
head(high.low.from2Mhigh.sigBoth.ME)

write.table(high.low.from2Mhigh.sigBoth.ME,'high.low.from2Mhigh.sigBoth.ME',sep="\t",row.names=F,quote=F)











#########################
#significant either  pom#
#########################

high.low.from2Mhigh.sigAny.mat<-as.matrix(high.low.from2Mhigh.sigAny[-1])
high.low.from2Mhigh.sigAny.mat[high.low.from2Mhigh.sigAny.mat < -2] <- -2
high.low.from2Mhigh.sigAny.mat[high.low.from2Mhigh.sigAny.mat > 2] <- 2

heatmap.2(high.low.from2Mhigh.sigAny.mat,col=brewer.pal(10,"PiYG"),trace="none",Colv=F,labRow=F,margins=c(13,13),cexCol = 0.8)
high.low.from2Mhigh.sigAny.mat<-as.matrix(high.low.from2Mhigh.sigAny[-1])
heatmap.2(high.low.from2Mhigh.sigAny.mat,col=brewer.pal(10,"PiYG"),trace="none",Colv=F,labRow=F,margins=c(13,13))

allowWGCNAThreads(nThreads = 4)
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(t(high.low.from2Mhigh.sigAny.mat), powerVector = powers, verbose = 5)

sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")
#20? in this case (softpower)

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

softPower = 20; #change from above
adjacency = adjacency(t(high.low.from2Mhigh.sigAny.mat), power = softPower);

TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM
geneTree = hclust(as.dist(dissTOM), method = "average");

sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);

minModuleSize = 50;
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)


#pdf('~/Documents/Cerasi/ComparisonsBetweenSpecies/AllSamplesInEdgeRNoHighRemovedPom/wgcnaClusterDrosPomCerDec2016.pdf')
#tiff("~/Documents/Cerasi/ComparisonsBetweenSpecies/AllSamplesInEdgeRNoHighRemovedPom/wgcnaClusterPomCerDec2016.tiff")

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
#sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
#calculate eigengenes and then cluster modules eigengenes
MEList = moduleEigengenes(t(high.low.from2Mhigh.sigAny.mat), colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average");
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

MEDissThres = 0.25
abline(h=MEDissThres, col = "red")
merge = mergeCloseModules(t(high.low.from2Mhigh.sigAny.mat), dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors;
mergedMEs = merge$newMEs;
mergedMEs #to see newMEs merged

out<-data.frame(gene_id=high.low.from2Mhigh.sigAny$gene_id,module=mergedColors)
write.table(out,'WgcnaModules.pom.2Mhaw.significantEither',sep="\t",row.names=F,quote=F)

sizeGrWindow(12, 9)
#tiff(file = "~/Documents/Cerasi/ComparisonsBetweenSpecies/AllSamplesInEdgeRNoHighRemovedPom/wgcnaClusterPomCerDec2016.eigengenes.tiff")
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#dev.off()


#to plot what each of the eigengenes do across the time series
xMat<-matrix(rep(1:9,1),nrow=9,ncol=ncol(mergedMEs),byrow=F)
matplot(xMat,as.matrix(mergedMEs),type="l")

#n.

test<-as.data.frame(cbind(c(2,3,4,5,6,3,4,5,6),as.numeric(mergedMEs$MEbrown4),c(rep("Apple",5),rep("Haw",4))))
#dev.off()
test$V2=as.numeric(levels(test$V2))[test$V2]

ggplot(test, aes(V1, V2,group = V3,colour =V3)) +
  geom_line() +
  scale_y_continuous("LogFC",limits = c(-2, 2)) +
  scale_x_discrete("Month") +
  scale_color_manual(values=c("blue", "red"),guide = guide_legend(title = "Population")) +
  theme_bw()

#which are largest ME's
head(out)
library(plyr)
count (out$module)

#> count (out$module)
#x freq
#1           black  479 #both step  down in similar pattern
#2          brown4   69 #haw big drop at 4M else flat
#3            cyan  245 #flatish apple goes down at 6m
#4        darkgrey  284 #?weird
#5  darkolivegreen  877 #weird kinda go down except for 4m in both and 6m in apple
#6         darkred  128 #could be out of step with apple going first~
#7   darkturquoise  516 #flatish then haw goes down before apple 
#8     greenyellow  639 #both go steadily down
#9            grey   74 #apple drop at 4M else flat
#10     lightgreen 2246 #both go steadily up
#11        magenta  155 #flat then up at end in both
#12         purple  152 #apple flat haw down then up
#13            red  172 #apple up then down at end haw different
#14            tan  513 #both go up but apple goes slightly down at 6m


#do it with one that has flybase IDS in it
high.low.from2Mhigh.sigAny.ME<-merge(out,high.low.from2Mhigh.sigAny,by="gene_id")

head(high.low.from2Mhigh.sigAny.ME)

flybase.geneid<-apple[,1:2]

high.low.from2Mhigh.sigAny.ME<-merge(flybase.geneid,high.low.from2Mhigh.sigAny.ME,by="gene_id")
head(high.low.from2Mhigh.sigAny.ME)

write.table(high.low.from2Mhigh.sigAny.ME,'high.low.from2Mhigh.sigAny.ME',sep="\t",row.names=F,quote=F)













#just all in no significance
high.low.from2Mhigh<-merge(pom.apple.from2Mhaw,pom.haw.from2Mhaw,by="gene_id")

high.low.from2Mhigh.mat<-as.matrix(high.low.from2Mhigh[-1])

high.low.from2Mhigh.mat[high.low.from2Mhigh.mat < -2] <- -2
high.low.from2Mhigh.mat[high.low.from2Mhigh.mat > 2] <- 2

heatmap.2(high.low.from2Mhigh.mat,col=brewer.pal(10,"PiYG"),trace="none",Colv=F,labRow=F,margins=c(13,13))

high.low.from2Mhigh.mat<-as.matrix(high.low.from2Mhigh[-1])
heatmap.2(high.low.from2Mhigh.mat,col=brewer.pal(10,"PiYG"),trace="none",Colv=F,labRow=F,margins=c(13,13))

allowWGCNAThreads(nThreads = 4)
#soft threshold powers (defaults on manual)
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(t(high.low.from2Mhigh.mat), powerVector = powers, verbose = 5)

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
#12 in this case

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#calculate (correlation or distance) network adjacency
softPower = 12;
adjacency = adjacency(t(pomcer.flybase.mat), power = softPower);

#?adjacency

# Turn adjacency into topological overlap matrix (TOM)
#this is to minimise noise and spurious association calculate the topolocial overlap matrix and calculate dissimilarity
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM


#cluster using TOM
# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
#each leaf corresponds to a gene module identification identifies clades of similarity 

#if you want smaller modules (clades) you can set this smaller (visa versa)
minModuleSize = 50;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
#here this returns 24 modules (clades), 0 is reseved for unassigned genes then largest to smallest in module number
table(dynamicMods)

#then we can plot module assignment

#pdf('~/Documents/Cerasi/ComparisonsBetweenSpecies/AllSamplesInEdgeRNoHighRemovedPom/wgcnaClusterDrosPomCerDec2016.pdf')
#tiff("~/Documents/Cerasi/ComparisonsBetweenSpecies/AllSamplesInEdgeRNoHighRemovedPom/wgcnaClusterPomCerDec2016.tiff")

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
#sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
#dev.off()

#now we can caculate eigengenes (genes that expression pattern is similar which can be merged together)
#quantifies expression of whole module (clade) 
#calculate eigengenes and then cluster modules eigengenes
MEList = moduleEigengenes(t(high.low.from2Mhigh.mat), colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

#We choose a height cut of 0.25, corresponding to correlation of 0.75, to merge eigengenes
MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(t(high.low.from2Mhigh.mat), dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;
mergedMEs #to see newMEs merged

#write out new merged MEs

out<-data.frame(gene_id=high.low.from2Mhigh$gene_id,module=mergedColors)
write.table(out,'WgcnaModules.pom.2Mhaw.allgenes',sep="\t",row.names=F,quote=F)

#to see how that grouped our initial clusters
sizeGrWindow(12, 9)
#tiff(file = "~/Documents/Cerasi/ComparisonsBetweenSpecies/AllSamplesInEdgeRNoHighRemovedPom/wgcnaClusterPomCerDec2016.eigengenes.tiff")
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#dev.off()


#to plot what each of the eigengenes do across the time series
xMat<-matrix(rep(1:18,1),nrow=18,ncol=ncol(mergedMEs),byrow=F)
matplot(xMat,as.matrix(mergedMEs),type="l")

#n.
xMat<-matrix(rep(1:19,3),nrow=18,ncol=3,byrow=F)
matplot(xMat,as.matrix(mergedMEs[,12:14]),type="l")








############################################################################################################################################
#########################################################gregs method of clustering#########################################################
############################################################################################################################################

chop<-function(x,thresh=2) {
  elchop<-function(y) {
    if (is.na(y)==FALSE & abs(y) > thresh) {
      if (y > thresh) {y<-thresh}
      if (y < -thresh) {y<--thresh}
    }
    return(y)
  }
  sapply(x,elchop)
}

setwd('~/Documents/Cerasi/Cerosi/salmonResults2018MainTranscript/BackToHigh2Mremove3Mand25sample/')
High<-read.table('Table.cerosi.high.fixannot',header=T,row.names=NULL,sep='\t')
low<-read.table('Table.cerosi.low.fixannot',header=T,row.names=NULL,sep='\t')



colnames(High)
pom.high.from2Mhigh<-High[c(1,8:11)]
colnames(low)
pom.low.from2Mhigh<-low[c(1,37:41)] 

#appledat<-read.table('Table.apple.JuneEJD',header=T,row.names=NULL)
#hawdat<-read.table('Table.haw.JuneEJD',header=T,row.names=NULL)
#pomAll<-merge(appledat[,c(1,12:22)],hawdat[,1:10])
high.low.from2Mhigh<-merge(pom.high.from2Mhigh,pom.low.from2Mhigh,by="Name")
#high.low.from2Mhigh.sigAny<-subset(high.low.from2Mhigh,Lowest_FDR_TimeSeriesHighvs2 < 0.05 | Lowest_FDR_SeriesLowvsHigh2M < 0.05)
#high.low.from2Mhigh.sigBoth<-subset(high.low.from2Mhigh,Lowest_FDR_TimeSeriesHighvs2< 0.05 & Lowest_FDR_SeriesLowvsHigh2M < 0.05)



library(WGCNA)
library(RColorBrewer)
library(gplots)

#select only genes DE between host races within months 
#ind<-appledat$Lowest_FDR_FruitBetweenMonths < 0.05
#ind<-appledat$ApplevsHaw6Month_FDR < 0.05

ind<-low$Lowest_FDR_SeriesLowvsHigh2M < 0.05

#mat<-as.matrix(pomAll[ind,c(13:16,2:6)]) #here he switches haw and apple so haw is now first
#rownames(mat)<-pomAll$gene_id[ind]
mat<-as.matrix(high.low.from2Mhigh[ind,c(2:10)])
rownames(mat)<-high.low.from2Mhigh$Name[ind]

clus<-heatmap.2(mat,col=brewer.pal(10,"PiYG"),trace="none",Colv=F)
choppedData<-apply(mat,2,chop,thresh=2)
heatmap.2(choppedData,col=brewer.pal(10,"PiYG"),trace="none",Colv=F,Rowv=clus$rowDendrogram,cexCol=0.2)
#heatmap.2(high.low.from2Mhigh.sigBoth.mat,col=brewer.pal(10,"PiYG"),trace="none",Colv=F,labRow=F,margins=c(13,13),cexCol=0.8)

minModuleSize = 30;
#minModuleSize = 10 # modsize for just genes DE at 6 months
dendro<-as.hclust(clus$rowDendrogram)
plot(dendro)
a<-dist(mat, method = "euclidean")
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = dendro, distM = as.matrix(a),
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);

modNums<-length(unique(dynamicMods))

dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
#sizeGrWindow(8,6)
plotDendroAndColors(dendro, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

#find average patterns within each cluster and plot
eigenGenes<-matrix(nrow=ncol(mat),ncol=modNums)
for (j in unique(dynamicMods)) {
  for (i in 1:ncol(mat)) {
    eigenGenes[i,j]<-mean(mat[dynamicMods==j,i])
  }
}
test<-eigenGenes[,1]
test
plot(c(2,2.5,3.5,4,4.5),ylim=c(-2,2),c(0,test[1:4]),col='green4',type='l',xlab='Month',ylab='LFC',las=1,main=paste('n =',sep=' ')) #haw
lines(c(2,2.5,3.5,4,4.5),test[5:9],col='purple') #apple
lines(c(2,2.5,3.5,4,4.5),rep(0,5),lty=2)
legend("topleft",c("High","Low"),lty=c(1,1),col=c("green4","purple"))


plotTraj<-function(x,main,n) {
  range<-c(min(x,na.rm=T),max(x,na.rm=T))
  if (range[1] > 0) {range[1]<-0}
  if (range[2] < 0) {range[2]<-0}
  plot(c(2,2.5,3.5,4,4.5),c(0,x[1:4]),ylim=range,col='green4',type='l',xlab='Month',ylab='LFC',las=1,main=paste(main,'n =',n,sep=' ')) #haw
  lines(c(2,2.5,3.5,4,4.5),x[5:9],col='purple') #apple
  lines(c(2,2.5,3.5,4,4.5),rep(0,5),lty=2)
  legend("topleft",c("High","Low"),lty=c(1,1),col=c("green4","purple"))
}
plotTraj(eigenGenes[,1],'1',as.vector(table(dynamicMods))[1])
#plot all
pdf('CerClustersDE.WithinMonths.gregmethod.pdf')
par(mfcol=c(3,3))
for (i in 1:ncol(eigenGenes)) {
  plotTraj(eigenGenes[,i],i,as.vector(table(dynamicMods))[i])
}
dev.off()

##flyIdmap<-data.frame(gene_id=appledat$gene_id,flyid=appledat$FlyBase_FBgn_tophit.x.x,stringsAsFactors = F)
#flyIdmap$flyid[is.na(flyIdmap$flyid)]<-'noHit'
#flyIdmap["flyid"][is.na(flyIdmap["flyid"])] <- 'noHit'
#out<-merge(flyIdmap,data.frame(gene_id=rownames(mat),mat))
#out<-data.frame(out,module=dynamicMods)
out<-data.frame(data.frame(gene_id=rownames(mat),mat),module=dynamicMods)
write.csv(out,'CerModsDEwithinMonth.csv') # for all genes DE in at least one month

#counts within each month
colnames(appledat)
colnames(hawdat)

head(appledat$ApplevsHaw2Month_FDR)
length(appledat[appledat$ApplevsHaw2Month_FDR < 0.05,]) #45
head(appledat[appledat$ApplevsHaw2Month_FDR < 0.05,])
test<-appledat[appledat$ApplevsHaw2Month_FDR < 0.05,]
nrow(test) #1134
test<-appledat[appledat$ApplevsHaw3Month_FDR < 0.05,]
nrow(test) #45
test<-appledat[appledat$ApplevsHaw4Month_FDR < 0.05,]
nrow(test) #6
test<-appledat[appledat$ApplevsHaw5Month_FDR < 0.05,]
nrow(test) #1755
test<-appledat[appledat$ApplevsHaw6Month_FDR < 0.05,]
nrow(test) #329



################################################## DE across (any) time points #######################
colnames(High)
pom.high.from2Mhigh<-High[c(1,8:11,16)]
colnames(low)
pom.low.from2Mhigh<-low[c(1,37:41,47)] 

#appledat<-read.table('Table.apple.JuneEJD',header=T,row.names=NULL)
#hawdat<-read.table('Table.haw.JuneEJD',header=T,row.names=NULL)
#pomAll<-merge(appledat[,c(1,12:22)],hawdat[,1:10])
high.low.from2Mhigh<-merge(pom.high.from2Mhigh,pom.low.from2Mhigh,by="Name")


ind<-high.low.from2Mhigh$Lowest_FDR_TimeSeriesHighvs2 < 0.05 & high.low.from2Mhigh$Lowest_FDR_SeriesLowvsHigh2M < 0.05
#ind<-high.low.from2Mhigh$Lowest_FDR_from.TimeSeries.haw < 0.05 & high.low.from2Mhigh$Lowest_FDR_from2Mhaw.TimeSeries < 0.05

colnames(high.low.from2Mhigh)
#mat<-as.matrix(pomAll[ind,c(13:16,2:6)])
mat<-as.matrix(high.low.from2Mhigh[ind,c(2:5,7:11)])

rownames(mat)<-high.low.from2Mhigh$Name[ind]

clus<-heatmap.2(mat,col=brewer.pal(10,"PiYG"),trace="none",Colv=F)
choppedData<-apply(mat,2,chop,thresh=2)
heatmap.2(choppedData,col=brewer.pal(10,"PiYG"),trace="none",Colv=F,Rowv=clus$rowDendrogram,cexCol = 0.3)


minModuleSize = 30;
#minModuleSize = 10 # modsize for just genes DE at 6 months
dendro<-as.hclust(clus$rowDendrogram)
a<-dist(mat, method = "euclidean")
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = dendro, distM = as.matrix(a),
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);

modNums<-length(unique(dynamicMods))

#find average patterns within each cluster and plot
eigenGenes<-matrix(nrow=ncol(mat),ncol=modNums)
for (j in unique(dynamicMods)) {
  for (i in 1:ncol(mat)) {
    eigenGenes[i,j]<-mean(mat[dynamicMods==j,i])
  }
}
pdf('CerClustersDE.AcrossMonths.pdf')
par(mfcol=c(3,3))
for (i in 1:ncol(eigenGenes)) {
  plotTraj(eigenGenes[,i],i,as.vector(table(dynamicMods))[i])
}
dev.off()

#out<-merge(flyIdmap,data.frame(gene_id=rownames(mat),mat))
#out<-data.frame(out,module=dynamicMods)
out<-data.frame(data.frame(Name=rownames(mat),mat),module=dynamicMods)
write.csv(out,'CerModsDEacrossMonth.csv') # for all genes DE across at least one month

#counts
test<-appledat[appledat$Apple2monthvs2Haw_FDR < 0.05,]
nrow(test) #1134
test<-appledat[appledat$Apple3monthvs2Haw_FDR < 0.05,]
nrow(test) #41
test<-appledat[appledat$Apple4monthvs2Haw_FDR < 0.05,]
nrow(test) #1088
test<-appledat[appledat$Apple5monthvs2Haw_FDR < 0.05,]
nrow(test) #1836
test<-appledat[appledat$Apple6monthvs2Haw_FDR < 0.05,]
nrow(test) #3152

test<-hawdat[hawdat$X3monthvs2_haw_FDR.x < 0.05,]
nrow(test) #106
test<-hawdat[hawdat$month4vs2_haw_FDR < 0.05,]
nrow(test) #556
test<-hawdat[hawdat$month5vs2_haw_FDR < 0.05,]
nrow(test) #3125
test<-hawdat[hawdat$month6vs2_haw_FDR < 0.05,]
nrow(test) #3732

#linear

test<-appledat[appledat$X3monthvs2_apple_FDR.x < 0.05,]
nrow(test) #555
test<-appledat[appledat$month4vs3_apple_FDR < 0.05,]
nrow(test) #330
test<-appledat[appledat$month5vs4_apple_FDR < 0.05,]
nrow(test) #452
test<-appledat[appledat$month6vs5_apple_FDR < 0.05,]
nrow(test) #1207

test<-hawdat[hawdat$X3monthvs2_haw_FDR.x < 0.05,]
nrow(test) #106
test<-hawdat[hawdat$month4vs3_haw_FDR < 0.05,]
nrow(test) #293
test<-hawdat[hawdat$month5vs4_haw_FDR < 0.05,]
nrow(test) #13
test<-hawdat[hawdat$month6vs5_haw_FDR < 0.05,]
nrow(test) #1











