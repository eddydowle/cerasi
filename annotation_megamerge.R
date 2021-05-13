#merging it all together

library(tidyverse)
setwd("~/Dropbox/cerasi/annotation/")

#entap
entap<-read.table('entap/final_results/final_annotations_no_contam_lvl1_reformatnames.tsv',header=T,row.names=NULL,sep="\t",quote="")
colnames(entap)
#cleanup entap
entap_reduced<-entap %>% select(Query.Sequence,Subject.Sequence,E.Value,Description, Species,EggNOG.Predicted.Gene,EggNOG.Description, EggNOG.Seed.Ortholog, EggNOG.GO.Biological,EggNOG.GO.Cellular,EggNOG.GO.Molecular)
head(entap_reduced)
colnames(entap_reduced) <- paste("EnTAP", colnames(entap_reduced), sep = "_")
write.table(entap_reduced,'entap/final_results/final_annotations_no_contam_lvl1_reformatnames_clean.tsv',row.names=F,quote=F,sep='\t')

#trinotate
trinotate<-read.table('trinotate/trinotate_cerasi_evigenes_annotation_report_reformatnames_onerow_cleaned.txt',header=T,row.names=NULL,sep='\t',quote="")
colnames(trinotate) <- paste("Trinotate", colnames(trinotate), sep = "_")
head(trinotate)

#blasts
blast_reciprocal_dros<-read.table('~/Documents/Cerasi/Cerosi/annotating/2021/evigenes.Oct6.incAllSampleTrin.OkayAlt.picardformat.fa.cleanedNames.maintranscriptonly_annotation_Dros.reciprocal.cutdown.txt',header=T,row.names=NULL,sep='\t',quote="")
blast_reciprocal_rpom<-read.table('~/Documents/Cerasi/Cerosi/annotating/2021/evigenes.Oct6.incAllSampleTrin.OkayAlt.picardformat.fa.cleanedNames.maintranscriptonly_annotation_Rpom.reciprocal.cutdown.txt',header=T,row.names=NULL,sep='\t',quote="")


colnames(blast_reciprocal_dros) <- paste("RB_dros", colnames(blast_reciprocal_dros), sep = "_")
colnames(blast_reciprocal_rpom) <- paste("RB_rpom", colnames(blast_reciprocal_rpom), sep = "_")

colnames(blast_reciprocal_dros)
colnames(blast_reciprocal_rpom)
reciprocal_blast<-full_join(blast_reciprocal_dros,blast_reciprocal_rpom,by=c('RB_dros_blastx_Cer_transcript'='RB_rpom_blastx_Cer_transcript'))
head(reciprocal_blast)

reciprocal_blast_clean<-reciprocal_blast %>% select(-RB_dros_blastx_blastx_yn,-RB_dros_tblastn_tblastn_yn,-RB_rpom_tblastn_tblastn_yn,-RB_rpom_blastx_blastx_yn,-RB_rpom_Recprical_blast,-RB_dros_Recprical_blast)

head(reciprocal_blast_clean)
names(reciprocal_blast_clean)[names(reciprocal_blast_clean) == "RB_dros_blastx_Cer_transcript"] <- "Cer_transcript"
write.table(reciprocal_blast_clean,'Reciprocal_blast_drosophila_rpomonella_cerasi.txt',row.names=F,quote=F,sep='\t')

#full blast results

blast_dros<-read.table('~/Documents/Cerasi/Cerosi/annotating/2021/evigenes.Oct6.incAllSampleTrin.OkayAlt.picardformat.fa.cleanedNames.maintranscriptonly_annotation_Dros.merge.cutdown.txt',header=T,row.names=NULL,sep='\t',quote="")
blast_rpom<-read.table('~/Documents/Cerasi/Cerosi/annotating/2021/evigenes.Oct6.incAllSampleTrin.OkayAlt.picardformat.fa.cleanedNames.maintranscriptonly_annotation_Rpom.merge.cutdown.txt',header=T,row.names=NULL,sep='\t',quote="")


colnames(blast_dros) <- paste("dros", colnames(blast_dros), sep = "_")
colnames(blast_rpom) <- paste("rpom", colnames(blast_rpom), sep = "_")

colnames(blast_dros)
colnames(blast_rpom)
full_blast<-full_join(blast_dros,blast_rpom,by=c('dros_blastx_Cer_transcript'='rpom_blastx_Cer_transcript'))
head(full_blast)
colnames(full_blast)
full_blast_clean<-full_blast %>% select(-dros_blastx_blastx_yn,-dros_tblastn_tblastn_yn,-rpom_tblastn_tblastn_yn,-rpom_blastx_blastx_yn,-rpom_Recprical_blast,-dros_Recprical_blast)

head(full_blast_clean)
names(full_blast_clean)[names(full_blast_clean) == "dros_blastx_Cer_transcript"] <- "Cer_transcript"
write.table(reciprocal_blast_clean,'Full_blast_drosophila_rpomonella_cerasi.txt',row.names=F,quote=F,sep='\t')
