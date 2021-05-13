#annotation 2021
#cerasi
drosophila<-read.table('../../../../annotation/final_clean_files/evigenes.Oct6.incAllSampleTrin.OkayAlt.picardformat.fa.cleanedNames.maintranscriptonly_annotation_Dros.merge.cutdown.txt',header=TRUE, row.names =  NULL, quote ="", sep='\t')
rpomonella<-read.table('../../../../annotation/final_clean_files/evigenes.Oct6.incAllSampleTrin.OkayAlt.picardformat.fa.cleanedNames.maintranscriptonly_annotation_Rpom.merge.cutdown.txt',header=TRUE, row.names =NULL,quote="",sep='\t')
colnames(drosophila) <- paste("Dros", colnames(drosophila), sep = "_")

colnames(rpomonella) <- paste("Rpom", colnames(rpomonella), sep = "_")


pom_dros<-full_join(drosophila,rpomonella,by=c('Dros_blastx_Cer_transcript'='Rpom_blastx_Cer_transcript'))
colnames(pom_dros)
pom_dros_cutdown<-pom_dros %>% select(-Dros_blastx_evalue,-Dros_tblastn_evalue,-Dros_blastx_blastx_yn,-Dros_tblastn_tblastn_yn,-Rpom_blastx_evalue,-Rpom_tblastn_evalue,-Rpom_blastx_blastx_yn,-Rpom_tblastn_tblastn_yn)

head(pom_dros_cutdown)


write.table(pom_dros_cutdown,'../../../../annotation/final_clean_files/evigenes.Oct6.incAllSampleTrin.OkayAlt.picardformat.fa.cleanedNames.maintranscriptonly_annotation_Rpom_Dros.merge.cutdown.txt',row.names=F,quote=F,sep='\t')
