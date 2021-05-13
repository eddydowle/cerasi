#moving annotations together

setwd('~/Documents/Cerasi/Cerosi/annotating/2021/')

#to rphom
pom_map<-read.table('GCF_013731165.1_Rhpom_1.0_feature_table.txt',header=T,sep='\t',quote="")

pom_map_xp<-pom_map %>% filter(str_detect(product_accession, "^XP_"))

pom_map_xp%>% group_by(product_accession) %>% summarize(count=n())

#30758 uniq proteins

#blastx
cer_pom_blastx<-read.table('evigenes.Oct6.incAllSampleTrin.OkayAlt.picardformat.fa.cleanedNames.maintranscriptonly_blastx_Rhpom_tab',header=F,sep='\t',quote="")
head(cer_pom_blastx)

names(cer_pom_blastx) <- c("Cer_transcript","Pom_transcript","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")
head(cer_pom_blastx)

library(tidyverse)
names(pom_map_xp)
cer_pom_blastx_loc<-left_join(cer_pom_blastx,pom_map_xp,by=c('Pom_transcript'='product_accession'))

head(cer_pom_blastx_loc)
tail(cer_pom_blastx_loc)
cer_pom_blastx_loc%>% group_by(Cer_transcript) %>% summarize(count=n())
#out of 109867 transcripts 102666 are annotated to pomonella
cer_pom_blastx_loc%>% group_by(Cer_transcript,symbol) %>% summarize(count=n())
#40087071

#so we need to remove all the things not the top 


filter_lowestevalue<-cer_pom_blastx_loc %>% 
  group_by(Cer_transcript) %>% 
  summarise(code = min(evalue))


cer_pom_blastx_loc_lowestevalue<-right_join(cer_pom_blastx_loc,filter_lowestevalue,by=c('Cer_transcript'='Cer_transcript', 'evalue'='code'))
head(test)
head(test)

cer_pom_blastx_loc_lowestevalue%>% group_by(Cer_transcript,symbol) %>% summarize(count=n())
139206-102666
#36540K have two genes with identical evalues


#collapse into one line per gene
#take useful columns
#Cer_transcript, symbol, name, evalue, bitscore

cer_pom_blastx_loc_lowestevalue$blastx_yn<-'yes'
#test<-cer_pom_blastx_loc_lowestevalue %>% select(.,Cer_transcript, symbol, name, evalue) %>% distinct(Cer_transcript,symbol,name,evalue)
head (test)


#cer_pom_blastx_loc_lowestevalue_collapse <- test %>% group_by(Cer_transcript) %>%
#  summarize(symbol = paste(unique(symbol[!is.na(symbol)]), collapse = ","),
#            name = paste(unique(name[!is.na(name)]), collapse = ","),        
#            evalue = paste(unique(evalue[!is.na(evalue)]), collapse = ","))

###annotationsCleanSingle <- annotationsCleanSingle %>% mutate_all(na_if,"")

#now lets look at the tblastn

cer_pom_tblastn<-read.table('evigenes.Oct6.incAllSampleTrin.OkayAlt.picardformat.fa.cleanedNames.maintranscriptonly_tblastn_Rhpom_tab',header=F,sep='\t',quote="")
head(cer_pom_tblastn)

names(cer_pom_tblastn) <- c("Pom_transcript","Cer_transcript","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")
head(cer_pom_tblastn)
cer_pom_tblastn_loc<-left_join(cer_pom_tblastn,pom_map_xp,by=c('Pom_transcript'='product_accession'))

filter_lowestevalue<-cer_pom_tblastn_loc %>% 
  group_by(Pom_transcript) %>% 
  summarise(code = min(evalue))

cer_pom_tblastn_loc_lowestevalue<-right_join(cer_pom_tblastn_loc,filter_lowestevalue,by=c('Pom_transcript'='Pom_transcript', 'evalue'='code'))
head(cer_pom_tblastn_loc_lowestevalue)
cer_pom_tblastn_loc_lowestevalue%>% group_by(Pom_transcript,symbol) %>% summarize(count=n())
#30532

#test<-cer_pom_tblastn_loc_lowestevalue %>% select(.,Cer_transcript, symbol, name, evalue) %>% distinct(Cer_transcript,symbol,name,evalue)
head (test)


#cer_pom_tblastn_loc_lowestevalue_collapse <- test %>% group_by(Cer_transcript) %>%
#  summarize(symbol = paste(unique(symbol[!is.na(symbol)]), collapse = ","),
#            name = paste(unique(name[!is.na(name)]), collapse = ","),        
#            evalue = paste(unique(evalue[!is.na(evalue)]), collapse = ","))

####annotationsCleanSingle <- annotationsCleanSingle %>% mutate_all(na_if,"")


head(cer_pom_blastx_loc_lowestevalue_collapse)
head(cer_pom_blastx_loc_lowestevalue)
head(cer_pom_tblastn_loc_lowestevalue)

cer_pom_blastx_loc_lowestevalue$blastx_yn<-'yes'
colnames(cer_pom_blastx_loc_lowestevalue) <- paste("blastx", colnames(cer_pom_blastx_loc_lowestevalue), sep = "_")
head(cer_pom_blastx_loc_lowestevalue)

colnames(cer_pom_tblastn_loc_lowestevalue) <- paste("tblastn", colnames(cer_pom_tblastn_loc_lowestevalue), sep = "_")

colnames(cer_pom_blastx_loc_lowestevalue)
colnames(cer_pom_tblastn_loc_lowestevalue)

test<-full_join(cer_pom_blastx_loc_lowestevalue,cer_pom_tblastn_loc_lowestevalue,by=c('blastx_Cer_transcript'='tblastn_Cer_transcript','blastx_symbol'='tblastn_symbol'))
head(test)

tail(test)
colnames(test)
test<-test %>% mutate(Recprical_blast = if_else(tblastn_tblastn_yn =='yes' & blastx_blastx_yn =='yes', 1, 0))
head(test)
tail(test)
colnames(test)
test_cutdown<-test %>% select(.,blastx_Cer_transcript,blastx_symbol,blastx_name,blastx_evalue,tblastn_name,tblastn_evalue,blastx_blastx_yn,tblastn_tblastn_yn,Recprical_blast)


#so we need to make them uniq
head(test_cutdown)
colnames(test_cutdown)
test_cutdown_uniq<-test_cutdown %>%  distinct()

#there is duplicates of gene description under the same loci
#test_cutdown_uniq_b<-test_cutdown %>%  distinct(blastx_Cer_transcript,blastx_symbol,tblastn_name,blastx_name)

test_cutdown_uniq_collapse <- test %>% group_by(blastx_Cer_transcript) %>%
  summarize(symbol = paste(unique(blastx_symbol[!is.na(blastx_symbol)]), collapse = ","),
            blastx_name = paste(unique(blastx_name[!is.na(blastx_name)]), collapse = ","),        
            blastx_evalue = paste(unique(blastx_evalue[!is.na(blastx_evalue)]), collapse = ","),
            tblastn_name = paste(unique(tblastn_name[!is.na(tblastn_name)]), collapse = ","),
            tblastn_evalue = paste(unique(tblastn_evalue[!is.na(tblastn_evalue)]), collapse = ","),
            blastx_blastx_yn = paste(unique(blastx_blastx_yn[!is.na(blastx_blastx_yn)]), collapse = ","),
            tblastn_tblastn_yn = paste(unique(tblastn_tblastn_yn[!is.na(tblastn_tblastn_yn)]), collapse = ","),
            Recprical_blast = paste(unique(Recprical_blast[!is.na(Recprical_blast)]), collapse = ","))

head(test_cutdown_uniq_collapse)
nrow(test_cutdown_uniq_collapse)
length(unique(test_cutdown_uniq_collapse$blastx_Cer_transcript))


write.table(test_cutdown_uniq_collapse ,'evigenes.Oct6.incAllSampleTrin.OkayAlt.picardformat.fa.cleanedNames.maintranscriptonly_annotation_Rpom.merge.cutdown.txt',row.names=F,quote=F,sep="\t")


#need another one where I jjust take reciprocal blast hits

head(test_cutdown_uniq)

#so where Recprical_blast == 1

test_cutdown_reciprocal<-test_cutdown_uniq %>% filter(Recprical_blast==1)

#still need to likely reduce down to single gene
test_cutdown_reciprocal_collapse <- test_cutdown_reciprocal %>% group_by(blastx_Cer_transcript) %>%
  summarize(symbol = paste(unique(blastx_symbol[!is.na(blastx_symbol)]), collapse = ","),
            blastx_name = paste(unique(blastx_name[!is.na(blastx_name)]), collapse = ","),        
            blastx_evalue = paste(unique(blastx_evalue[!is.na(blastx_evalue)]), collapse = ","),
            tblastn_name = paste(unique(tblastn_name[!is.na(tblastn_name)]), collapse = ","),
            tblastn_evalue = paste(unique(tblastn_evalue[!is.na(tblastn_evalue)]), collapse = ","),
            blastx_blastx_yn = paste(unique(blastx_blastx_yn[!is.na(blastx_blastx_yn)]), collapse = ","),
            tblastn_tblastn_yn = paste(unique(tblastn_tblastn_yn[!is.na(tblastn_tblastn_yn)]), collapse = ","),
            Recprical_blast = paste(unique(Recprical_blast[!is.na(Recprical_blast)]), collapse = ","))

head(test_cutdown_reciprocal_collapse)
write.table(test_cutdown_reciprocal_collapse ,'evigenes.Oct6.incAllSampleTrin.OkayAlt.picardformat.fa.cleanedNames.maintranscriptonly_annotation_Rpom.reciprocal.cutdown.txt',row.names=F,quote=F,sep="\t")

#########################
#doing it for drosophila#
#########################

#to droso
droso_gene_map<-read.table('fbgn_annotation_ID_fb_2020_06.tsv',header=T,sep='\t',quote="")

droso_pp_map<-read.table('fbgn_fbtr_fbpp_fb_2020_06.tsv',header=T,sep='\t',quote="")

head(droso_gene_map)
head(droso_pp_map)

drosophila_gene_map<-left_join(droso_pp_map,droso_gene_map,by=c('FlyBase_FBgn'='primary_FBgn'))
head(drosophila_gene_map)


#blastx
cer_dros_blastx<-read.table('evigenes.Oct6.incAllSampleTrin.OkayAlt.picardformat.fa.cleanedNames.maintranscriptonly_blastx_dmel_6.37_tab',header=F,sep='\t',quote="")
head(cer_dros_blastx)

names(cer_dros_blastx) <- c("Cer_transcript","dros_protein","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")
head(cer_dros_blastx)


library(tidyverse)
names(drosophila_gene_map)
cer_dros_blastx_loc<-left_join(cer_dros_blastx,drosophila_gene_map,by=c('dros_protein'='FlyBase_FBpp'))

head(cer_dros_blastx_loc)
tail(cer_dros_blastx_loc)
cer_dros_blastx_loc%>% group_by(Cer_transcript) %>% summarize(count=n())
#out of 109867 transcripts 92685 are annotated to drosophila
cer_dros_blastx_loc%>% group_by(Cer_transcript,FlyBase_FBgn) %>% summarize(count=n())
#567634

#so we need to remove all the things not the top 


filter_lowestevalue<-cer_dros_blastx_loc %>% 
  group_by(Cer_transcript) %>% 
  summarise(code = min(evalue))

head(cer_dros_blastx_loc)

cer_dros_blastx_loc_lowestevalue<-right_join(cer_dros_blastx_loc,filter_lowestevalue,by=c('Cer_transcript'='Cer_transcript', 'evalue'='code'))
head(cer_dros_blastx_loc_lowestevalue)
head(cer_dros_blastx_loc)
cer_dros_blastx_loc_lowestevalue%>% group_by(Cer_transcript,FlyBase_FBgn) %>% summarize(count=n())
175299-94757
#80542K have two genes with identical evalues


#collapse into one line per gene
#take useful columns
#Cer_transcript, flybase_fbgn

cer_dros_blastx_loc_lowestevalue$blastx_yn<-'yes'


#now lets look at the tblastn

cer_dros_tblastn<-read.table('evigenes.Oct6.incAllSampleTrin.OkayAlt.picardformat.fa.cleanedNames.maintranscriptonly_tblastn_dmel_6.37_tab',header=F,sep='\t',quote="")
head(cer_dros_tblastn)

names(cer_dros_tblastn) <- c("dros_protein","Cer_transcript","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")
head(cer_dros_tblastn)
head(drosophila_gene_map)
cer_dros_tblastn_loc<-left_join(cer_dros_tblastn,drosophila_gene_map,by=c('dros_protein'='FlyBase_FBpp'))

filter_lowestevalue<-cer_dros_tblastn_loc %>% 
  group_by(dros_protein) %>% 
  summarise(code = min(evalue))



cer_dros_tblastn_loc_lowestevalue<-right_join(cer_dros_tblastn_loc,filter_lowestevalue,by=c('dros_protein'='dros_protein', 'evalue'='code'))
head(cer_dros_tblastn_loc_lowestevalue)
cer_dros_tblastn_loc_lowestevalue%>% group_by(dros_protein,FlyBase_FBgn) %>% summarize(count=n())
#30063

#
head(cer_dros_tblastn_loc_lowestevalue)

cer_dros_tblastn_loc_lowestevalue$tblastn_yn<-'yes'
colnames(cer_dros_blastx_loc_lowestevalue) <- paste("blastx", colnames(cer_dros_blastx_loc_lowestevalue), sep = "_")
head(cer_dros_blastx_loc_lowestevalue)

colnames(cer_dros_tblastn_loc_lowestevalue) <- paste("tblastn", colnames(cer_dros_tblastn_loc_lowestevalue), sep = "_")

colnames(cer_dros_blastx_loc_lowestevalue)
colnames(cer_dros_tblastn_loc_lowestevalue)

test<-full_join(cer_dros_blastx_loc_lowestevalue,cer_dros_tblastn_loc_lowestevalue,by=c('blastx_Cer_transcript'='tblastn_Cer_transcript','blastx_FlyBase_FBgn'='tblastn_FlyBase_FBgn'))
head(test)

tail(test)
colnames(test)
test<-test %>% mutate(Recprical_blast = if_else(tblastn_tblastn_yn =='yes' & blastx_blastx_yn =='yes', 1, 0))
head(test)
tail(test)
colnames(test)
test_cutdown<-test %>% select(.,blastx_Cer_transcript,blastx_FlyBase_FBgn,blastx_gene_symbol,blastx_evalue,tblastn_gene_symbol,tblastn_evalue,blastx_blastx_yn,tblastn_tblastn_yn,Recprical_blast)

#so we need to make them uniq
head(test_cutdown)
colnames(test_cutdown)
test_cutdown_uniq<-test_cutdown %>%  distinct()

#there is duplicates of gene description under the same loci
#test_cutdown_uniq_b<-test_cutdown %>%  distinct(blastx_Cer_transcript,blastx_symbol,tblastn_name,blastx_name)
colnames(test_cutdown_uniq)
test_cutdown_uniq_collapse <- test_cutdown_uniq %>% group_by(blastx_Cer_transcript) %>%
  summarize(symbol = paste(unique(blastx_FlyBase_FBgn[!is.na(blastx_FlyBase_FBgn)]), collapse = ","),
            blastx_name = paste(unique(blastx_gene_symbol[!is.na(blastx_gene_symbol)]), collapse = ","),        
            blastx_evalue = paste(unique(blastx_evalue[!is.na(blastx_evalue)]), collapse = ","),
            tblastn_name = paste(unique(tblastn_gene_symbol[!is.na(tblastn_gene_symbol)]), collapse = ","),
            tblastn_evalue = paste(unique(tblastn_evalue[!is.na(tblastn_evalue)]), collapse = ","),
            blastx_blastx_yn = paste(unique(blastx_blastx_yn[!is.na(blastx_blastx_yn)]), collapse = ","),
            tblastn_tblastn_yn = paste(unique(tblastn_tblastn_yn[!is.na(tblastn_tblastn_yn)]), collapse = ","),
            Recprical_blast = paste(unique(Recprical_blast[!is.na(Recprical_blast)]), collapse = ","))

head(test_cutdown_uniq_collapse)
nrow(test_cutdown_uniq_collapse)
length(unique(test_cutdown_uniq_collapse$blastx_Cer_transcript))


write.table(test_cutdown_uniq_collapse ,'evigenes.Oct6.incAllSampleTrin.OkayAlt.picardformat.fa.cleanedNames.maintranscriptonly_annotation_Dros.merge.cutdown.txt',row.names=F,quote=F,sep="\t")



#need another one where I jjust take reciprocal blast hits

head(test_cutdown_uniq)

#so where Recprical_blast == 1

test_cutdown_reciprocal<-test_cutdown_uniq %>% filter(Recprical_blast==1)

#still need to likely reduce down to single gene
test_cutdown_reciprocal_collapse <- test_cutdown_reciprocal %>% group_by(blastx_Cer_transcript) %>%
  summarize(symbol = paste(unique(blastx_FlyBase_FBgn[!is.na(blastx_FlyBase_FBgn)]), collapse = ","),
            blastx_name = paste(unique(blastx_gene_symbol[!is.na(blastx_gene_symbol)]), collapse = ","),        
            blastx_evalue = paste(unique(blastx_evalue[!is.na(blastx_evalue)]), collapse = ","),
            tblastn_name = paste(unique(tblastn_gene_symbol[!is.na(tblastn_gene_symbol)]), collapse = ","),
            tblastn_evalue = paste(unique(tblastn_evalue[!is.na(tblastn_evalue)]), collapse = ","),
            blastx_blastx_yn = paste(unique(blastx_blastx_yn[!is.na(blastx_blastx_yn)]), collapse = ","),
            tblastn_tblastn_yn = paste(unique(tblastn_tblastn_yn[!is.na(tblastn_tblastn_yn)]), collapse = ","),
            Recprical_blast = paste(unique(Recprical_blast[!is.na(Recprical_blast)]), collapse = ","))

head(test_cutdown_reciprocal_collapse)
write.table(test_cutdown_reciprocal_collapse ,'evigenes.Oct6.incAllSampleTrin.OkayAlt.picardformat.fa.cleanedNames.maintranscriptonly_annotation_Dros.reciprocal.cutdown.txt',row.names=F,quote=F,sep="\t")
