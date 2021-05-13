setwd("~/Dropbox/cerasi/annotation/trinotate/")

library(tidyverse)

#cleaning up trinotate annotation file so that it is sane
#annotation
annotation<-read.table('trinotate_cerasi_evigenes_annotation_report_reformatnames.txt', sep="\t", header=TRUE, row.names=NULL, quote="", na.strings=".",stringsAsFactors = F) 
head(annotation)
colnames(annotation)
colnames(annotation)[which(names(annotation) == "X.nesi.project.uoo00108.bin.blast_db.dmel.all.translation.r6.37.fasta_BLASTX")] <- "dmel_BLASTX"
colnames(annotation)[which(names(annotation) == "X.nesi.project.uoo00108.bin.blast_db.dmel.all.translation.r6.37.fasta_BLASTP")] <- "dmel_BLASTP"
colnames(annotation)[which(names(annotation) == "X.nesi.project.uoo00108.bin.blast_db.GCF_013731165.1_Rhpom_1.0_protein.faa_BLASTX")] <- "rpom_BLASTX"
colnames(annotation)[which(names(annotation) == "X.nesi.project.uoo00108.bin.blast_db.GCF_013731165.1_Rhpom_1.0_protein.faa_BLASTP")] <- "rpom_BLASTP"

#so the extras must be duplicates in the annotation file
duplicate_rows<-annotation %>% group_by(gene_id) %>% filter(n()>1)

#its where transdecoder has annotation >1 protein for a gene and its returning more than one match.

#lets simplyfy the annotations and then collapse into one annotation per gene
annotationsSimplified = annotation %>% mutate(sprot_BLASTX_gene=gsub('_.*',"",sprot_Top_BLASTX_hit))
annotationsSimplified$sprot_BLASTX_gene = gsub('"','',annotationsSimplified$sprot_BLASTX_gene)

annotationsSimplified = annotationsSimplified %>% mutate(sprot_BLASTX_name = str_extract(sprot_Top_BLASTX_hit, "Full=.*?;"))
annotationsSimplified$sprot_BLASTX_name = gsub('Full=','',annotationsSimplified$sprot_BLASTX_name)
annotationsSimplified$sprot_BLASTX_name = gsub(';','',annotationsSimplified$sprot_BLASTX_name)
annotationsSimplified$sprot_BLASTX_name = gsub(' \\{.*\\}','',annotationsSimplified$sprot_BLASTX_name)

annotationsSimplified = annotationsSimplified %>% mutate(sprot_BLASTP_gene=gsub('_.*',"",sprot_Top_BLASTP_hit))
annotationsSimplified$sprot_BLASTP_gene = gsub('"','',annotationsSimplified$sprot_BLASTP_gene)

annotationsSimplified = annotationsSimplified %>% mutate(sprot_BLASTP_name = str_extract(sprot_Top_BLASTP_hit, "Full=.*?;"))
annotationsSimplified$sprot_BLASTP_name = gsub('Full=','',annotationsSimplified$sprot_BLASTP_name)
annotationsSimplified$sprot_BLASTP_name = gsub(';','',annotationsSimplified$sprot_BLASTP_name)
annotationsSimplified$sprot_BLASTP_name = gsub(' \\{.*\\}','',annotationsSimplified$sprot_BLASTP_name)

annotationsSimplified = annotationsSimplified %>% mutate(dmel_Top_BLASTX=gsub('\\^.*',"",dmel_BLASTX))
annotationsSimplified$dmel_BLASTX = gsub('"','',annotationsSimplified$dmel_BLASTX)

annotationsSimplified = annotationsSimplified %>% mutate(dmel_Top_BLASTP=gsub('\\^.*',"",dmel_BLASTP))
annotationsSimplified$dmel_BLASTP = gsub('"','',annotationsSimplified$dmel_BLASTP)

annotationsSimplified = annotationsSimplified %>% mutate(rpom_Top_BLASTX=gsub('\\^.*',"",rpom_BLASTX))
annotationsSimplified$rpom_BLASTX = gsub('"','',annotationsSimplified$rpom_BLASTX)

annotationsSimplified = annotationsSimplified %>% mutate(rpom_Top_BLASTP=gsub('\\^.*',"",rpom_BLASTP))
annotationsSimplified$rpom_BLASTP = gsub('"','',annotationsSimplified$rpom_BLASTP)

head(annotationsSimplified) 
colnames(annotationsSimplified)
annotationsSimplified = annotationsSimplified %>% mutate(GO_code=gsub("\\^.*","",gene_ontology_blast))

annotationsSimplified = annotationsSimplified %>%mutate(GO_function=gsub(".*\\^","",gene_ontology_blast))

annotationsSimplified = annotationsSimplified %>% mutate(Kegg_KO=gsub(".*`","",Kegg))

annotationsSimplified = annotationsSimplified %>% mutate(Kegg_rno=gsub("`.*","",Kegg))
annotationsSimplified$Kegg_rno = gsub('KEGG:','',annotationsSimplified$Kegg_rno)

annotationsSimplified = annotationsSimplified %>% mutate(RNAMMER_clean=gsub("\\^.*","",RNAMMER))

head(annotationsSimplified)
colnames(annotationsSimplified)
annotationsClean<-annotationsSimplified[c(1,21:33,5:6)]
head(annotationsClean)

#need to get common names of dmel and rpom in for annotaions
dmel_info<-read.table('../fbgn_fbtr_fbpp_fb_2020_06.tsv',header=T,row.names=NULL,sep='\t',quote='',stringsAsFactors = F)
#rows we need
head(dmel_info)

dmel_info$FlyBase_FBtr<-NULL
annotationsClean<-left_join(annotationsClean,dmel_info,by=c('dmel_Top_BLASTX'='FlyBase_FBpp'))
colnames(annotationsClean)[colnames(annotationsClean)=="FlyBase_FBgn"] <- "dmel_Top_BLASTX_FBgn"
colnames(annotationsClean)[colnames(annotationsClean)=="dmel_Top_BLASTX"] <- "dmel_Top_BLASTX_FBpp"
colnames(annotationsClean)

annotationsClean<-left_join(annotationsClean,dmel_info,by=c('dmel_Top_BLASTP'='FlyBase_FBpp'))
colnames(annotationsClean)[colnames(annotationsClean)=="FlyBase_FBgn"] <- "dmel_Top_BLASTP_FBgn"
colnames(annotationsClean)[colnames(annotationsClean)=="dmel_Top_BLASTP"] <- "dmel_Top_BLASTP_FBpp"
colnames(annotationsClean)
head(annotationsClean)
#Rpom
rpom_info<-read.table('../GCF_013731165.1_Rhpom_1.0_feature_table.txt',header=T,row.names=NULL,sep='\t',quote='',stringsAsFactors = F)
head(rpom_info)
rpom_info_useful<-rpom_info %>% filter(str_detect(product_accession, "^XP_"))
length(unique(rpom_info_useful$product_accession))
length(rpom_info_useful$product_accession)

rpom_info_useful<-rpom_info_useful %>% select(product_accession,symbol)

head(rpom_info_useful)
colnames(rpom_info_useful)[colnames(rpom_info_useful)=="product_accession"] <- "rpom_protien"
colnames(rpom_info_useful)[colnames(rpom_info_useful)=="symbol"] <- "rpom_LOC"

annotationsClean<-left_join(annotationsClean,rpom_info_useful,by=c('rpom_Top_BLASTX'='rpom_protien'))
colnames(annotationsClean)[colnames(annotationsClean)=="rpom_LOC"] <- "rpom_blastX_LOC"
colnames(annotationsClean)
head(annotationsClean)

annotationsClean<-left_join(annotationsClean,rpom_info_useful,by=c('rpom_Top_BLASTP'='rpom_protien'))
colnames(annotationsClean)[colnames(annotationsClean)=="rpom_LOC"] <- "rpom_blastP_LOC"
colnames(annotationsClean)
head(annotationsClean)

#okay so now I think we just have to deal with the fact there is around 1000 genes that have more than one protein prediction for them

#tis a bit of a shit way of doing this but Im just going to collapse into a string
#this means you cant blindly use columns in enrichment analyses as there maybe more than one code in a cell
colnames(annotationsClean)
annotationsCleanSingle <- annotationsClean %>% group_by(gene_id) %>%
  summarize(sprot_BLASTX_gene = paste(unique(sprot_BLASTX_gene[!is.na(sprot_BLASTX_gene)]), collapse = ","),
            sprot_BLASTX_name = paste(unique(sprot_BLASTX_name[!is.na(sprot_BLASTX_name)]), collapse = ","),        
            sprot_BLASTP_gene = paste(unique(sprot_BLASTP_gene[!is.na(sprot_BLASTP_gene)]), collapse = ","),
            sprot_BLASTP_name = paste(unique(sprot_BLASTP_name[!is.na(sprot_BLASTP_name)]), collapse = ","),
            dmel_Top_BLASTX_FBpp = paste(unique(dmel_Top_BLASTX_FBpp[!is.na(dmel_Top_BLASTX_FBpp)]), collapse = ","),
            dmel_Top_BLASTP_FBpp = paste(unique(dmel_Top_BLASTP_FBpp[!is.na(dmel_Top_BLASTP_FBpp)]), collapse = ","),
            rpom_Top_BLASTX = paste(unique(rpom_Top_BLASTX[!is.na(rpom_Top_BLASTX)]), collapse = ","),
            rpom_Top_BLASTP = paste(unique(rpom_Top_BLASTP[!is.na(rpom_Top_BLASTP)]), collapse = ","),
            GO_code = paste(unique(GO_code[!is.na(GO_code)]), collapse = ","),
            GO_function = paste(unique(GO_function[!is.na(GO_function)]), collapse = ","),
            Kegg_KO = paste(unique(Kegg_KO[!is.na(Kegg_KO)]), collapse = ","),
            Kegg_rno = paste(unique(Kegg_rno[!is.na(Kegg_rno)]), collapse = ","),
            RNAMMER_clean = paste(unique(RNAMMER_clean[!is.na(RNAMMER_clean)]), collapse = ","),
            protien_id = paste(unique(prot_id[!is.na(prot_id)]), collapse = ","),
            protien_coords = paste(unique(prot_coords[!is.na(prot_coords)]), collapse = ","),
            dmel_Top_BLASTX_FBgn = paste(unique(dmel_Top_BLASTX_FBgn[!is.na(dmel_Top_BLASTX_FBgn)]), collapse = ","),
            dmel_Top_BLASTP_FBg = paste(unique(dmel_Top_BLASTP_FBgn[!is.na(dmel_Top_BLASTP_FBgn)]), collapse = ","),
            rpom_blastX_LOC = paste(unique(rpom_blastX_LOC[!is.na(rpom_blastX_LOC)]), collapse = ","),
            rpom_blastP_LOC = paste(unique(rpom_blastP_LOC[!is.na(rpom_blastP_LOC)]), collapse = ","))

annotationsCleanSingle <- annotationsCleanSingle %>% mutate_all(na_if,"")
write.table(annotationsCleanSingle,'trinotate_cerasi_evigenes_annotation_report_reformatnames_onerow_cleaned.txt',row.names=F,sep='\t',quote=F,na='NA')

