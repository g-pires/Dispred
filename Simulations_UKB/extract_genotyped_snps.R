suppressMessages(library('Xmisc'))
suppressMessages(library(bigsnpr))
suppressMessages(library(bigstatsr))
suppressMessages(library(foreach))
suppressMessages(library(bigreadr))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))


PATH = ''

SNP_LISTE_PATH = ''

snps_imp<-read.table(SNP_LISTE_PATH, fill=T)

dim(snps_imp)

imputed<-snps_imp %>% filter(V9==0)

genotyped<-snps_imp %>% filter(V9==1)

dim(imputed)
dim(genotyped)

sample_snps<-snps_imp[sample(nrow(snps_imp), size=20000),]

rownames(imputed)<-NULL
rownames(genotyped)<-NULL
rownames(sample_snps)<-NULL


liste_all_snp<-snps_imp
liste_imputed<-imputed
liste_genotyped<-genotyped
liste_sample_snps<-sample_snps


write.table(snps_imp, file=paste0(PATH,'all_snps'), quote=TRUE, append = FALSE, dec=".", row.names=FALSE, col.names=FALSE, sep ="\t", qmethod = c("escape"))
write.table(imputed, file=paste0(PATH,'imputed_snps'), quote=TRUE, append = FALSE, dec=".", row.names=FALSE, col.names=FALSE, sep ="\t", qmethod = c("escape"))
write.table(genotyped, file=paste0(PATH,'genotyped_snps'), quote=TRUE, append = FALSE, dec=".", row.names=FALSE, col.names=FALSE, sep ="\t", qmethod = c("escape"))
write.table(sample_snps, file=paste0(PATH,'sample_snps'), quote=TRUE, append = FALSE, dec=".", row.names=FALSE, col.names=FALSE, sep ="\t", qmethod = c("escape"))

