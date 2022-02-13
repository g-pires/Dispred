suppressMessages(library('Xmisc'))
suppressMessages(library(bigsnpr))
suppressMessages(library(bigstatsr))
suppressMessages(library(foreach))
suppressMessages(library(bigreadr))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))

##########################################
#Code inspired by : https://github.com/privefl/UKBiobank/blob/31a563fecda17d270d7e6cc6310f59a9c99b960c/0-download-genotype-data.R
#
#This code is to extract the genotypes only, 
#with the use of plink to merge all bed files
#into one to be read by the function snp_readBed2() 
#
#
##########################################

BGEN_PATH_ext = "ukb_imp_chr{chr}_v3.bgen"

IND_PATH='new_ind.txt'

SAMPLE_PATH = 'ukb42260_imp_chr1_v3_s487298.sample'

BGEN_OUTPUT_PATH = ''

PATH = ""

ind<-read.csv(IND_PATH, sep=',')

final_ind<-na.omit(ind$eid[ind$X22021.0.0==0])

sample <- fread2(SAMPLE_PATH)

ind_keep<-sample[sample$ID_2 %in% final_ind,]

geno_table<-read.table(paste0(PATH, 'all_genotyped_snps'))

liste_geno_snp_id<-split(with(geno_table, paste(V1, V3, V4, V5, sep='_')), factor(geno_table$V1, levels=1:22, ordered = T))

bgen_file<-snp_readBGEN(bgenfiles = glue::glue(BGEN_PATH_ext, chr=1:22), backingfile = paste0(BGEN_OUTPUT_PATH, 'all_chr_geno_v3'), read_as="dosage", list_snp_id = liste_geno_snp_id, ind_row=rows_along(ind_keep))

bgen<-snp_attach(bgen_file)

G<-bgen$genotypes
dim(G)
