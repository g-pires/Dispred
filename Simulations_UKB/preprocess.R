suppressMessages(library('Xmisc'))
suppressMessages(library(bigsnpr))
suppressMessages(library(bigstatsr))
suppressMessages(library(foreach))
suppressMessages(library(bigreadr))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(dbplyr))
suppressMessages(library(biganalytics))


CutBySize <- function(m, block.size, nb = ceiling(m / block.size)) {
  int <- m / nb
  
  upper <- round(1:nb * int)
  lower <- c(1, upper[-nb] + 1)
  size <- c(upper[1], diff(upper))
  
  cbind(lower, upper, size)
}

seq2 <- function(lims) seq(lims["lower"], lims["upper"])

big_aggregate <- function(X, FUN, .combine, block.size = 1e6) {
  intervals <- CutBySize(ncol(X), block.size)
  
  foreach(k = 1:nrow(intervals), .combine = .combine) %do% {
    FUN(X[, seq2(intervals[k, ])])
  }
}


PATH = ''

OUTPUT_PATH=''

BGEN_PATH_ext = 'ukb_imp_chr{chr}_v3.bgen'

SAMPLE_PATH = 'ukb42260_imp_chr1_v3_s487298.sample'

IND_PATH='new_ind.txt'

LDH='LD-Hub_study_information_and_SNP_heritability.csv'

Traits<-c('Type 2 Diabetes', 'Coronary artery disease', 'Inflammatory Bowel Disease (Euro)', 'Alzheimers disease', 'Bipolar disorder', 'Rheumatoid Arthritis', 'Schizophrenia', 'Lung cancer (all)', 'Asthma')

LD_table<-read.csv(LDH)

H2_traits<-sapply(Traits, function(X){
  LD_table$H2[LD_table$Trait==X]
})

heritability = 0.6

k = 0.03

t<-qnorm(1-k)

ind<-read.csv(IND_PATH, sep=',')

final_ind<-na.omit(ind$eid[ind$X22021.0.0==0])

unlink(paste0(paste0(OUTPUT_PATH, '*'), c(".bk", ".rds")))

sample <- fread2(SAMPLE_PATH)

ind_keep<-sample[sample$ID_2 %in% final_ind,]

s<-read.table(paste0(PATH, 'sample_snps'))

sub<-s[sample(nrow(s), size=10000),]

liste_snp_id<-split(with(sub, paste(V1, V3, V4, V5, sep='_')), factor(sub$V1, levels=1:22, ordered = T))

bgen_file_100<-snp_readBGEN(bgenfiles = glue::glue(BGEN_PATH_ext, chr=1:22), backingfile = paste0(OUTPUT_PATH, '10000_imp_chr_v3'), read_as="dosage", list_snp_id = liste_snp_id, ind_row=rows_along(ind_keep))

bgen_100<-snp_attach(paste0(OUTPUT_PATH,'10000_imp_chr_v3.rds'))

G<-bgen_100$genotypes
dim(G)

bgen_100$fam<-sample

beta <- rnorm(ncol(G), 0, sqrt(heritability/ncol(G)))

E<-rnorm(nrow(G), 0, sqrt(1-heritability))

stats<-big_colstats(G)

sd_G<-sqrt(stats$var)

new_beta<-beta*sd_G

n<-nrow(G)
m<-ncol(G)

ind.row<-rows_along(G)
ind.col<-cols_along(G)

length(t(new_beta))

A<-big_prodVec(G, t(new_beta))

P<-A+E

P01<-(P>t)*1

write.table(P01, file=paste0(PATH,'10000_t2d_phenotypes_unrelated_individuals.txt'))

