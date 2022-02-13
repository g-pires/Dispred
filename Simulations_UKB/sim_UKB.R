suppressMessages(library('Xmisc'))
suppressMessages(library(bigsnpr))
suppressMessages(library(bigstatsr))
suppressMessages(library(foreach))
suppressMessages(library(bigreadr))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))


PATH = ''

GENO_PATH = ''

OUTPUT_PATH=''

PHENO_PATH = '10000_t2d_phenotypes_unrelated_individuals.txt'

BGEN_PATH='all_chr_geno_v3.rds'

SAVE_PATH=''

NBCORES=12

P01<-read.table(PHENO_PATH)$x

geno<-snp_attach(BGEN_PATH)

G<-geno$genotypes
dim(G)

geno$fam<-sample

set.seed(1)
ind_train <- sample(nrow(G), 320000)
ind_test <- setdiff(rows_along(G), ind_train)

beta_hat<-big_univLogReg(G, P01[ind_train], ind.train = ind_train, ncores = NBCORES)

beta_hat$p.value <- predict(beta_hat, log10=F)

geno$map$beta<-beta_hat$estim
geno$map$pvalue<-round(beta_hat$p.value, digits=1000)

ind.keep <- snp_clumping(G, infos.chr = as.integer(geno$map$chromosome),
                         ind.row = ind_train,
						 thr.r2 = 0.01,
                         S = abs(beta_hat$score), 
                         ncores = NBCORES/2)

summary(lpS.keep <- -predict(beta_hat)[ind.keep])

#définition des valeurs seuils 
thrs <- seq(0, 4, by = 0.5)

#définition du nombre de prédicteurs en fonction des valeurs seuils 
nb.pred <- sapply(thrs, function(thr) sum(lpS.keep > thr))

prs <- snp_PRS(G, betas.keep = beta_hat$estim[ind.keep], 
               ind.test = ind_test, 
               ind.keep = ind.keep, 
               lpS.keep = lpS.keep, 
               thr.list = thrs)

aucs<-apply(prs, 2, AUC, target = P01[ind_test])

auc_boot<-apply(prs, 2, AUCBoot, target=P01[ind_test])

auc_plot<-qplot(nb.pred, aucs) +
  geom_line() +
  scale_x_log10(breaks = nb.pred) +
  labs(x = "Number of predictors", y = "AUC") +
  theme_bigstatsr()

aucs

write.table(beta_hat, paste0(SAVE_PATH, 'beta_hat'), col.names=T, row.names = T)
write.table(ind.keep, paste0(SAVE_PATH, 'ind.keep'), col.names=T, row.names = T)
write.table(lpS.keep, paste0(SAVE_PATH, 'lpS.keep'), col.names=T, row.names = T)
write.table(thrs, paste0(SAVE_PATH, 'thrs'), col.names=T, row.names = T)
write.table(nb.pred, paste0(SAVE_PATH, 'nb.pred'), col.names=T, row.names = T)
write.table(prs, paste0(SAVE_PATH, 'prs'), col.names=T, row.names = T)
write.table(aucs, paste0(SAVE_PATH, 'aucs'), col.names=T, row.names = T)
write.table(auc_boot, paste0(SAVE_PATH, 'auc_boot'), col.names=T, row.names = T)
ggsave('auc_plot', plot=auc_plot, path=SAVE_PATH, device = 'png')





