suppressMessages(library('Xmisc'))
suppressMessages(library(bigsnpr))
suppressMessages(library(bigstatsr))
suppressMessages(library(foreach))

parser <- ArgumentParser$new()
parser$add_argument('--n',type='numeric',  help='Number of individuals')
parser$add_argument('--s',type='numeric',  help='Number of snps')
parser$add_argument('--p',type='numeric',  help='Proportion of individuals in the training set')
parser$add_argument('--h2',type='numeric', help='Heritability')
parser$add_argument('--k', type='numeric', help='Disease prevalence')


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


#############################################################
#	This code is building a Liability Threshold Model.
#	
#	Usage : 
#			Rscript R_LTM2.R --n NUM --s NUM --p NUM --h2 NUM --k NUM
#
#	Input : 
#			- n : number of individuals
#			- s : number of SNPs
#     - p : number of individuals in the test set 
#     - h2 : heritability
#			- k : disease prevalence
#	
#	Output :
#			- AUCmax_train and mean_AUC_train and tests over 10 simulations  
############################################################

OUTPUT_PATH = ''

nb_individuals = n

nb_snps = s

heritability = h2

k = k

test_prop = p

PATH = paste0('/tmp_dir/', nb_individuals, '_', nb_snps, '_', heritability, '_', k, '_', p, '_')

print(nb_individuals)
print(nb_snps)
print(heritability)
print(k)

liste_AUC_test = c()
liste_cas_total = c()
liste_control_total = c()

for (j in 1:10){
  mafs <- runif(nb_snps, 0.05,0.95)
  beta <- rnorm(nb_snps, 0, sqrt(heritability/nb_snps))
  nb_run=ceiling((nb_individuals*0.25)/(nb_individuals*k))
  bind_df = data.frame()
  liste_control = c()
  liste_cas = c()
  for (i in 1:nb_run){  
    r <- FBM(nrow=nb_individuals, ncol=nb_snps, backingfile=paste0(PATH, 'r'))
    
    r[]<-rbinom(nb_individuals*nb_snps, 2, mafs)

    G <- as_FBM(matrix(r[], nrow=nb_individuals, ncol=nb_snps, byrow=T), backingfile = paste0(PATH, 'G'))
    
    Gstd<-big_aggregate(G, function(X){ 
      as_FBM(scale(X, scale=T, center=T), backingfile = paste0(PATH, 'Gstd'))
    }) 
    
    t <- qnorm(1-k)
    
    A<-big_aggregate(Gstd[[1]], function(X){  
      tcrossprod(beta, X)
    })
    
    A<-t(A[[1]])
    
    E <- rnorm(nb_individuals, 0, sqrt(1-heritability))
    
    P <- A+E
    
    df<-as_FBM(matrix(c(A, E, P), nrow=nb_individuals, ncol=3), backingfile=paste0(PATH, 'df'))
    
    P01_train <- (df[][,3]>t)*1
    
    bind_df1<-FBM(nrow = nb_individuals, ncol = nb_snps+1, backingfile = paste0(PATH, 'bind_df1'))
    bind_df1[]<-data.frame(Gstd[[1]][], P01_train)
    cas=bind_df1[P01_train==1,]
    nb_cas=nrow(cas)
    liste_cas[i]<-nb_cas
    all_control=bind_df1[P01_train==0,]
    control=all_control[sample(nrow(all_control), size=nb_cas),]
    nb_control=nrow(control)
    liste_control[i]<-nb_control
    bind_df2=rbind(cas, control)
    bind_df=rbind(bind_df, bind_df2)
    
    print(i)
    
  }    
  
  nb_cas_total = sum(liste_cas)
  nb_control_total = sum(liste_control)
  
  liste_cas_total[j] <- nb_cas_total
  liste_control_total[j] <- sum(liste_control)
  
  
  beta_hat<-big_aggregate(bind_df[1:nb_snps], function(X){
    big_univLogReg(as_FBM(X, backingfile = paste0(PATH, 'bind_df')), bind_df[,ncol(bind_df)])$estim
  })
  
  r_test <- FBM(nrow=test_prop, ncol=nb_snps, backingfile=paste0(PATH, 'r_test'))
  
  r_test[]<-rbinom(test_prop*nb_snps, 2, mafs)
  
  G_test <- as_FBM(matrix(r_test[], nrow=test_prop, ncol=nb_snps, byrow=T), backingfile = paste0(PATH, 'G_test'))
  
  Gstd_test<-big_aggregate(G_test, function(X){ 
    as_FBM(scale(X, scale=T, center=T), backingfile = paste0(PATH, 'Gstd_test'))
  }) 
  
  A_test<-big_aggregate(Gstd_test[[1]], function(X){  
    tcrossprod(beta, X)
  })
  
  A_test<-t(A_test[[1]])

  E_test <- rnorm(test_prop, 0, sqrt(1-heritability))

  P_test <- A_test+E_test

  df_test<-as_FBM(matrix(c(A_test, E_test, P_test), nrow=test_prop, ncol=3), backingfile=paste0(PATH, 'df_test'))

  P01_test <- (df_test[][,3]>t)*1

  PRS_test<-big_aggregate(Gstd_test[[1]], function(X){
    tcrossprod(beta_hat[[1]], X)
  })

  AUC_test<-big_aggregate(PRS_test[[1]], function(X){
    AUC(X, P01_test)
  })

  liste_AUC_test[j]<-AUC_test[[1]]
  
  print(AUC_test)
  print(nb_cas_total)
  print(nb_control_total)
  
}
z = dnorm(qnorm(1-k), 0, 1)
i = z/k
v = -i*k/(1-k)
AUCmax = pnorm(((i-v)*heritability)/sqrt(heritability*(1-heritability*i*(i-t)+(1-heritability*v*(v-t)))))
print(AUCmax)


mean_AUC_test<-mean(liste_AUC_test)
mean_cas_total = mean(liste_cas_total)
mean_control_total = mean(liste_control_total)
print(mean_AUC_test)
sd_AUC_test<-sd(liste_AUC_test)
var_AUC_test<-var(liste_AUC_test)
df<-data.frame(mean_AUC_test, sd_AUC_test, var_AUC_test, nb_individuals, nb_snps, heritability, k, mean_cas_total, mean_control_total)

write.table(df, file=paste0(OUTPUT_PATH,'v12_AUCs.txt'), quote=TRUE, append = TRUE, dec=".", row.names=FALSE, col.names=FALSE, sep ="\t", qmethod = c("escape"))

rm(AUC_test, t)
gc()
print('removal done')
unlink(paste0(paste0(PATH, '*'), c(".bk", ".rds")))
