suppressMessages(library('Xmisc'))

parser <- ArgumentParser$new()
parser$add_argument('--n',type='numeric',  help='Number of individuals')
parser$add_argument('--s',type='numeric',  help='Number of snps')
parser$add_argument('--h2',type='numeric', help='Heritability')
parser$add_argument('--k', type='numeric', help='Disease prevalence')



#############################################################
#	This code is building a Threshold Liability Model.
#	
#	Usage : 
#			Rscript R_LTM2.R --n NUM --s NUM --h2 NUM --k NUM
#
#	Input : 
#			- n : number of individuals
#			- s : number of SNPs
#			- h2 : heritability
#			- k : disease prevalence
#	
#	Output :
#			- AUCmax and mean_AUC over 100 simulations  
############################################################


nb_individuals = n

nb_snps = s

heritability = h2

k = k

liste_AUC = c()

for (i in 0:100){
  mafs<-runif(nb_snps, 0.05,0.95)
  
  r <- rbinom(nb_individuals*nb_snps, 2, mafs)
  
  G <- matrix(r, nrow = nb_individuals, ncol=nb_snps, byrow=T)
  
  t<-qnorm(1-k)
  
  Gstd<-scale(G, center=T, scale=T)
  
  b<-rnorm(nb_snps, 0, sqrt(heritability/nb_snps))
  
  A<-b%*%t(Gstd)
  
  A<-t(A)
  
  E<-rnorm(nb_individuals, 0, sqrt(1-heritability))
  
  P<-A+E
  
  P01<-(P>t)*1
  
  df<-data.frame(A, P01)
  
  ordered_df<-df[order(A),]
  
  liste<-(ordered_df$P01==1)*1
  
  mean_rank<-mean(which(liste==1))
  
  nb_diseased = sum(liste)
  
  nb_not_diseased = nb_individuals-nb_diseased
  
  AUC<-(1/nb_not_diseased)*(mean_rank-(nb_diseased/2)-(1/2))
  
  liste_AUC[i]<-AUC
  
  print(AUC)
  print(i)
} 


print(mean(liste_AUC))





