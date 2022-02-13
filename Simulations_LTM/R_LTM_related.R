suppressMessages(library('Xmisc'))


parser <- ArgumentParser$new()
parser$add_argument('--n',type='numeric',  help='Number of individuals')
parser$add_argument('--h2',type='numeric', help='Heritability')
parser$add_argument('--k', type='numeric', help='Disease prevalence')



######################################################
#	This code is building a Threshold Liability Model.
#	
#	Usage : 
#			Rscript R_LTM2.R --n NUM --h2 NUM --k NUM
#
#	Input : 
#			- n : number of individuals
#			- h2 : heritability
#			- k : disease prevalence
#	
#	Output :
#			- AUCmax and mean_AUC over 100 simulations  
######################################################


nb_individuals = n

heritability = h2

k = k


for (i in 0:100){
  P1<-rnorm(nb_individuals, 0, sqrt(heritability))
  
  P2<-rnorm(nb_individuals, 0, sqrt(heritability))
  
  mendel<-rnorm(nb_individuals, 0, sqrt(1/2*heritability))
  
  E1<-1/2*P1+1/2*P2+mendel
  
  E2<-1/2*P1+1/2*P2+mendel
  
  t = qnorm(1-k)
  
  E<-rnorm(nb_individuals, 0, sqrt(1-heritability))
  
  A<-c(P1, P2, E1, E2)
  
  P<-A+E
  
  P01<-(P>t)*1
  
  df<-data.frame(A, P01)
  
  ordered_df<-df[order(A),]
  
  liste<-(ordered_df$P01==1)*1
  
  mean_rank<-mean(which(liste==1))
  
  nb_diseased = sum(liste)
  
  nb_not_diseased = nb_individuals*4-nb_diseased
  
  AUC<-(1/nb_not_diseased)*(mean_rank-(nb_diseased/2)-(1/2))
  
  liste_AUC<-AUC
  
  print(AUC)
  print(i)
  
}

print(mean(liste_AUC))










