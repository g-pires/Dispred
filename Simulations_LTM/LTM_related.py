#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from math import *
import pandas as pd
from sklearn import preprocessing
from sklearn.linear_model import LogisticRegression
import statsmodels.api as sm
import argparse
from statistics import *
from sklearn import metrics
import seaborn as sns
import time

parser=argparse.ArgumentParser(prog='Liability Threshold Model')
parser.add_argument('--n', help='nb_individuals', type=int)
parser.add_argument('--h2', help='heritability', type=float)
parser.add_argument('--k', help='prevalence', type=float)
parser.add_argument('--l', help='sibling recurrence risk (lambda_s)', type=float)
args=parser.parse_args()

"""
    This module is building a Liability Threshold Model.

Usage::

    python3 LTM.py --n INT --h2 FLOAT --k FLOAT --l FLOAT --d STR --r STR

Input: - n : number of individuals
       - h2 : heritability
	   - k : disease prevalence
	   - l : recurrence sibling risk

Output: - distribution plot of A, P01 = 0 and P01 = 1
		- ROC curve plot
		- AUC and AUCmax
		- mean AUC for 100 simulations

"""

nb_individuals = args.n

heritability = args.h2

k = args.k

lambda_s = args.l

distribution_plot = 'Figures/distrib_'+str(nb_individuals)+'_'+str(heritability)+'_'+str(k)+'_'+str(lambda_s)+'.png'

roc_plot = 'Figures/roc_'+str(nb_individuals)+'_'+str(heritability)+'_'+str(k)+'_'+str(lambda_s)+'.png'


def generate_families(mu, heritability):
	'''
		Function used to generate nuclear families, 2 parents (P1 and P2), 2 children (E1 and E2). 
		Return dataframe A from the equation P=A+E. P the liability, A the genetic variation and E the environmental factors.
		params: mu: mean of the distribution.
		params: heritability: proportion of variance explained by genetic variation.
	'''
	P1 = np.random.normal(mu, sqrt(heritability), int(nb_individuals/4))
	P2 = np.random.normal(mu, sqrt(heritability), int(nb_individuals/4))
	mendel = np.random.normal(mu, sqrt(1/2*heritability), int(nb_individuals/4))
	E1 = 1/2*P1+1/2*P2+mendel
	E2 = 1/2*P1+1/2*P2+mendel
	liste = list(map(list, zip(P1, P2, E1, E2)))
	liste_final=sum(liste, [])
	df = pd.DataFrame(liste_final, columns=['A'])
	return df


def compute_threshold(k, lambda_s):
	'''
		Function used to compute the threshold from the prevalence (K)
		params: k: prevalence of the disease.
		params: lambda_s: sibling recurrence risk.
	'''
	t = norm.ppf(1-k)
	#t = norm.ppf(1-lambda_s*k)
	return t


def compute_P(table):
	'''
		Function used to compute the E and then the P (liability) from the equation P=A+E.
		params: table: table with A values.
	'''
	E = np.random.normal(0, sqrt(1-heritability), nb_individuals) #environmental effect
	table['E']=E
	table['P']=table['A']+table['E']
	return table


def infer_status(table, t):
	'''
		Function used to compute P01 from P. If P > t then P01=1 else P01=0.
		params: t: threshold computed from the prevalence.
	'''
	table['P01']=table['P'].gt(t)
	table['P01']=table['P01']*1
	return table


def compute_AUC_2(table):
	'''
		Function used to compute the AUC from the equation 2 of Wray paper : AUC = 1/n_d_prime*(r_d-n_d/2-1/2).
		n_d_prime : non diseased individuals, n_d : diseased individuals, r_d : mean rank of diseased individuals.
		params: table: the table with the A and P01 status.
	'''
	df=pd.DataFrame(table['A'], columns=['A'])
	df['P01']=table['P01']
	new_df=df.sort_values(by=['A']).reset_index()
	liste=new_df.index[new_df['P01']==1] + 1
	liste_final=liste.tolist()
	mean_rank=mean(liste_final)
	nb_diseased=len(liste_final)
	nb_not_diseased = nb_individuals - nb_diseased
	AUC = (1/nb_not_diseased)*(mean_rank-(nb_diseased/2)-(1/2))
	return AUC


def compute_AUC_3(t, k, heritability):
	'''
		Function used to compute the AUC from the equation 3 of Wray paper : AUCmax = phi(((i-v)*heritability)/(sqrt(heritability*(1-heritability*i*(i-t)+(1-heritability*v*(v-t)))))
		params: t: the thresold computed by the function compute_threshold().
		params: k: the disease prevalence.
		params: heritability: proportion of phenotypic variation explaine by genetic variations. 
	'''
	z = norm.pdf(norm.ppf(1-k), 0, 1)
	i = z / k
	v = -i*k/(1-k)
	AUCmax = norm.cdf(((i-v)*heritability)/(sqrt(heritability*(1-heritability*i*(i-t)+(1-heritability*v*(v-t))))))
	return AUCmax


def compute_ROC(A, P01):
	'''
		Function used to compute the ROC curve.
		params: PRS: Polygenic Risk Score table of size nb_individuals.
		params: Y: Disease status (0 or 1) of case controls. Size nb_individuals.
	'''
	fpr, tpr, thresholds = metrics.roc_curve(np.array(P01), np.array(A))
	auc = metrics.roc_auc_score(P01, A)
	plt.plot(fpr,tpr,label="data 1, auc="+str(auc))
	plt.legend(loc=4)
	plt.savefig(roc_plot)
	plt.show()


def create_distribution_plot(table, mu, sigma, t):
	'''
		Function used to create a distribution plot for a normal distribution.
		params: table: the data that need to be plotted.
		params: mu: the mean (usually 0).
		params: sigma: the variance (usually 1).
		params: t: the threshold for a liability threshold model.
	'''
	sns.kdeplot(table.loc[(table['P01']==1), 
            'A'], color='r', shade=True, Label='1', bw=0.5, legend=True) 
  
	sns.kdeplot(table.loc[(table['P01']==0),  
            'A'], color='b', shade=True, Label='0', bw=0.5, legend=True) 
	plt.xlabel('A') 
	plt.ylabel('Probability Density') 
	plt.savefig(distribution_plot)
	plt.show()


def run():
	'''
		Function used to run the model.
	'''
	f=generate_families(0, heritability)
	
	t=compute_threshold(k, lambda_s)
	
	t_final = compute_P(f)	

	t_final2=infer_status(t_final, t)

	AUC = compute_AUC_2(t_final2)

	AUCmax = compute_AUC_3(t, k, heritability)
	
	print('AUC : '+str(AUC))
	
	print('AUCmax : '+str(AUCmax))
	
	df = pd.DataFrame(t_final2['A'], columns=['A'])

	df['P01'] = t_final2['P01']

	A = list(t_final2['A'].values)
	
	P01 = list(t_final2['P01'].values)

	#compute_ROC(A, P01)

	#create_distribution_plot(df, mu=0, sigma=1, t=t)

	return AUC


def main():
	'''
		Function used to run the model and simulations.
	'''
	liste=[]
	for i in range(0, 100):
		run()
		AUC=run()
		print('Iteration : '+str(i))
		liste.append(AUC)
		print(time.ctime())
	print('Mean AUC : '+str(mean(liste)))

if __name__=="__main__":
	main()
