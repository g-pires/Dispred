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


parser=argparse.ArgumentParser(prog='Liability Threshold Model')
parser.add_argument('--n', help='nb_individuals', type=int)
parser.add_argument('--s', help='nb_SNPs', type=int)
parser.add_argument('--h2', help='heritability', type=float)
parser.add_argument('--k', help='prevalence', type=float)
parser.add_argument('--l', help='sibling recurrence risk (lambda_s)', type=float)
args=parser.parse_args()

"""
    This module is building a Liability Threshold Model.

Usage::

    python3 LTM.py --n INT --s INT --h2 FLOAT --k FLOAT --l FLOAT

Input: - n : number of individuals
       - s : number of SNPs
       - h2 : heritability
	   - k : disease prevalence

Output: - boxplot showing the distribution of Y = 0 or 1 in function of the PRS
		- distribution plot of Y predicted
		- ROC curve plot
		- AUC and AUCmax

"""


nb_individuals = args.n

nb_SNPs = args.s

heritability = args.h2

k = args.k #prevalence

lambda_s = args.l


box_plot = 'Figures/boxplot_'+str(nb_individuals)+'_'+str(nb_SNPs)+'_'+str(heritability)+'_'+str(k)+'_'+str(lambda_s)+'.png'

distribution_plot = 'Figures/distrib_'+str(nb_individuals)+'_'+str(nb_SNPs)+'_'+str(heritability)+'_'+str(k)+'_'+str(lambda_s)+'.png'

roc_plot = 'Figures/roc_'+str(nb_individuals)+'_'+str(nb_SNPs)+'_'+str(heritability)+'_'+str(k)+'_'+str(lambda_s)+'.png'


def generate_probas():
	'''
		Function used to generate the p in the binomial distribution. The p follows an uniform distribution
	'''
	p = np.random.uniform(0.05, 0.95, nb_SNPs)
	return p


def generate_sample(n, p):
	'''
		Function used to generate sample folowing a binomial distribution.
		params: n: the number of events.
		params: p: the probability for each event.
	'''
	#np.random.seed(12)
	s = np.random.binomial(n, p, [nb_individuals, nb_SNPs])
	table = pd.DataFrame(s, columns=['SNP'+str(i) for i in range(0,nb_SNPs)], index=['Indiv'+str(j) for j in range(0,nb_individuals)])
	return table


def compute_threshold(k):
	'''
		Function used to compute the threshold from the prevalence (K)
		params: k: prevalence of the disease.
	'''
	t = norm.ppf(1-k) 
	return t


def compute_threshold1(k, lambda_s):
	t1 = norm.ppf(1-lambda_s*k)
	return t1


def transform_scaler(table):
	'''
		Function used to standardize the value from a table. Standardization of the genotypes with 0 mean and unit variance.
		params: table: the table with the genotypes (indiv in rows and SNPs in columns).
	'''
	scaler = preprocessing.StandardScaler().fit(table)
	data_transformed = scaler.transform(table)
	table_trans = pd.DataFrame(data_transformed, columns=['SNP'+str(i) for i in range(0,nb_SNPs)], index=['Indiv'+str(j) for j in range(0,nb_individuals)])
	return table_trans


def generate_betas(mu, sigma):
	'''
		Function used to generate the betas from a random normal distribution.
		params: mu: the mean.
		params: sigma: the variance.	
	'''
	#np.random.seed(12)
	b = np.random.normal(mu, sigma, [nb_SNPs])
	return b


def Y_pred(table, b, heritability):
	'''
		Function used to compute the Y*. Y*=beta*Gi+noise. Y_type = vector. Size = number of individuals.
		params: table: the table with the genotypes (Gi).
		params: b: beta generated randomly.
	'''
	liste = []
	for columns, row in table.iterrows():
		Y = np.dot(b, row.values) 
		liste.append(Y)
	bruit = np.random.normal(0, sqrt(1-heritability), nb_individuals)
	liste_final = liste + bruit
	return liste_final


def infer_status(Y_predicted, t):
	'''
		Function used to compute Y from Y*. If Y* > t then Y=1 else Y=0.
		params: Y_predicted: Y* from Y*=beta*Gi+noise.
		params: t: threshold computed from the prevalence.
	'''
	Y = []
	for i in Y_predicted:
		if i > t:
			Y.append(1)
		else:
			Y.append(0)
	return Y


def logistic_regression(X, Y):
	'''
		Function used to do a logistic regression on two tables.
		params: X: the genotype table, (indiv, SNPs)
		params: Y: the table with the Y* for each SNPs (0 or 1).
	'''
	liste = []
	params = []
	for i in range(0, nb_SNPs):
		for columns, row in X.iterrows():
			liste.append(row['SNP'+str(i)])
		model = sm.Logit(np.array(Y).astype(float), np.array(liste).astype(float))
		res = model.fit(disp=0)
		liste = []
		params.append(float(res.params))
	return params


def compute_odds_ratio(beta_i):
	'''
		Function used to compute the odds_ratios from the coefficient of the logistic regression.
		params: beta_i: coefficient from the Logistic Regression of Y on the genotypes.
	'''
	odds_ratio=[]
	for i in beta_i:
		odds_ratio.append(log(exp(i)))
	return odds_ratio


def compute_PRS(table, OR):
	'''
		Function used to compute the PRS from the dot product of the odds_ratio and the genotypes table (standardized).
		params: table: standardized genotypes table.
		params: OR: odds_ratio.
	'''
	liste=[]
	for columns, row in table.iterrows():
		liste.append(np.dot(OR, row.values))
	return liste



def compute_AUC_2(PRS, Y):
	'''
		Function used to compute the AUC. From the equation of Wray and al. : AUC = 1/n_d_prime(r_d-n_d/2-1/2). n_d_prime : non-diseased indivs, n_d : diseased indivs, 
		r_d : mean rank of diseased indivs.
		params: PRS: the table with the PRS.
		params: Y: the table with the case controls status (1 or 0).

	'''
	df=pd.DataFrame(PRS, columns=['PRS'])
	df['Y']=Y
	new_df=df.sort_values(by=['PRS']).reset_index()
	liste=new_df.index[new_df['Y'].values==1] + 1
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


def compute_ROC(PRS, Y):
	'''
		Function used to compute the ROC curve.
		params: PRS: Polygenic Risk Score table of size nb_individuals.
		params: Y: Disease status (0 or 1) of case controls. Size nb_individuals.
	'''
	fpr, tpr, thresholds = metrics.roc_curve(np.array(Y), np.array(PRS))
	auc = metrics.roc_auc_score(Y, PRS)
	plt.plot(fpr,tpr,label="data 1, auc="+str(auc))
	plt.legend(loc=4)
	plt.savefig(roc_plot)
	plt.show()
			
		
def create_boxplot(dataframe):
	'''
		Function used to create a boxplot from a dataframe.
		params: dataframe: data on which the boxplot is done.
	'''
	boxplot = dataframe.boxplot(column=['PRS'], by='Y')
	plt.savefig(box_plot)
	plt.show()


def create_distribution_plot(table, mu, sigma, t):
	'''
		Function used to create a distribution plot for a normal distribution.
		params: table: the data that need to be plotted.
		params: mu: the mean (usually 0).
		params: sigma: the variance (usually 1).
		params: t: the threshold for a liability threshold model.
	'''
	count, bins, ignored = plt.hist(table, density=True)
	plt.plot(bins, 1/(sigma * np.sqrt(2 * np.pi)) *
                np.exp( - (bins - mu)**2 / (2 * sigma**2) ),
          linewidth=2, color='r')
	plt.axvline(x=t,linewidth=1, color='k')
	plt.savefig(distribution_plot)
	plt.show()



def main():
	p = generate_probas()

	table = generate_sample(2, p)

	t = compute_threshold(k)

	t1 = compute_threshold1(k, lambda_s)

	table_trans = transform_scaler(table)

	b = generate_betas(0, sqrt(heritability/nb_SNPs))

	Y_predicted = Y_pred(table_trans, b, heritability)

	Y = infer_status(Y_predicted, t)
	
	beta_i = logistic_regression(table_trans, Y)
	
	odds_ratio = compute_odds_ratio(beta_i)

	PRS = compute_PRS(table_trans, odds_ratio)

	df = pd.DataFrame(PRS, columns=['PRS'])

	df['Y'] = Y
	
	AUC_2 = compute_AUC_2(PRS, Y)
	
	print('AUC_2 : '+str(AUC_2))

	AUCmax = compute_AUC_3(t, k, heritability)
	
	print('AUCmax : '+str(AUCmax))

	compute_ROC(PRS, Y)

	create_boxplot(df)

	create_distribution_plot(np.array(Y_predicted), mu=0, sigma=1, t=t)

if __name__ == "__main__":
	main()
