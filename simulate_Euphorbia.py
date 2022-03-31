#!/usr/bin/python3

## in order to use this code you have to have ms installed on your computer
## ms can be freely downloaded from:
## http://home.uchicago.edu/rhudson1/source/mksamples.html

import random
import os
import math
import shlex, subprocess
import numpy as np

def ms2nparray(xfile):
	g = list(xfile)
	k = [idx for idx,i in enumerate(g) if len(i) > 0 and i.startswith(b'//')]
	f = []
	for i in k:
		L = g[i+4:i+nDNANsam+4]
		q = []
		for i in L:
			i = i = [int(j) for j in list(i.decode('utf-8'))]
			i = np.array(i)
			q.append(i)
		q = np.array(q)
		q = q.astype("int8")
		f.append(np.array(q))
	return f

### variable declarations

#define the number of simulations
Priorsize = 10000

## Sample size of Fuerteventura+Lanzarote (15 individuos)
nFL = 15*2
## Sample size of Gran Canaria (19 individuos)
nGC = 19*2
## Sample size of Tenerife (Oriental)  (21 individuos)
nTor = 21*2
## Sample size of Tenerife (Occidental) (7 individuos)
nToc = 7*2
## Sample size of La Gomera (29 individuos)
nLG = 29*2
## Sample size of La Palma (5 individuos)
nLP = 5*2
## Sample size of El Hierro (5 individuos)
nEH = 5*2
## Sample size of Continent (15 individuos)
nC = 15*2

## Sample sample size of all lineages.
nDNANsam = nFL + nGC + nTor + nToc + nLG + nLP + nEH +nC

Mod_SSCont = []
Mod_SSBackCol = []
Mod_SSHBackCol = []
Mod_SSHclim = []
Mod_EastWest = []

## create a file to store parameters and one to store the models
parSSCont = open("parSSCont.txt","w")
parSSBackCol = open("parSSBackCol.txt","w")
parSSHBackCol = open("parSSHBackCol.txt","w")
parSSHclim = open("parSSHclim.txt","w")
parEastWest = open("parEastWest.txt","w")

### MODEL SSCo
for i in range(Priorsize):

	### Define parameters
	## number of years per generation
	genlen = random.uniform(5,6)
	##mutation rate
	mutrate = 1E-8*genlen
	## Theta values from 1 to 5. 
	Theta = random.uniform(1,5)

	##Calculate Ne for each simulation according to the given Theta value
	Ne = Theta/(4*mutrate*25)

	##relative values of Ne in each island.
	#Values are proportional to the Ne of FU_LZ
	NGC = random.uniform(0.2,2)
	NTor = random.uniform(0.2,2)
	NToc = random.uniform(0.2,2)
	NLG = random.uniform(0.2,2)
	NLP = random.uniform(0.2,2)
	NEH = random.uniform(0.2,2)
	NC = random.uniform(0.2,2)
	
	## migration prior set to 0 in this model.
	m=0

	## divergence time prior in years, following uniform distributions.
	#T8 is not used in this model (assign default value).
	Tm=0
	Tb=random.uniform(14000,1250000)
	T7=random.uniform(Tb,3500000)
	T6=random.uniform(Tb,T7)
	T5=random.uniform(Tb,T6)
	T4=random.uniform(Tb,T5)
	T3=random.uniform(Tb,T4)
	T2=random.uniform(Tb,T3)
	T1=random.uniform(Tb,T2)

	## Transform to coalescent units
	coalT7=T7/(genlen*4.0*Ne)
	coalT6=T6/(genlen*4.0*Ne)
	coalT5=T5/(genlen*4.0*Ne)
	coalT4=T4/(genlen*4.0*Ne)
	coalT3=T3/(genlen*4.0*Ne)
	coalT2=T2/(genlen*4.0*Ne)
	coalT1=T1/(genlen*4.0*Ne)
	coalTb=Tb/(genlen*4.0*Ne)

	##Size of founded populations
	FoundedSizeRatio = random.uniform(0.0, 0.05)
	##Calculates the growth rate
	GrowthT1 = -(1/T6)*math.log(1/(1/FoundedSizeRatio))
	GrowthT2 = -(1/T5)*math.log((1/NGC)/(1/FoundedSizeRatio))
	GrowthT3 = -(1/T4)*math.log((1/NTor)/(1/FoundedSizeRatio))
	GrowthT4 = -(1/T3)*math.log((1/NToc)/(1/FoundedSizeRatio))
	GrowthT5 = -(1/T2)*math.log((1/NLP)/(1/FoundedSizeRatio))
	GrowthT6 = -(1/T1)*math.log((1/NLG)/(1/FoundedSizeRatio))
	GrowthT7 = -(1/T1)*math.log((1/NEH)/(1/FoundedSizeRatio))	
	GrowthT8 = -(1/T7-Tb)*math.log((1/NC)/(1/FoundedSizeRatio))
	
	## simulate SNPs
	com=subprocess.Popen("./ms %d 1000 -s 1 -t %f -I 8 %d %d %d %d %d %d %d %d -n 2 %f -n 3 %f -n 4 %f -n 5 %f -n 6 %f -n 7 %f -n 8 %f -g 1 %f -g 2 %f -g 3 %f -g 4 %f -g 5 %f -g 6 %f -g 7 %f -g 8 %f -eg %f 8 %f -en %f 7 %f -en %f 6 %f -ej %f 7 6 -en %f 5 %f -ej %f 6 5 -en %f 4 %f -ej %f 5 4 -en %f 3 %f -ej %f 4 3 -en %f 2 %f -ej %f 3 2 -en %f 1 %f -ej %f 2 1 -en %f 8 %f -ej %f 1 8" % (nDNANsam, Theta, nFL, nGC, nTor, nToc, nLP, nLG, nEH, nC, NGC, NTor, NToc, NLP, NLG, NEH, NC, GrowthT1, GrowthT2, GrowthT3, GrowthT4, GrowthT5, GrowthT6, GrowthT7, GrowthT8, coalTb, FoundedSizeRatio, coalT1, FoundedSizeRatio, coalT1, FoundedSizeRatio, coalT1, coalT2, FoundedSizeRatio, coalT2, coalT3, FoundedSizeRatio, coalT3, coalT4, FoundedSizeRatio, coalT4, coalT5, FoundedSizeRatio, coalT5, coalT6, FoundedSizeRatio, coalT6, coalT7, FoundedSizeRatio, coalT7), shell=True, stdout=subprocess.PIPE).stdout
	output = com.read().splitlines()
	Mod_SSCont.append(np.array(ms2nparray(output)).swapaxes(0,1).reshape(nDNANsam,-1).T)
	
	## save parameter values and models
	parSSCont.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (Theta, Tm, Tb, T1, T2, T3, T4, T5, T6, T7, Ne, NGC, NTor, NToc, NLP, NLG, NEH, NC, m, FoundedSizeRatio))
	print("Completed %d %% of Model 1 simulations" % (float(i)/Priorsize*100))

Mod_SSCont=np.array(Mod_SSCont)
np.save('Mod_SSCont.npy', Mod_SSCont)
del(Mod_SSCont)

### MODEL BackCo
for i in range(Priorsize):

	### Define parameters
	## number of years per generation
	genlen = random.uniform(5,6)
	##mutation rate
	mutrate = 1E-8*genlen
	## Theta values from 1 to 5. 
	Theta = random.uniform(1,5)

	##Calculate Ne for each simulation according to the given Theta value
	Ne = Theta/(4*mutrate*25)

	##relative values of Ne in each island.
	#Values are proportional to the Ne of FU_LZ
	NGC = random.uniform(0.2,2)
	NTor = random.uniform(0.2,2)
	NToc = random.uniform(0.2,2)
	NLG = random.uniform(0.2,2)
	NLP = random.uniform(0.2,2)
	NEH = random.uniform(0.2,2)
	NC = random.uniform(0.2,2)
	
	## migration prior set to 0 in this model.
	m=0

	## divergence time prior in years, following uniform distributions.
	#First, sample the divergence time between continent and FL (most recent)
	#T7 and T8 are not used in this model (assign default value).
	T7=0
	Tm=0
	Tb=random.uniform(14000,1250000)
	T6=random.uniform(Tb,3500000)
	T5=random.uniform(Tb,T6)
	T4=random.uniform(Tb,T5)
	T3=random.uniform(Tb,T4)
	T2=random.uniform(Tb,T3)
	T1=random.uniform(Tb,T2)

	## Transform to coalescent units
	coalT6=T6/(genlen*4.0*Ne)
	coalT5=T5/(genlen*4.0*Ne)
	coalT4=T4/(genlen*4.0*Ne)
	coalT3=T3/(genlen*4.0*Ne)
	coalT2=T2/(genlen*4.0*Ne)
	coalT1=T1/(genlen*4.0*Ne)
	coalTb=Tb/(genlen*4.0*Ne)

	##Size of founded populations
	FoundedSizeRatio = random.uniform(0.0, 0.05)
	##Calculates the growth rate
	GrowthT2 = -(1/T5)*math.log((1/NGC)/(1/FoundedSizeRatio))
	GrowthT3 = -(1/T4)*math.log((1/NTor)/(1/FoundedSizeRatio))
	GrowthT4 = -(1/T3)*math.log((1/NToc)/(1/FoundedSizeRatio))
	GrowthT5 = -(1/T2)*math.log((1/NLP)/(1/FoundedSizeRatio))
	GrowthT6 = -(1/T1)*math.log((1/NLG)/(1/FoundedSizeRatio))
	GrowthT7 = -(1/T1)*math.log((1/NEH)/(1/FoundedSizeRatio))	
	GrowthT8 = -(1/Tb)*math.log((1/NC)/(1/FoundedSizeRatio))
		
	## simulate SNPs
	com=subprocess.Popen("./ms %d 1000 -s 1 -t %f -I 8 %d %d %d %d %d %d %d %d -n 2 %f -n 3 %f -n 4 %f -n 5 %f -n 6 %f -n 7 %f -n 8 %f -g 2 %f -g 3 %f -g 4 %f -g 5 %f -g 6 %f -g 7 %f -g 8 %f -en %f 8 %f -ej %f 8 1 -en %f 7 %f -en %f 6 %f -ej %f 7 6 -en %f 5 %f -ej %f 6 5 -en %f 4 %f -ej %f 5 4 -en %f 3 %f -ej %f 4 3 -en %f 2 %f -ej %f 3 2 -en %f 1 %f -ej %f 2 1" % (nDNANsam, Theta, nFL, nGC, nTor, nToc, nLP, nLG, nEH, nC, NGC, NTor, NToc, NLP, NLG, NEH, NC, GrowthT2, GrowthT3, GrowthT4, GrowthT5, GrowthT6, GrowthT7, GrowthT8, coalTb, FoundedSizeRatio, coalTb, coalT1, FoundedSizeRatio, coalT1, FoundedSizeRatio, coalT1, coalT2, FoundedSizeRatio, coalT2, coalT3, FoundedSizeRatio, coalT3, coalT4, FoundedSizeRatio, coalT4, coalT5, FoundedSizeRatio, coalT5, coalT6, FoundedSizeRatio, coalT6), shell=True, stdout=subprocess.PIPE).stdout
	output = com.read().splitlines()
	Mod_SSBackCol.append(np.array(ms2nparray(output)).swapaxes(0,1).reshape(nDNANsam,-1).T)
	
	## save parameter values and models
	parSSBackCol.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (Theta, Tm, Tb, T1, T2, T3, T4, T5, T6, T7, Ne, NGC, NTor, NToc, NLP, NLG, NEH, NC, m, FoundedSizeRatio))
	print("Completed %d %% of Model 2 simulations" % (float(i)/Priorsize*100))

Mod_SSBackCol=np.array(Mod_SSBackCol)
np.save('Mod_SSBackCol.npy', Mod_SSBackCol)
del(Mod_SSBackCol)

### MODEL BackCo+m
for i in range(Priorsize):

	### Define parameters
	## number of years per generation
	genlen = random.uniform(5,6)
	##mutation rate
	mutrate = 1E-8*genlen
	## Theta values from 1 to 5. 
	Theta = random.uniform(1,5)

	##Calculate Ne for each simulation according to the given Theta value
	Ne = Theta/(4*mutrate*25)

	##relative values of Ne in each island.
	#Values are proportional to the Ne of FU_LZ
	NGC = random.uniform(0.2,2)
	NTor = random.uniform(0.2,2)
	NToc = random.uniform(0.2,2)
	NLG = random.uniform(0.2,2)
	NLP = random.uniform(0.2,2)
	NEH = random.uniform(0.2,2)
	NC = random.uniform(0.2,2)
	
	## migration prior from 0 to 10 in this model.
	m=random.uniform(0,10)

	## divergence time prior in years, following uniform distributions.
	#First, sample the divergence time between continent and FL (most recent)
	#T7 and T8 are not used in this model (assign default value).
	Tb=random.uniform(14000,1250000)
	Tm=random.uniform(13000,Tb)
	T7=0
	T6=random.uniform(Tb,3500000)
	T5=random.uniform(Tb,T6)
	T4=random.uniform(Tb,T5)
	T3=random.uniform(Tb,T4)
	T2=random.uniform(Tb,T3)
	T1=random.uniform(Tb,T2)

	## Transform to coalescent units
	coalT6=T6/(genlen*4.0*Ne)
	coalT5=T5/(genlen*4.0*Ne)
	coalT4=T4/(genlen*4.0*Ne)
	coalT3=T3/(genlen*4.0*Ne)
	coalT2=T2/(genlen*4.0*Ne)
	coalT1=T1/(genlen*4.0*Ne)
	coalTb=Tb/(genlen*4.0*Ne)
	coalTm=Tm/(genlen*4.0*Ne)

	##Size of founded populations
	FoundedSizeRatio = random.uniform(0.0, 0.05)
	##Calculates the growth rate
	GrowthT2 = -(1/T5)*math.log((1/NGC)/(1/FoundedSizeRatio))
	GrowthT3 = -(1/T4)*math.log((1/NTor)/(1/FoundedSizeRatio))
	GrowthT4 = -(1/T3)*math.log((1/NToc)/(1/FoundedSizeRatio))
	GrowthT5 = -(1/T2)*math.log((1/NLP)/(1/FoundedSizeRatio))
	GrowthT6 = -(1/T1)*math.log((1/NLG)/(1/FoundedSizeRatio))
	GrowthT7 = -(1/T1)*math.log((1/NEH)/(1/FoundedSizeRatio))	
	GrowthT8 = -(1/Tb)*math.log((1/NC)/(1/FoundedSizeRatio))
		
	## simulate SNPs
	com=subprocess.Popen("./ms %d 1000 -s 1 -t %f -I 8 %d %d %d %d %d %d %d %d -n 2 %f -n 3 %f -n 4 %f -n 5 %f -n 6 %f -n 7 %f -n 8 %f -g 2 %f -g 3 %f -g 4 %f -g 5 %f -g 6 %f -g 7 %f -g 8 %f -em %f 1 2 %f -em %f 2 3 %f -em %f 1 2 0 -em %f 2 3 0 -en %f 8 %f -ej %f 8 1 -en %f 7 %f -en %f 6 %f -ej %f 7 6 -en %f 5 %f -ej %f 6 5 -en %f 4 %f -ej %f 5 4 -en %f 3 %f -ej %f 4 3 -en %f 2 %f -ej %f 3 2 -en %f 1 %f -ej %f 2 1" % (nDNANsam, Theta, nFL, nGC, nTor, nToc, nLP, nLG, nEH, nC, NGC, NTor, NToc, NLP, NLG, NEH, NC, GrowthT2, GrowthT3, GrowthT4, GrowthT5, GrowthT6, GrowthT7, GrowthT8, coalTm, m, coalTm, m, coalTb, coalTb, coalTb, FoundedSizeRatio, coalTb, coalT1, FoundedSizeRatio, coalT1, FoundedSizeRatio, coalT1, coalT2, FoundedSizeRatio, coalT2, coalT3, FoundedSizeRatio, coalT3, coalT4, FoundedSizeRatio, coalT4, coalT5, FoundedSizeRatio, coalT5, coalT6, FoundedSizeRatio, coalT6), shell=True, stdout=subprocess.PIPE).stdout
	output = com.read().splitlines()
	Mod_SSHBackCol.append(np.array(ms2nparray(output)).swapaxes(0,1).reshape(nDNANsam,-1).T)
	
	## save parameter values and models
	parSSHBackCol.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (Theta, Tm, Tb, T1, T2, T3, T4, T5, T6, T7, Ne, NGC, NTor, NToc, NLP, NLG, NEH, NC, m, FoundedSizeRatio))
	print("Completed %d %% of Model 3 simulations" % (float(i)/Priorsize*100))

Mod_SSHBackCol=np.array(Mod_SSHBackCol)
np.save('Mod_SSHBackCol.npy', Mod_SSHBackCol)
del(Mod_SSHBackCol)

### MODEL SSCO+m
for i in range(Priorsize):

	### Define parameters
	## number of years per generation
	genlen = random.uniform(5,6)
	##mutation rate
	mutrate = 1E-8*genlen
	## Theta values from 1 to 5. 
	Theta = random.uniform(1,5)

	##Calculate Ne for each simulation according to the given Theta value
	Ne = Theta/(4*mutrate*25)

	##relative values of Ne in each island.
	#Values are proportional to the Ne of FU_LZ
	NGC = random.uniform(0.2,2)
	NTor = random.uniform(0.2,2)
	NToc = random.uniform(0.2,2)
	NLG = random.uniform(0.2,2)
	NLP = random.uniform(0.2,2)
	NEH = random.uniform(0.2,2)
	NC = random.uniform(0.2,2)
	
	## migration prior from 0 to 10 in this model.
	m=random.uniform(0,10)

	## divergence time prior in years, following uniform distributions.
	#T8 is not used in this model (assign default value).
	Tb=random.uniform(14000,1250000)
	Tm=random.uniform(13000,Tb)
	T7=random.uniform(Tb,3500000)
	T6=random.uniform(Tb,T7)
	T5=random.uniform(Tb,T6)
	T4=random.uniform(Tb,T5)
	T3=random.uniform(Tb,T4)
	T2=random.uniform(Tb,T3)
	T1=random.uniform(Tb,T2)

	## Transform to coalescent units
	coalT7=T7/(genlen*4.0*Ne)
	coalT6=T6/(genlen*4.0*Ne)
	coalT5=T5/(genlen*4.0*Ne)
	coalT4=T4/(genlen*4.0*Ne)
	coalT3=T3/(genlen*4.0*Ne)
	coalT2=T2/(genlen*4.0*Ne)
	coalT1=T1/(genlen*4.0*Ne)
	coalTb=Tb/(genlen*4.0*Ne)
	coalTm=Tm/(genlen*4.0*Ne)

	##Size of founded populations
	FoundedSizeRatio = random.uniform(0.0, 0.05)
	##Calculates the growth rate
	GrowthT1 = -(1/T6)*math.log(1/(1/FoundedSizeRatio))
	GrowthT2 = -(1/T5)*math.log((1/NGC)/(1/FoundedSizeRatio))
	GrowthT3 = -(1/T4)*math.log((1/NTor)/(1/FoundedSizeRatio))
	GrowthT4 = -(1/T3)*math.log((1/NToc)/(1/FoundedSizeRatio))
	GrowthT5 = -(1/T2)*math.log((1/NLP)/(1/FoundedSizeRatio))
	GrowthT6 = -(1/T1)*math.log((1/NLG)/(1/FoundedSizeRatio))
	GrowthT7 = -(1/T1)*math.log((1/NEH)/(1/FoundedSizeRatio))	
	GrowthT8 = -(1/T7-Tb)*math.log((1/NC)/(1/FoundedSizeRatio))
	
	## simulate SNPs
	com=subprocess.Popen("./ms %d 1000 -s 1 -t %f -I 8 %d %d %d %d %d %d %d %d -n 2 %f -n 3 %f -n 4 %f -n 5 %f -n 6 %f -n 7 %f -n 8 %f -g 1 %f -g 2 %f -g 3 %f -g 4 %f -g 5 %f -g 6 %f -g 7 %f -g 8 %f -em %f 8 1 %f -em %f 1 2 %f -em %f 2 3 %f -em %f 8 1 0 -em %f 1 2 0 -em %f 2 3 0 -eg %f 8 %f -en %f 7 %f -en %f 6 %f -ej %f 7 6 -en %f 5 %f -ej %f 6 5 -en %f 4 %f -ej %f 5 4 -en %f 3 %f -ej %f 4 3 -en %f 2 %f -ej %f 3 2 -en %f 1 %f -ej %f 2 1 -en %f 8 %f -ej %f 1 8" % (nDNANsam, Theta, nFL, nGC, nTor, nToc, nLP, nLG, nEH, nC, NGC, NTor, NToc, NLP, NLG, NEH, NC, GrowthT1, GrowthT2, GrowthT3, GrowthT4, GrowthT5, GrowthT6, GrowthT7, GrowthT8, coalTm, m, coalTm, m, coalTm, m, coalTb, coalTb, coalTb, coalTb, FoundedSizeRatio, coalT1, FoundedSizeRatio, coalT1, FoundedSizeRatio, coalT1, coalT2, FoundedSizeRatio, coalT2, coalT3, FoundedSizeRatio, coalT3, coalT4, FoundedSizeRatio, coalT4, coalT5, FoundedSizeRatio, coalT5, coalT6, FoundedSizeRatio, coalT6, coalT7, FoundedSizeRatio, coalT7), shell=True, stdout=subprocess.PIPE).stdout
	output = com.read().splitlines()
	Mod_SSHclim.append(np.array(ms2nparray(output)).swapaxes(0,1).reshape(nDNANsam,-1).T)
	
	## save parameter values and models
	parSSHclim.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (Theta, Tm, Tb, T1, T2, T3, T4, T5, T6, T7, Ne, NGC, NTor, NToc, NLP, NLG, NEH, NC, m, FoundedSizeRatio))
	print("Completed %d %% of Model 4 simulations" % (float(i)/Priorsize*100))

Mod_SSHclim=np.array(Mod_SSHclim)
np.save('Mod_SSHclim.npy', Mod_SSHclim)
del(Mod_SSHclim)

### MODEL CIH
for i in range(Priorsize):

	### Define parameters
	## number of years per generation
	genlen = random.uniform(5,6)
	##mutation rate
	mutrate = 1E-8*genlen
	## Theta values from 1 to 5. 
	Theta = random.uniform(1,5)

	##Calculate Ne for each simulation according to the given Theta value
	Ne = Theta/(4*mutrate*25)

	##relative values of Ne in each island.
	#Values are proportional to the Ne of FU_LZ
	NGC = random.uniform(0.2,2)
	NTor = random.uniform(0.2,2)
	NToc = random.uniform(0.2,2)
	NLG = random.uniform(0.2,2)
	NLP = random.uniform(0.2,2)
	NEH = random.uniform(0.2,2)
	NC = random.uniform(0.2,2)

	## migration prior set to 0 in this model.
	m=0
	
	## divergence time prior in years, following uniform distributions.
	Tb=0
	Tm=0
	T7=random.uniform(1250000,3500000)
	T6=random.uniform(0,T7)
	T5=random.uniform(0,T6)
	T4=random.uniform(0,T5)
	T3=random.uniform(0,T4)
	T2=random.uniform(0,T3)
	T1=random.uniform(0,T2)

	## Transform to coalescent units
	coalT7=T7/(genlen*4.0*Ne)
	coalT6=T6/(genlen*4.0*Ne)
	coalT5=T5/(genlen*4.0*Ne)
	coalT4=T4/(genlen*4.0*Ne)
	coalT3=T3/(genlen*4.0*Ne)
	coalT2=T2/(genlen*4.0*Ne)
	coalT1=T1/(genlen*4.0*Ne)

	##Size of founded populations
	FoundedSizeRatio = random.uniform(0.0, 0.05)
	##Calculates the growth rate
	GrowthT1 = -(1/T1)*math.log(1/(1/FoundedSizeRatio))
	GrowthT2 = -(1/T4)*math.log((1/NGC)/(1/FoundedSizeRatio))
	GrowthT3 = -(1/T6)*math.log((1/NTor)/(1/FoundedSizeRatio))
	GrowthT4 = -(1/T5)*math.log((1/NToc)/(1/FoundedSizeRatio))
	GrowthT5 = -(1/T2)*math.log((1/NLP)/(1/FoundedSizeRatio))
	GrowthT6 = -(1/T3)*math.log((1/NLG)/(1/FoundedSizeRatio))
	GrowthT7 = -(1/T2)*math.log((1/NEH)/(1/FoundedSizeRatio))	
	GrowthT8 = -(1/T1)*math.log((1/NC)/(1/FoundedSizeRatio))	
		
	## simulate SNPs
	com=subprocess.Popen("./ms %d 1000 -s 1 -t %f -I 8 %d %d %d %d %d %d %d %d -n 2 %f -n 3 %f -n 4 %f -n 5 %f -n 6 %f -n 7 %f -n 8 %f -g 1 %f -g 2 %f -g 3 %f -g 4 %f -g 5 %f -g 6 %f -g 7 %f -g 8 %f -en %f 8 %f -en %f 1 %f -ej %f 8 1 -en %f 7 %f -en %f 6 %f -ej %f 7 6 -en %f 5 %f -ej %f 6 5 -en %f 2 %f -ej %f 2 1 -en %f 4 %f -ej %f 5 4 -en %f 3 %f -ej %f 3 1 -ej %f 4 1" % (nDNANsam, Theta, nFL, nGC, nTor, nToc, nLG, nLP, nEH, nC, NGC, NTor, NToc, NLG, NLP, NEH, NC, GrowthT1, GrowthT2, GrowthT3, GrowthT4, GrowthT5, GrowthT6, GrowthT7, GrowthT8, coalT1, FoundedSizeRatio, coalT1, FoundedSizeRatio, coalT1, coalT2, FoundedSizeRatio, coalT2, FoundedSizeRatio, coalT2, coalT3, FoundedSizeRatio, coalT3, coalT4, FoundedSizeRatio, coalT4, coalT5, FoundedSizeRatio, coalT5, coalT6, FoundedSizeRatio, coalT6, coalT7), shell=True, stdout=subprocess.PIPE).stdout
	output = com.read().splitlines()
	Mod_EastWest.append(np.array(ms2nparray(output)).swapaxes(0,1).reshape(nDNANsam,-1).T)
	
	## save parameter values and models
	parEastWest.write("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (Theta, Tm, Tb, T1, T2, T3, T4, T5, T6, T7, Ne, NGC, NTor, NToc, NLP, NLG, NEH, NC, m, FoundedSizeRatio))
	print("Completed %d %% of Model 5 simulations" % (float(i)/Priorsize*100))

Mod_EastWest=np.array(Mod_EastWest)
np.save('Mod_EastWest.npy', Mod_EastWest)
del(Mod_EastWest)
