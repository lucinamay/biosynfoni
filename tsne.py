#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct  1 21:53:10 2023

@author: lucina-may
"""
from sklearn.manifold import TSNE
#from keras.datasets import mnist
from sklearn.datasets import load_iris
import numpy as np
from numpy import reshape
import seaborn as sns
import pandas as pd  

#parameters

tsne = TSNE(n_components=2, verbose=0, perplexity=40, n_iter=300)
                             
fingerprint_file = "../arch/0816_BioSynFoni_fpsfull_3.csv"
data = []
with open (fingerprint_file, 'r') as f:
    for line in f:
         data.append(line.strip().split(','))
                          
data = np.array(data)

#tsne = TSNE(n_components=2, verbose=1, random_state=123) 
tsne = TSNE(n_components=2, verbose=0, perplexity=40, n_iter=300)
z = tsne.fit_transform(data[:10000])

df = pd.DataFrame()
df["comp-1"] = z[:,0]
df["comp-2"] = z[:,1]           


npcs = []
with open('../arch/0914_COCONUT_DB_rdk_npcs.tsv') as cl:
    for line in cl:
        npcs.append(line.strip().split('\t')[0].split(',')[0])
df["npcs"] = npcs
df["npcs"] = df["npcs"].replace('', "No NP-Classifier prediction")             

subset_df = df.sample(n=10000)
print(subset_df.index().tolist)

