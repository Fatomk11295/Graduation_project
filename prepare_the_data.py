# -*- coding: utf-8 -*-
"""Prepare_the_data

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1uBnAHJNwqd3sB7D1bhLZxxdTqbJ3xFrs
"""

#Import essential packages
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

"""# Control dataset

GTEX phenotype dataset
"""

df_pheno = pd.read_csv('/content/drive/MyDrive/Bioinformatics Spe /Graduation_project_12Jul/Datasets/GTEX_phenotype.tsv', delimiter='\t')

#Take out samples collected from Bone marrow tissue
df_pheno_BM = df_pheno.loc[df_pheno['_primary_site'] == 'Bone Marrow']

df_pheno_BM.head()

"""GTEX raw counts dataset (gene expression)"""

df_GTEXraw_counts = pd.read_csv('/content/drive/MyDrive/Bioinformatics Spe /Graduation_project_12Jul/Datasets/gtex_gene_expected_count.tsv', delimiter='\t')

df_GTEXraw_counts.head()

#Set Genes as index (column name is sample)
df_GTEXraw_counts = df_GTEXraw_counts.set_index('sample')

print(f'The GTEX raw counts dataset contains {df_GTEXraw_counts.shape[0]} genes and {df_GTEXraw_counts.shape[1]} samples.')

"""Subset the raw counts dataset to samples of BM """

#find samples which are obtained from bone marrow
p = df_pheno_BM['Sample'].to_list()
c = list(df_GTEXraw_counts.columns)
#convert list to set
p_set = set(p)
c_set = set(c)
#find shared id
shared = p_set.intersection(c_set)

len(shared)

#Subset dataframe based on columns (samples from blood tissue)
BM_samples = list(shared)
controlled_df = df_GTEXraw_counts[df_GTEXraw_counts.columns & BM_samples]

print(f'The controlled dataset contains {controlled_df.shape[0]} genes and {controlled_df.shape[1]} samples.')

controlled_df.head()

"""Check the data"""

#Check null values
for val in controlled_df.isnull().sum():
  if val != 0 :
    print(val)

#drop duplicated genes
controlled_df_WithoutDuplicates = controlled_df.drop_duplicates()

controlled_df_WithoutDuplicates.to_csv('Control_data.csv')

"""# Case dataset"""

#Dataset: gene expression RNAseq - HTSeq - FPKM
cases_df = pd.read_csv('/content/drive/MyDrive/Bioinformatics Spe /Graduation_project_12Jul/Datasets/TCGA-LAML.htseq_counts.tsv',delimiter='\t')

print(f'The Cases dataset contains {cases_df.shape[0]} genes and {cases_df.shape[1]-1} samples.')

#Set gene as index (column name Ensembl_ID)
cases_df = cases_df.set_index('Ensembl_ID')

#Check null values
for val in cases_df.isnull().sum():
  if val != 0 :
    print(val)

#check duplicated genes
cases_df_WithoutDuplicates = cases_df.drop_duplicates()

#Check duplicated samples for same patient
pt_idList = []
for col in cases_df_WithoutDuplicates.columns:
  pt_idList.append(col[:-3])

if cases_df_WithoutDuplicates.shape[1] == len(set(pt_idList)) :
  print('There is only one sample for each patient')

cases_df_WithoutDuplicates.to_csv('cases_df_WithoutDuplicates.csv')

"""# Merge the two dataframes"""

#Merge cases and control dataset based on gene ID
df = pd.merge(controlled_df_WithoutDuplicates, cases_df_WithoutDuplicates,
                    left_index=True, right_index=True)

df.head()

print(f'The final dataset contains {df.shape[0]} genes and {df.shape[1]} samples.')

"""# Process the dataframe"""

#Add 1 to all values so when convert back to raw counts (-1) we wont have negative values
df = df.add(1)

def convert_unit(df):
  '''This function convert normalized counts log2(x+1) to raw counts'''
  #The equation is as following:
  # 2 ** value (right) = x + 1 (left)
  '''Step 1 : Find x value (Reads value before normalization)'''
  right_array = pow(df.values,2)
  '''Step 2 : substract all values in right array by 1'''
  raw_array = right_array - 1
  '''Step 3 : return the dataframe with raw counts'''
  rawCounts_df = pd.DataFrame(raw_array, index=df.index,
                        columns=df.columns)
  return rawCounts_df

rawCounts_df = convert_unit(df)

print(f'The dataset with raw counts contains {rawCounts_df.shape[0]} genes and {rawCounts_df.shape[1]} samples.')

#Remove the transcript number our of the gene id
rawCounts_df.index = [gene[:gene.find('.')]  for gene in rawCounts_df.index ]

rawCounts_df.head()

#check duplicated genes
rawCounts_df_WithoutDuplicates = rawCounts_df.drop_duplicates()

print(f'The dataset with raw counts after dropping duplicates contains {rawCounts_df_WithoutDuplicates.shape[0]} genes and {rawCounts_df_WithoutDuplicates.shape[1]} samples.')

rawCounts_df_WithoutDuplicates.to_csv('counts_data_Ready.csv')

"""Visualize the data distribution"""

plt.figure(figsize=(15,15))
rawCounts_df.plot(kind = 'density' ,legend=False, title='Gene expression distribution among all samples')
plt.show()

"""Create the sample file (Contains a label of the sample type)"""

#Transpose the dataframe
df_Transposed = rawCounts_df.T

df_Transposed.head()

#Create a list to label sample as case or control 
kind= []
for x in df_Transposed.index:
  if x[0] == 'K':
    kind.append('Control')
  if x[0] == 'T':
    kind.append('Case')

#Define the new column Type based on kind list
df_Transposed['Type'] = kind

df_Transposed.head()

#Create a file which contains only sample ID and their type
df_Transposed.iloc[:,-1].to_csv('Sample_list.csv')

"""Visualize genes distribution among each group of samples"""

#Group the samples based on their type by mean value
grouped_df = df_Transposed.groupby(by="Type").mean()

#Retranspose the grouped dataframe
grouped_df_Retranspose = grouped_df.T

grouped_df_Retranspose.boxplot(figsize=(12,8),showfliers=False,fontsize=12)  
#showfliers Show the outliers beyond the caps if set true