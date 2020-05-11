import data as dt
import Modify_data as mdt
from Bio import Entrez
import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn2
import pandas as pd
import os
import glob

def load_df(bfNumber):
  '''Load dataframe from new files with gene names.'''
  os.chdir("/home/carlier/Documents/Stage/Interactome_proteins_ecoli/") # where is the MS data
  path_batch = "MS data csv reportbuilder-20200408T212019Z-001/MS data csv reportbuilder/"
  letter = bfNumber[0]
  number = bfNumber[1]
  full_path = path_batch+"batch "+letter+"/"+letter+number+'.csv'
  df = pd.read_csv(full_path, sep=',', header=0)
  df = df.set_index("Accession_number")
  return df

pd_samples = dt.load_based_screen_samples()
pd_controls = dt.load_based_screen_controls()
used_prot_tuples = dt.good_proteins() 

def used_bnames(used_prot_tuples):
  '''List of batch names used'''
  bnames = []
  for my_tuple in used_prot_tuples:
    bprot = []
    batch_names = dt.get_batches(all_samples,my_tuple[0], my_tuple[1])
    for brep in batch_names:
      bprot.append(brep)
  bnames.append(bprot)
  return bnames

def intersect(my_list):
  intersect = set(my_list[0]['Gene_name'])
  for i in my_list:
    intersect = intersect.intersection(set(i['Gene_name']))
  return intersect

controls = {'LB log':['S1','S2','S3'], 'LB O/N':['O9', 'R4', 'R5'], 'M9 0.2% ac O/N':['R1', 'R2', 'R3']}
prot1 = used_prot_tuples[0]
prot_batches = dt.get_batches(pd_samples, prot1[0], prot1[1])
rep = [load_df(i) for i in prot_batches]
ctr = [load_df(i) for i in controls[prot1[1]]]

#print(len(intersect([rep[0], rep[1]])))

def venn_diagram(data):
  set_array = []
  for rep in data:
    set_array.append(set(rep['Gene_name']))
  if len(data) == 3:
    set_names = ['Rep1', 'Rep2', 'Rep3']
    venn3(set_array, set_names)   # venn3 works for three sets
    plt.show()
  elif len(data) == 2:
    set_names = ['Rep1', 'Rep2']
    venn2(set_array, set_names)   # venn3 works for three sets
    plt.show()
  else: print('error, please change data length')

#print(used_prot_tuples)

def file_result(filename):
  myfile = open(filename, 'w')
  for prot1 in used_prot_tuples:
    prot_batches = dt.get_batches(pd_samples, prot1[0], prot1[1])
    rep = [load_df(i) for i in prot_batches]
    ctr = [load_df(i) for i in controls[prot1[1]]]
    intR = intersect(rep)
    intC = intersect(ctr)
    pattern = '|'.join(intersect(ctr))
#    print(ctr[0].Protein_name.values[ctr[0].Gene_name.str.contains(pattern)])
    myfile.write('(Protein) '+ prot1[0]+ ' (Condition) '+ prot1[1] + ':\n')
    myfile.write('\tLen rep intersect '+ str(len(intR))+ '\n')
    pNamesR = rep[0].Protein_name.values[rep[0].Gene_name.str.contains('|'.join(intR))]
    myfile.write('\tCommon rep proteins :\n')    
    for i in pNamesR:
      myfile.write(i+'\n')
    myfile.write('\tLen ctr intersect '+ str(len(intC))+'\n')
    pNamesC = ctr[0].Protein_name.values[ctr[0].Gene_name.str.contains('|'.join(intC))]
    myfile.write('\tCommon ctr proteins :\n') 
    for i in pNamesC:
      myfile.write(i+'\n')   
    myfile.write('\n')   

file_result('Common_proteins')

#print(len(set(rep[0]['Protein_name']).intersection(set(rep[1]['Protein_name'])).intersection(set(rep[2]['Protein_name']))))
#index_inter = dt.intersect_index([rep[0],rep[1], rep[2], ctr[0], ctr[1], ctr[2]])
#print(len(index_inter))
#print(len(intersect([rep[0],rep[1], rep[2], ctr[0], ctr[1], ctr[2]])))
#print(len(rep[0]['Protein_name']))
