import Load_PPI_screen as dt
import Data_processing as dp
from Bio import Entrez
import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn2
import pandas as pd
import os
import glob
from matplotlib.ticker import FormatStrFormatter
from scipy import stats

def load_df_unique_gene(bfNumber):
  '''Load dataframe from new files with all gene names.'''
  os.chdir("/home/carlier/Documents/Stage/Interactome_proteins_ecoli/") # where is the MS data
  path_batch = "MS data csv reportbuilder-20200408T212019Z-001/MS data csv reportbuilder/"
  letter = bfNumber[0]
  number = bfNumber[1:]
  full_path = path_batch+"batch "+letter+"/unique_gene_"+letter+number+'.csv'
  df = pd.read_csv(full_path, sep=',', header=0)
  df = df.set_index("Accession_number")
  return df

pd_samples = dt.load_based_screen_samples()
pd_controls = dt.load_based_screen_controls()
used_prot_tuples = dt.good_proteins() 
controls = {'LB log':['S1','S2','S3'], 'LB O/N':['O9', 'R4', 'R5'], 'M9 0.2% ac O/N':['R1', 'R2', 'R3']}
controls_typeA = {'LB log':['L1', 'T7', 'T8', 'T9'], 'LB O/N':['C13', 'P5', 'U10'], 'M9 0.2% ac O/N':['A10', 'T5', 'T6']}

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
  '''Returns a set '''
  intersect = set(my_list[0]['Gene_name'])
  for i in my_list:
    intersect = intersect.intersection(set(i['Gene_name']))
  return intersect

def df_intersect(my_list):
  '''Returns a df pandas '''
  df = my_list[0]
  return df[df.Gene_name.isin(list(intersect(my_list)))]

def file_result(filename):
  '''Outputs all proteins in control intersect and test intersect. '''
  myfile = open(filename, 'w')
  for prot1 in used_prot_tuples:
    prot_batches = dt.get_batches(pd_samples, prot1[0], prot1[1])
    rep = [load_df_unique_gene(i) for i in prot_batches]
    ctr = [load_df_unique_gene(i) for i in controls[prot1[1]]]
    intR = intersect(rep)
    intC = intersect(ctr)
    pattern = '|'.join(intersect(ctr))
    myfile.write('(Protein) '+ prot1[0]+ ' (Condition) '+ prot1[1] + ':\n')
    myfile.write('\tLen rep intersect '+ str(len(intR))+ '\n')
    pNamesR = rep[0].Protein_name.values[rep[0].Gene_name.isin(list(intR))]
    myfile.write('\tCommon rep proteins :\n')    
    for i in pNamesR:
      myfile.write(i+'\n')
    myfile.write('\tLen ctr intersect '+ str(len(intC))+'\n')
    pNamesC = ctr[0].Protein_name.values[ctr[0].Gene_name.isin(list(intC))]
    myfile.write('\tCommon ctr proteins :\n')  
    for i in pNamesC:
      myfile.write(i+'\n')   
    myfile.write('\n')   

def venn_diagram(data):
  set_array = []
  for rep in data:
    set_array.append(set(rep['Gene_name']))
  if len(data) == 3:
    set_names = ['Rep1', 'Rep2', 'Rep3']
    venn3(set_array, set_names)   # venn3 works for three sets
  elif len(data) == 2:
    set_names = ['Rep', 'Ctr']
    venn2(set_array, set_names)   # venn3 works for three sets
  elif len(data) == 4:
    set_names = ['Rep1', 'Rep2', 'Rep3']
    venn3(set_array[:3], set_names)   # venn3 works for three sets
  else : print('error, please change data length')

def venn_rep(used_prot_tuples):
  plt.suptitle('Venn diagrams in replicates')
  for (i,prot1) in enumerate(used_prot_tuples):
    prot_batches = dt.get_batches(pd_samples, prot1[0], prot1[1])
    rep = [load_df_unique_gene(i) for i in prot_batches]
    plt.subplot(3,3, i+1)
    plt.title(prot1[0]+' in '+prot1[1])
    venn_diagram(rep)
  plt.show()

def venn_ctr(used_prot_tuples, controls_dico):
  plt.suptitle('Venn diagrams in controls')
  for (i,prot1) in enumerate(used_prot_tuples):
    prot_batches = dt.get_batches(pd_samples, prot1[0], prot1[1])
    ctr = [load_df_unique_gene(i) for i in controls_dico[prot1[1]]]
    plt.subplot(3,3, i+1)
    plt.title(prot1[0]+' in '+prot1[1])
    venn_diagram(ctr)
  plt.show()

def venn_inter(used_prot_tuples, controls_dico):
  plt.suptitle('Venn diagrams of intersection of controls and replicates')
  for (i,prot1) in enumerate(used_prot_tuples):
    prot_batches = dt.get_batches(pd_samples, prot1[0], prot1[1])
    rep = [load_df_unique_gene(i) for i in prot_batches]
    interR = df_intersect(rep)
    ctr = [load_df_unique_gene(i) for i in controls_dico[prot1[1]]]
    interC = df_intersect(ctr)
    plt.subplot(3,3, i+1)
    plt.title(prot1[0]+' in '+prot1[1])
    venn_diagram([interR, interC])
  plt.show()

#venn_rep(used_prot_tuples)
#venn_ctr(used_prot_tuples, controls_typeA)
#venn_inter(used_prot_tuples, controls)

def create_table(used_prot_tuples):
  for prot1 in used_prot_tuples:
    prot_batches = dt.get_batches(pd_samples, prot1[0], prot1[1])
    rep = [load_df_unique_gene(i) for i in prot_batches]
    ctrC = [load_df_unique_gene(i) for i in controls[prot1[1]]]
    ctrA = [load_df_unique_gene(i) for i in controls_typeA[prot1[1]]]
    all_genes = []
    for i in rep:
      all_genes += list(i.Gene_name)
    for i in ctrC:
      all_genes += list(i.Gene_name)
    for i in ctrA:
      all_genes += list(i.Gene_name)
    all_genes = list(set(all_genes))
    all_genes.sort(key=str.lower)
    feature_list = ['Rep1', 'Rep2', 'Rep3', 'CtrC1', 'CtrC2', 'CtrC3', 'CtrA1', 'CtrA2', 'CtrA3']
    if len(ctrA) == 4:
      feature_list.append('CtrA4')
    df = pd.DataFrame(0, index=all_genes, columns=feature_list)
    path_batch = "MS data csv reportbuilder-20200408T212019Z-001/Used_proteins_csv/"
    for (i,repi) in enumerate(rep):
      for index, row in repi.iterrows():
        df.loc[row.Gene_name, feature_list[i]] = row.emPAI
    for (i,repi) in enumerate(ctrC):
      for index, row in repi.iterrows():
        df.loc[row.Gene_name, feature_list[3+i]] = row.emPAI
    for (i,repi) in enumerate(ctrA):
      for index, row in repi.iterrows():
        df.loc[row.Gene_name, feature_list[6+i]] = row.emPAI
    df.to_csv(path_batch+ prot1[0]+"_"+prot1[1].replace('/', '_')+'.csv')

def load_df_table(prot1, string = ''):
  '''Load dataframe from new files with all gene names.'''
  path_batch = "MS data csv reportbuilder-20200408T212019Z-001/Used_proteins_csv/"
  full_path = path_batch+ prot1[0]+"_"+prot1[1].replace('/', '_')+string+'.csv'
  df = pd.read_csv(full_path, sep=',', header=0, index_col = 0)
  return df

def create_table_condensed_table(used_prot_tuples):
  for prot1 in used_prot_tuples:
    df = load_df_table(prot1)
    indexes = []
    for index, row in df.iterrows():
      counter = 0
      if row.Rep1 == 0:
        counter += 1
      if row.Rep2 == 0:
        counter += 1
      if row.Rep3 == 0:
        counter += 1
      if counter > 1:
        indexes.append(index)
    df = df.drop(indexes)
    path_batch = "MS data csv reportbuilder-20200408T212019Z-001/Used_proteins_csv/"
    df.to_csv(path_batch+ prot1[0]+"_"+prot1[1].replace('/', '_')+'_multipleRep.csv')

#create_table_condensed_table(used_prot_tuples)

def plot_emPAI(prot1, control = 'AC'):
  df = load_df_table(prot1, '_multipleRep')
  empRep = []; empCtr = []
  df['max_empai'] = df[['Rep1', 'Rep2', 'Rep3']].max(axis=1)
  df = df.sort_values(by = 'max_empai', ascending = False)
  a = max(df.max(axis = 1)) # max value of the table
  #  empRep.append(np.mean([row.Rep1, row.Rep2, row.Rep3]))
  #  empCtr.append(np.mean([row.Ctr1, row.Ctr2, row.Ctr3]))
  plt.scatter(df.index, df.Rep1, label="Protein test", color='royalblue')
  plt.scatter(df.index, df.Rep2, color='royalblue')
  plt.scatter(df.index, df.Rep3, color='royalblue')
  plt.title(prot1[0]+' in '+prot1[1])
  if 'A' in control:
    plt.scatter(df.index, df.CtrA1, label = "CtrA (without SPA tag)", color='red')
    plt.scatter(df.index, df.CtrA2, color='red')
    plt.scatter(df.index, df.CtrA3, color='red')
    if 'CtrA4' in df :
      plt.scatter(df.index, df.CtrA4, color='red')
  if 'C' in control:
    plt.scatter(df.index, df.CtrC1, label = "CtrC (with SPA tag)", color='chartreuse')
    plt.scatter(df.index, df.CtrC2, color='chartreuse')
    plt.scatter(df.index, df.CtrC3, color='chartreuse')
  plt.xticks(rotation=90)
  plt.ylim(-0.02, 100)
  plt.xlabel('Gene name')
  plt.ylabel('emPAI value')
  plt.yscale('symlog')
  plt.grid(axis = 'x') # vertical lines
  plt.legend()
  manager = plt.get_current_fig_manager()
  manager.window.showMaximized()
  plt.tight_layout()
  plt.show()

def plot_log2_emPAI(prot1, control = 'AC'):
  df = load_df_table(prot1, '_multipleRep')
  empRep = []; empCtr = []
  df['max_empai'] = df[['Rep1', 'Rep2', 'Rep3']].max(axis=1)
  maxval = max(df.max(axis = 1)) # max value of the table
  min
  df = df.sort_values(by = 'max_empai', ascending = False)
  df = df.replace(0, 0.01)
  #  empRep.append(np.mean([row.Rep1, row.Rep2, row.Rep3]))
  #  empCtr.append(np.mean([row.Ctr1, row.Ctr2, row.Ctr3]))
  plt.scatter(df.index, np.log2(df.Rep1), label="Protein test", color='royalblue')
  plt.scatter(df.index, np.log2(df.Rep2), color='royalblue')
  plt.scatter(df.index, np.log2(df.Rep3), color='royalblue')
  plt.title(prot1[0]+' in '+prot1[1])
  if 'A' in control:
    plt.scatter(df.index, np.log2(df.CtrA1), label = "CtrA (without SPA tag)", color='red')
    plt.scatter(df.index, np.log2(df.CtrA2), color='red')
    plt.scatter(df.index, np.log2(df.CtrA3), color='red')
    if 'CtrA4' in df :
      plt.scatter(df.index, np.log2(df.CtrA4), color='red')
  if 'C' in control:
    plt.scatter(df.index, np.log2(df.CtrC1), label = "CtrC (with SPA tag)", color='chartreuse')
    plt.scatter(df.index, np.log2(df.CtrC2), color='chartreuse')
    plt.scatter(df.index, np.log2(df.CtrC3), color='chartreuse')
  plt.xticks(rotation=90)
  plt.ylim(np.log2(0.01)-0.5, np.log2(maxval)+1)
  plt.xlabel('Gene name')
  plt.ylabel('log2(emPAI) value')
#  plt.yscale('symlog')
#  plt.yticks([0, 1, -1, np.log2(maxval), np.log2(0.01)])
#  plt.gca().yaxis.set_major_formatter(FormatStrFormatter('%.d'))
  plt.grid(axis = 'x') # vertical lines
  plt.legend()
  manager = plt.get_current_fig_manager()
  manager.window.showMaximized()
  plt.tight_layout()
  plt.show()
  os.chdir("/home/carlier/Documents/Stage/Images/emPAI/") # where is the MS data
  os.chdir("/home/carlier/Documents/Stage/Interactome_proteins_ecoli/") # where is the MS data

#plot_emPAI(used_prot_tuples[0], 'AC')
#plot_log2_emPAI(used_prot_tuples[0], 'AC')

def alone_rep():
  prot1 = used_prot_tuples[0]
  df = load_df_table(prot1)
  indexes = []
  for index, row in df.iterrows():
    counter = 0
    if row.Rep1 == 0:
      counter += 1
    if row.Rep2 == 0:
      counter += 1
    if row.Rep3 == 0:
      counter += 1
    if counter > 1:
      indexes.append(index)
  prot_batches = dt.get_batches(pd_samples, prot1[0], prot1[1])
  rep = [load_df_unique_gene(i) for i in prot_batches]
  df_select = rep[0].loc[rep[0].Gene_name.isin(indexes)]
  path_batch = "MS data csv reportbuilder-20200408T212019Z-001/Used_proteins_csv/"
  df_select.to_csv(path_batch+ prot1[0]+"_"+prot1[1].replace('/', '_')+'_alone.csv')

def set_pval(df):
  prot1 = used_prot_tuples[0]
  df = load_df_table(prot1, '_multipleRep')
  testC = stats.ttest_ind(df[['Rep1', 'Rep2', 'Rep3']], df[['CtrC1', 'CtrC2', 'CtrC3']], axis = 1)
  if 'CtrA4' in df:
    testA = stats.ttest_ind(df[['Rep1', 'Rep2', 'Rep3']], df[['CtrA1', 'CtrA2', 'CtrA3', 'CtrA4']], axis = 1)
  else :
    testA = stats.ttest_ind(df[['Rep1', 'Rep2', 'Rep3']], df[['CtrA1', 'CtrA2', 'CtrA3']], axis = 1)
  df['tstat_C'] = testC[0]
  df['pval_C'] = testC[1]
  df['tstat_A'] = testA[0]
  df['pval_A'] = testA[1]
  path_batch = "MS data csv reportbuilder-20200408T212019Z-001/Used_proteins_csv/"
  df.to_csv(path_batch+ prot1[0]+"_"+prot1[1].replace('/', '_')+'_multipleRep.csv')

prots = []
name_prots = []
for prot1 in used_prot_tuples:
  df = load_df_table(prot1, '_multipleRep')
  sums = df[['Rep1', 'Rep2', 'Rep3', 'CtrC1', 'CtrC2', 'CtrC3', 'CtrA1', 'CtrA2', 'CtrA3']].sum(axis=0)
  name = prot1[0]+'_'+prot1[1][:4]
  name_prots.append(name)
  prots.append(sums)
  print(sums[0])
for i in range(len(prots)):
  plt.scatter([name_prots[i]]*3, prots[i][:3] , label="Protein test" if i == 0 else "", color='royalblue')
  plt.scatter([name_prots[i]]*3, prots[i][3:6] , label="CtrC" if i == 0 else "", color='chartreuse')
  plt.scatter([name_prots[i]]*3, prots[i][6:] , label="CtrA" if i == 0 else "", color='red')
plt.legend()
plt.title('Differences of emPAI sums between files')
plt.xlabel('Protein studied')
plt.ylabel('sum(emPAI) value')
plt.xticks(rotation=90)
plt.grid(axis = 'x') # vertical lines
manager = plt.get_current_fig_manager()
manager.window.showMaximized()
plt.tight_layout()
plt.show()
