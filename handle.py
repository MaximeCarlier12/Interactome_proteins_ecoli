import Load_PPI_screen as dt
import Data_processing as dp
import csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn2
import pandas as pd
import os
#from matplotlib.ticker import FormatStrFormatter
from scipy import stats
import statsmodels.stats.multitest

params = {
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
   'axes.labelsize': 11,
   'font.size': 11,
   'legend.fontsize': 9,
   'xtick.labelsize': 10,
   'ytick.labelsize': 10,
   'text.usetex': False,
   'axes.linewidth':1.5, #0.8
   'axes.titlesize':11,
   'axes.spines.top':True,
   'axes.spines.right':True,
   'font.family': "Arial"
   }
plt.rcParams.update(params)

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
#contaminant_genes = dt.load_contaminant_list()
controls_typeC = {'LB log':['S1','S2','S3'], 'LB O/N':['O9', 'R4', 'R5'], 'M9 0.2% ac O/N':['R1', 'R2', 'R3']}
controls_typeA = {'LB log':['L1', 'T7', 'T8', 'T9'], 'LB O/N':['C13', 'P5', 'U10'], 'M9 0.2% ac O/N':['A10', 'T5', 'T6']}

def used_bnames(used_prot_tuples):
  '''List of batch names used'''
  bnames = []
  for my_tuple in used_prot_tuples:
    bprot = []
    batch_names = dt.get_batches(pd_samples,my_tuple[0], my_tuple[1])
    for brep in batch_names:
      bprot.append(brep)
    bnames.append(bprot)
  return bnames

def intersect(my_list):
  '''Returns a set, the intersection is now made on 'Gene_name' column.'''
  intersect = set(my_list[0]['Gene_name'])
  for i in my_list:
    intersect = intersect.intersection(set(i['Gene_name']))
  return intersect

def df_intersect(my_list):
  '''Returns a df pandas where each protein is present in all replicates.'''
  df = my_list[0]
  return df[df.Gene_name.isin(list(intersect(my_list)))]

def file_result(filename):
  '''Outputs all proteins in control intersect and test intersect. '''
  myfile = open(filename, 'w')
  for prot1 in used_prot_tuples:
    prot_batches = dt.get_batches(pd_samples, prot1[0], prot1[1])
    rep = [load_df_unique_gene(i) for i in prot_batches]
    ctr = [load_df_unique_gene(i) for i in controls_typeC[prot1[1]]]
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
  '''Plot a venn2 or venn3 diagram.'''
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
  '''Venn diagram of replicate of protein test for all used proteins.'''
  plt.suptitle('Venn diagrams in replicates of a protein test')
  for (i,prot1) in enumerate(used_prot_tuples):
    prot_batches = dt.get_batches(pd_samples, prot1[0], prot1[1])
    rep = [load_df_unique_gene(i) for i in prot_batches]
    plt.subplot(3,3, i+1)
    plt.title(prot1[0]+' in '+prot1[1])
    venn_diagram(rep)
  plt.show()

def venn_ctr(used_prot_tuples, controls_dico, control_type):
  '''Venn diagram of replicates of a given control for all used proteins.'''
  plt.suptitle('Venn diagrams in controls type'+control_type)
  for i, cond in enumerate(controls_dico.keys()):
    ctr = [load_df_unique_gene(bname) for bname in controls_dico[cond]]
    plt.subplot(1,3, i+1)
    plt.title('Condition '+cond)
    venn_diagram(ctr)
  plt.show()

def venn_inter(used_prot_tuples, controls_dico):
  '''Venn diagram between intersect control and test.'''
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
#venn_ctr(used_prot_tuples, controls_typeA, 'A')
#venn_ctr(used_prot_tuples, controls_typeC, 'C')
#venn_inter(used_prot_tuples, controls_typeC)

def create_table(used_prot_tuples):
  '''Create a file containing emPAI values for each gene present at list once for test and controls replicates. Absence equals to 0.'''
  for prot1 in used_prot_tuples:
    prot_batches = dt.get_batches(pd_samples, prot1[0], prot1[1])
    rep = [load_df_unique_gene(i) for i in prot_batches]
    ctrC = [load_df_unique_gene(i) for i in controls_typeC[prot1[1]]]
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

def load_df_table(prot1, multiple = False):
  '''Load dataframe from new files with all gene names or unique gene name by adding unique = True.'''
  path_batch = "MS data csv reportbuilder-20200408T212019Z-001/Used_proteins_csv/"
  if multiple == False:
    full_path = path_batch+ prot1[0]+"_"+prot1[1].replace('/', '_')+'.csv'
  else :
    full_path = path_batch+ prot1[0]+"_"+prot1[1].replace('/', '_')+'_multipleRep'+'.csv'
  df = pd.read_csv(full_path, sep=',', header=0, index_col = 0)
  return df

def create_condensed_table(used_prot_tuples, threshold):
  '''For each protein used : creates a file containing emPAI values for each gene present in 2 replicates of test protein files for test and controls.
  threshold can be a specific value, 0.01 or min_value'''
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
    min_value = df[df > 0].min().min()
#    df = df.replace(0, min_value)
    df = df.apply(lambda x: np.where(x < threshold ,threshold,x))
    indexes = []
    for index, row in df.iterrows():
      counter = 0
      if row.Rep1 == threshold:
        counter += 1
      if row.Rep2 == threshold:
        counter += 1
      if row.Rep3 == threshold:
        counter += 1
      if counter < 2:
        indexes.append(index)
    df = df.loc[df.index.isin(indexes)]
    path_batch = "MS data csv reportbuilder-20200408T212019Z-001/Used_proteins_csv/"
    df.to_csv(path_batch+ prot1[0]+"_"+prot1[1].replace('/', '_')+'_multipleRep.csv')
  return threshold

def alone_rep():
  '''Create a file containing proteins that are present in only one replicate.'''
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

def nbr_prot_file():
  '''Print number of prot per file '''
  allconditions = ['LB log', 'LB O/N','M9 0.2% ac O/N']
  feature_list = ['Rep1', 'Rep2', 'Rep3', 'CtrC1', 'CtrC2', 'CtrC3', 'CtrA1', 'CtrA2', 'CtrA3']
  nb_prot = pd.DataFrame(0, index=allconditions, columns=feature_list)
  dt.header('tests')
  for prot in used_prot_tuples:
    bname_prot = dt.get_batches(pd_samples, prot[0], prot[1])
    for (i, bname) in enumerate(bname_prot):
      df = load_df_unique_gene(bname)
      print(df.shape[0], end = ', ')
      print(i)
      nb_prot.loc[prot[1], 'Rep'+str(i+1)] = df.shape[0]
  print()
  dt.header('ctrC')
  for dico in controls_typeC:
    for i,bname in enumerate(controls_typeC[dico]):
      df = load_df_unique_gene(bname)
      print(df.shape[0], end = ', ')
      print(i)
      nb_prot.loc[dico, 'CtrC'+str(i+1)] = df.shape[0]
  print()
  dt.header('ctrA')
  for dico in controls_typeA:
    for i,bname in enumerate(controls_typeA[dico]):
      df = load_df_unique_gene(bname)
      print(df.shape[0], end = ', ')
      if i <= 2:
        nb_prot.loc[dico, 'CtrA'+str(i+1)] = df.shape[0]
  print()
  print(nb_prot)

def get_repeated_proteins(used_prot_tuples):
  '''Get proteins that are present in different samples.'''
  dico_common_proteins = {}
  for j in (0,2,5,6,7):
    prot = used_prot_tuples[j]
    df = load_df_table(prot, True)
    for i in df.index.values:
      if i in dico_common_proteins:
        dico_common_proteins[i] +=1
      else: dico_common_proteins[i] = 1
  for i in range(5):
    dt.header(str(i+1))
    for key, value in dico_common_proteins.items():
      if value == i+1:
        print(key, end = ', ')
    print()

def load_df_equal_test():
  os.chdir("/home/carlier/Documents/Stage/Interactome_proteins_ecoli/") # where is the MS data
  path = "MS data csv reportbuilder-20200408T212019Z-001/Significant_proteins/contaminant_or_not.csv" 
#  df = pd.read_csv(path, header=0, sep =';', usecols = lambda column : column not in ['Potential contaminant','Not contaminant'])
  df = pd.read_csv(path, header=0, sep =';')
  df = df.set_index(['Bait protein', 'Condition'])
  return df

#create_table(used_prot_tuples)
#threshold = create_condensed_table(used_prot_tuples, 0.25)
#set_pval(used_prot_tuples)
#for prot1 in used_prot_tuples:
#  df = load_df_table(prot1, True)
#  print(prot1[0], prot1[1][:6]+'\t', 'TestC', df[df.C_is ==True].shape[0],'\t', 'TestA', df[df.A_is ==True].shape[0],'\t', 'TestC & TestA', df[(df.A_is ==True) & (df.C_is == True)].shape[0])

#for prot in used_prot_tuples:
#  if prot[0] == "HolD":
##    if prot[1] == "LB log":
#    print(dt.get_batches(pd_samples, prot[0], prot[1]))

#plot_log2_emPAI(used_prot_tuples[1], 0.25,contaminant_genes, 'AC')

#for i in used_prot_tuples[0:1]:
#  plot_log2_emPAI(i, 0.25,contaminant_genes, 'AC')
#plt.close('all')
