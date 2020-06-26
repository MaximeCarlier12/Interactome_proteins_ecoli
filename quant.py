import numpy as np
import pandas as pd
import os
import Load_PPI_screen as dt
import qualitatif_stats as st
import glob
import matplotlib.pyplot as plt

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

os.chdir("/home/carlier/Documents/Stage/Interactome_proteins_ecoli/")

CONDITION = ['LB log', 'LB O/N' ,'M9 0.2% ac O/N']
CONTROL = ['MG1655 (TYPE A)', 'MG1655 (placI)mVenus-SPA-pUC19 (TYPE C2)', 'MG1655']
PROTEINS = ['DnaA', 'DiaA', 'Hda', 'SeqA', 'HolD','DnaB', 'DnaG', 'NrdB']
controls_typeC = {'LB log':['S1','S2','S3'], 'LB O/N':['O9', 'R4', 'R5'], 'M9 0.2% ac O/N':['R1', 'R2', 'R3']}
controls_typeA = {'LB log':['L1', 'T7', 'T8', 'T9'], 'LB O/N':['C13', 'P5', 'U10'], 'M9 0.2% ac O/N':['A10', 'T5', 'T6']}
pd_samples = dt.load_based_screen_samples()
pd_controls = dt.load_based_screen_controls()

def load_df(bfNumber):
  '''From a batch file number, it gets the right corresponding dataframe.'''
  path_batch = "maxQ/SEGREGATED-20200619T092017Z-001/SEGREGATED/single files/"
  # full_path get the only good file and can be in two different ways (containing 'A_5' or 'A5').
  if len(bfNumber) == 1:
    letter = ''
    number = bfNumber
    full_path = path_batch+"1/"+number+".xlsx"
  else :
    letter = bfNumber[0]
    number = bfNumber[1:]
    full_path = path_batch+letter+'/'+letter+number+".xlsx"
  df = pd.read_excel(full_path, sheet_name = 0, header = 0)
  df = df[df['Gene names'].notnull()]
  df = df.set_index("Gene names")
  dico = {}
  for col in df.columns:
    if (bfNumber in col or letter+'_'+number in col) and 'MS' in col:
      dico[col] = 'MS/MS_count_sample'
    if (bfNumber in col or letter+'_'+number in col) and 'LFQ' in col:
      dico[col] = 'LFQ_intensity_sample'
    if (bfNumber in col or letter+'_'+number in col) and 'Intensity' in col:
      dico[col] = 'Intensity_sample'
    if (bfNumber in col or letter+'_'+number in col) and 'Identification' in col:
      dico[col] = 'Identification_type_sample'
  df = df.rename(columns = dico)
  return df

def select_intensity(df, LFQ):
  if LFQ ==True:
    return df[df['LFQ_intensity_sample'] != 0]
  else:
    return df[df['Intensity_sample'] != 0]

def create_table(prot1, threshold, LFQ):
  '''Create a file containing emPAI values for each gene present at list once for test and controls replicates. Absence equals to 0.'''
  prot_batches = dt.get_batches(pd_samples, prot1[0], prot1[1])
  rep = [select_intensity(load_df(i), LFQ) for i in prot_batches]
  ctrC = [select_intensity(load_df(i), LFQ) for i in controls_typeC[prot1[1]]]
  ctrA = [select_intensity(load_df(i), LFQ) for i in controls_typeA[prot1[1]]]
  print(len(dt.intersect_index(rep)))
  indexes = list(dt.intersect_index(rep))
  feature_list = ['Rep1', 'Rep2', 'Rep3', 'CtrC1', 'CtrC2', 'CtrC3', 'CtrA1', 'CtrA2', 'CtrA3']
  if len(ctrA) == 4:
    feature_list.append('CtrA4')
  df = pd.DataFrame(0, index=indexes, columns=feature_list)
  print(df.shape)
  path_batch = "maxQ/SEGREGATED-20200619T092017Z-001/Protein_table/"
  for (i,repi) in enumerate(rep):
    for my_index, row in repi.iterrows():
      if my_index in indexes:
        if LFQ == True:
          df.loc[my_index, feature_list[i]] = row.LFQ_intensity_sample
        else:
          df.loc[my_index, feature_list[i]] = row.Intensity_sample
  for (i,repi) in enumerate(ctrC):
    repi = select_intensity(repi, LFQ)
    for my_index, row in repi.iterrows():
      if my_index in indexes:
        if LFQ == True:
          df.loc[my_index, feature_list[3+i]] = row.LFQ_intensity_sample
        else:
          df.loc[my_index, feature_list[3+i]] = row.Intensity_sample
  for (i,repi) in enumerate(ctrA):
    repi = select_intensity(repi, LFQ)
    for my_index, row in repi.iterrows():
      if my_index in indexes:
        if LFQ == True:
          df.loc[my_index, feature_list[6+i]] = row.LFQ_intensity_sample
        else:
          df.loc[my_index, feature_list[6+i]] = row.Intensity_sample
  df = df.apply(lambda x: np.where(x < threshold ,threshold,x))
  if LFQ == True:
    df.to_csv(path_batch+ prot1[0]+"_"+prot1[1].replace('/', '_')+'LFQ.csv')
  else:
    df.to_csv(path_batch+ prot1[0]+"_"+prot1[1].replace('/', '_')+'Int.csv')

def load_df_table(prot1, LFQ):
  '''Load dataframe from new files with all gene names or unique gene name by adding unique = True.'''
  path_batch = "maxQ/SEGREGATED-20200619T092017Z-001/Protein_table/"
  if LFQ == True:
    full_path = path_batch+ prot1[0]+"_"+prot1[1].replace('/', '_')+'LFQ.csv'
  else:
    full_path = path_batch+ prot1[0]+"_"+prot1[1].replace('/', '_')+'Int.csv'
  df = pd.read_csv(full_path, sep=',', header=0, index_col = 0)
  return df

def log_calculation(prot1, LFQ):
  '''Calculation of emPAI logs and for a specific protein/condition.'''
  path_batch = "maxQ/SEGREGATED-20200619T092017Z-001/Protein_table/"
  df = load_df_table(prot1, LFQ)
  if 'CtrC4' in df:
    df = df[['Rep1', 'Rep2', 'Rep3', 'CtrC1', 'CtrC2', 'CtrC3', 'CtrA1','CtrA2', 'CtrA3', 'CtrA4']]
  else:
    df = df[['Rep1', 'Rep2', 'Rep3', 'CtrC1', 'CtrC2', 'CtrC3', 'CtrA1','CtrA2', 'CtrA3']]
  for i in df.columns:
    newname = i+'_log10'
    df[newname] = np.log10(df[i])
  if LFQ == True:
    df.to_csv(path_batch+ prot1[0]+"_"+prot1[1].replace('/', '_')+'LFQ.csv')
  else:
    df.to_csv(path_batch+ prot1[0]+"_"+prot1[1].replace('/', '_')+'Int.csv')


def plot_log10_emPAI(prot1, threshold, contaminant_genes, control = 'AC'):
  '''Plot log10(emPAI) value for each gene for controls and test.'''

  fig,ax = plt.subplots()
  width = 12.5 ; height = 6.35 # taille finale de ta figure png
#  width = 14.4 ; height = 7.15 # taille finale de ta figure svg
  fig.set_size_inches(width, height)
#  var_test = hd.load_df_equal_test()
  df = load_df_table(prot1, True)
  indexes = df.index
  df['max_intensity'] = df[['Rep1', 'Rep2', 'Rep3']].max(axis=1)
  minval = df[['Rep1', 'Rep2', 'Rep3']].min().min() # min value of the table
  maxval = df[['Rep1', 'Rep2', 'Rep3']].max().max() # max value of the table
  df = df.sort_values(by = 'max_intensity', ascending = False)
  for i,rep in enumerate(['Rep1_log10', 'Rep2_log10', 'Rep3_log10']):
    ax.scatter(df.index, df[rep], label="Protein test with LFQ_Intensity" if i == 0 else "", color='royalblue', alpha=0.5, marker = 'o', s=40)
  plt.title(prot1[0]+' in '+prot1[1]+' with intensity = '+str(threshold)+' as threshold ')
  if 'A' in control:
    for i,rep in enumerate(['CtrA1_log10', 'CtrA2_log10', 'CtrA3_log10']):
      ax.scatter(df.index, df[rep], label="CtrA (without SPA tag) with LFQ" if i == 0 else "", color='red', alpha=0.6, marker = 0, s=40)
    if 'CtrA4' in df :
      ax.scatter(df.index, np.log10(df.CtrA4), color='red', alpha=0.6, marker = 0, s=40)
  if 'C' in control:
    for i,rep in enumerate(['CtrC1_log10', 'CtrC2_log10', 'CtrC3_log10']):
      ax.scatter(df.index, df[rep], label="CtrC (with SPA tag) with LFQ" if i == 0 else "", color='yellowgreen', alpha=0.6, marker = 1, s=40)

  df_all = load_df_table(prot1, False)
  df = df_all[df_all.index.isin(indexes)]
  df_remove = df_all[~df_all.index.isin(indexes)]
  df_remove['max_intensity'] = df_remove[['Rep1', 'Rep2', 'Rep3']].max(axis=1)
  df_remove = df_remove.sort_values(by = 'max_intensity', ascending = False)
  for i, rep in enumerate(['Rep1_log10', 'Rep2_log10', 'Rep3_log10']):
    ax.scatter(df.index, df[rep], label="Protein test with Intensity" if i == 0 else "", color='magenta', alpha=0.3, marker = 'o', s=40)
  list_others = ['others1', 'others2', 'others3']
  for i, rep in enumerate(['Rep1_log10', 'Rep2_log10', 'Rep3_log10']):
    ax.scatter(list_others*(df_remove.shape[0]//3)+list_others[0:df_remove.shape[0]%3], df_remove[rep], color='green', label="Protein only identified with Intensity for "+str(df_remove.shape[0])+' other genes' if i == 0 else "", alpha=0.3, marker = 'o', s=40)
  fig.canvas.draw()
  plt.ylim(np.log10(minval)-0.3, np.log10(maxval)+0.3)
  plt.xlabel('Gene name')
  plt.ylabel('log10(intensity) value')
  plt.grid(axis = 'x') # vertical lines
  plt.xticks(rotation=90)
  plt.legend()
  path_batch = "maxQ/Images/"
# get an appropriate plot and saved image.
  manager = plt.get_current_fig_manager() # get full screen
  manager.window.showMaximized() # get full screen
  fig.tight_layout()
  fig.subplots_adjust(left=.05, bottom=.2, right=.96, top=.93) # marges
  filename = path_batch+prot1[0]+'_'+prot1[1][:6].replace('/', '_')+'_'+str(threshold)+'_log10values.png'
#  plt.savefig(path_batch+'test.svg') # image vectorisée
  plt.savefig(filename, transparent=False, dpi = 300) # image pixelisée, dpi = résolution
#  plt.show()

df = load_df('F2')
df_int = select_intensity(df, False)
df_LFQ = select_intensity(df, True)

print('intensity', df_int.shape)
print('LFQ', df_LFQ.shape)
print('identif type', df[df['Identification_type_sample'].notnull()].shape)
#print(df.iloc[11])
#bn = dt.get_batches(pd_samples, 'DnaA', 'LB log') #batch names
#replicates = [load_df(bn[0]), load_df(bn[1]), load_df(bn[2])] # pandas replicates
used_prot_tuples = dt.good_proteins()
contaminant_genes = dt.load_contaminant_list()
threshold = 1000000
for prot in [used_prot_tuples[2], used_prot_tuples[3], used_prot_tuples[5]]:
  create_table(prot, threshold, False)
  create_table(prot, threshold, True)
  log_calculation(prot, False)
  log_calculation(prot, True)
  plot_log10_emPAI(prot, threshold, contaminant_genes)
