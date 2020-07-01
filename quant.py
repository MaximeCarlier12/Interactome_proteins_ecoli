import numpy as np
import pandas as pd
import os
import Load_PPI_screen as dt
import Data_processing as dp
import qualitatif_stats as st
import handle as hd
import glob
import matplotlib.pyplot as plt
from scipy import stats
import statsmodels.stats.multitest
import matplotlib.patches as mpatches

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
missing_files = ['A1', 'A2', 'M4', 'M5', 'A9', 'A10', 'I1', 'I2']
prey = {'DnaA':'diaA', 'DiaA':'dnaA', 'Hda' : '', 'DnaB':'dnaC', 'DnaG':'dnaB', 'NrdB':'nrdA', 'HolD':['dnaE', 'dnaN', 'dnaQ', 'dnaX', 'holA', 'holB', 'holC', 'holE'], 'SeqA':''}
interesting_prey = {'DnaA':['purR', 'eno'], 'DiaA':['gldA', 'ndh', 'wbbK', 'rfaF', 'rfaB', 'rfaG','rfaP','RfaD','rfaB','gmhA'], 'Hda':'', 'DnaB':'tdcB', 'DnaG':['nrdB', 'glgB', 'amyA', 'glgA', 'seqA'], 'NrdB':['dnaN', 'skp'], 'HolD':['topB'], 'SeqA':'rfaD'}
contaminant_genes = dt.load_contaminant_list()
pd_samples = dt.load_based_screen_samples()
pd_controls = dt.load_based_screen_controls()

def all_proteins():
  '''List of tuple containing the protein to be used and its condition. Each tuple contains the protein name and the growth condition.'''
  good_proteins = []
  for p in PROTEINS :
    for c in CONDITION :
      good_proteins.append((p,c))
  return good_proteins

def get_batches_missing_files(prot_batches):
  '''Some files should not be treated so we don't use them so far. '''
  missing_files = ['A1', 'A2', 'M4', 'M5', 'A9', 'A10', 'I1', 'I2']
  for bn in prot_batches: 
    if bn in missing_files:
      prot_batches.remove(bn)
  return prot_batches

def prot_two_rep():
  '''List of tuple containing the protein to be used and its condition. Each tuple contains the protein name and the growth condition. Here are only present proteins with two replicates.'''
  good_proteins = []
  for p in PROTEINS :
    for c in CONDITION :
      prot_batches = dt.get_batches(pd_samples, p, c)
      prot_batches = get_batches_missing_files(prot_batches)
      if len(prot_batches) == 2:
        good_proteins.append((p,c))
  return good_proteins

def prot_three_rep():
  '''List of tuple containing the protein to be used and its condition. Each tuple contains the protein name and the growth condition. Some files are missing and we only have 2 replicats so we remove this protein/condition from our analysis so far.'''
  good_proteins = []
  for p in PROTEINS :
    for c in CONDITION :
      prot_batches = dt.get_batches(pd_samples, p, c)
      if len(get_batches_missing_files(prot_batches)) >= 3 :
        good_proteins.append((p,c))
  return good_proteins

def load_df_maxQ(bfNumber):
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
  # df = df[df['Majority protein IDs'].notnull()]
  # df = df.set_index("Majority protein IDs") # we keep when gene_name == None so far.
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
  '''Either we can choose LFQ intensity or raw intensity.'''
  if LFQ ==True:
    return df[df['LFQ_intensity_sample'] != 0]
  else:
    return df[df['Intensity_sample'] != 0]

def create_table(prot1, threshold, LFQ):
  '''Create a file containing intensity values for each gene present in our 3 test replicates. Absence equals to threshold.'''
  prot_batches = dt.get_batches(pd_samples, prot1[0], prot1[1])
  prot_batches = get_batches_missing_files(prot_batches)
  rep = [select_intensity(load_df_maxQ(i), LFQ) for i in prot_batches]
  ctrC = [select_intensity(load_df_maxQ(i), LFQ) for i in controls_typeC[prot1[1]]]
  ctrA = [select_intensity(load_df_maxQ(i), LFQ) for i in controls_typeA[prot1[1]]]
  print(len(dp.intersect_index(rep)))
  indexes = list(dp.intersect_index(rep))
  feature_list = ['Rep1', 'Rep2', 'Rep3', 'CtrC1', 'CtrC2', 'CtrC3', 'CtrA1', 'CtrA2', 'CtrA3']
  if len(ctrA) == 4:
    feature_list.append('CtrA4')
  df = pd.DataFrame(0, index=indexes, columns=feature_list)
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
  for col in df[['Rep1', 'Rep2', 'Rep3']].columns:
    df = df.loc[df[col] != threshold] # remove each row where one replicate value equals to threshold
  if LFQ == True:
    df.to_csv(path_batch+ prot1[0]+"_"+prot1[1].replace('/', '_')+'_LFQ.csv')
  else:
    df.to_csv(path_batch+ prot1[0]+"_"+prot1[1].replace('/', '_')+'_Int.csv')

def log10_calculation(prot1, LFQ):
  '''Calculation of intensity logs10 and for a specific protein/condition.'''
  path_batch = "maxQ/SEGREGATED-20200619T092017Z-001/Protein_table/"
  df = hd.load_df_table_maxQ(prot1, LFQ)
  if 'CtrC4' in df:
    df = df[['Rep1', 'Rep2', 'Rep3', 'CtrC1', 'CtrC2', 'CtrC3', 'CtrA1','CtrA2', 'CtrA3', 'CtrA4']]
  else:
    df = df[['Rep1', 'Rep2', 'Rep3', 'CtrC1', 'CtrC2', 'CtrC3', 'CtrA1','CtrA2', 'CtrA3']]
  for i in df.columns:
    newname = i+'_log10'
    df[newname] = np.log10(df[i])
  if LFQ == True:
    df.to_csv(path_batch+ prot1[0]+"_"+prot1[1].replace('/', '_')+'_LFQ.csv')
  else:
    df.to_csv(path_batch+ prot1[0]+"_"+prot1[1].replace('/', '_')+'_Int.csv')

def normalize_median(prot, grand_med):
  df = hd.load_df_table_maxQ(prot, False)
  for col in df.columns:
    if 'log10' in col:
      med = df.loc[col].median()
      diff_med = med - grand_med
      df.loc[col] = df.loc[col] - diff_med

def get_samples_bartlett_test(prot, data, threshold, LFQ):
  '''Statistical test that checks if variance is similar or not between predatory proteins of a similar bait protein. H0 is equal variance.
  data can be 0 (test) or 1(ctrA) or 2(ctrC).
  Returns a list that contains a list of log10(intensity) values for each protein.'''
  if data == 0: #test
    data_type = 'test'
    df = hd.load_df_table_maxQ(prot, LFQ)
    df = df[['Rep1', 'Rep2', 'Rep3']]
  elif data == 1: #ctrA
    data_type = 'ctrlA'
    df = hd.load_df_table_maxQ(prot, LFQ)
    if 'CtrA4' in df.columns:
      df = df[['CtrA1', 'CtrA2', 'CtrA3', 'CtrA4']]
    else:
      df = df[['CtrA1', 'CtrA2', 'CtrA3']]
  elif data == 2: #ctrC
    data_type = 'ctrlC'
    df = hd.load_df_table_maxQ(prot, LFQ)
    df = df[['CtrC1', 'CtrC2', 'CtrC3']]
  else : print('error in data variable')
  df = df.apply(lambda x: np.where(x < threshold ,threshold,x))
  for i in df.columns:
    df = df.loc[df[i] != threshold]
  print('Number of proteins for '+df.columns[0]+' : '+str(df.shape[0]))
  df = np.log10(df)
  b_test = []
#  print(df.shape)
  nul_var_idx = []
  for i, my_index in enumerate(df.index): # remove row when var = 0.
    nul_var = [df.iloc[i,j] == df.iloc[i,j+1] for j in range(len(df.columns)-1)]
#    print(nul_var)
    if nul_var.count(True) == len(nul_var):
      nul_var_idx.append(my_index)
#      print('ok')
  df = df.drop(nul_var_idx)
#  print(df.shape)
#  print('Number of different proteins :', df.shape[0])
  df_norm = df.sub(df.mean(axis=1), axis=0)
#  print('variance :', np.mean(df_norm.var(axis=1)))
  my_array = [list(df.loc[i]) for i in df.index] #lists of 3 values 
  # print(my_array)
  return my_array

def test_bart(my_array):
  '''Realization of the test.'''
  return stats.bartlett(*my_array)

def test_bart_per_all_prots(used_prot_tuples, threshold):
  '''All prots from all condition with test and controls.'''
  for prot in used_prot_tuples:
    dt.header(prot[0]+' in '+prot[1])
    array = []
    # array += get_samples_bartlett_test(prot,1,threshold, False)
    # array += get_samples_bartlett_test(prot,2,threshold, False)
    array += get_samples_bartlett_test(prot,0,threshold, False)
    array += get_samples_bartlett_test(prot,1,threshold, False)
    array += get_samples_bartlett_test(prot,2,threshold, False)
    # if prot[0] == 'DnaG':
      # print(array)
    print('Bartlett test for test and ctrA and ctrC :', test_bart(array)[1])
    # print('Bartlett test for ctrA :', test_bart(get_samples_bartlett_test(prot,1,threshold, False))[1])
    # print('Bartlett test for ctrC :', test_bart(get_samples_bartlett_test(prot,2,threshold, False))[1])
    # print('Raw :', test_bart(array)[1])

def plot_var_repartition(used_prot_tuples, threshold):
  '''Shows variance repartition for each protein to see if a normal distribution suits well.'''
  step = 0
  fig = plt.figure()
  fig.suptitle('Variance repartition per bait protein for tests')
  for i,prot in enumerate(used_prot_tuples):
    if i == 9:
      step = 9
      plt.subplots_adjust(hspace=0.5, wspace=1.0)
      plt.show()
      fig = plt.figure()
      fig.suptitle('Variance repartition per bait protein for tests')
    array = []
    array += get_samples_bartlett_test(prot,0,threshold, False)
    array += get_samples_bartlett_test(prot,1,threshold, False)
    array += get_samples_bartlett_test(prot,2,threshold, False)
    list_var = []
    for j in array:
      list_var.append(np.var(j))
    ax = fig.add_subplot(3,3,1+i-step)
    ax.set_title(prot[0]+' in '+prot[1])
    # ax.scatter(labx, np.sort(list_var))
    ax.hist(list_var, bins = 20)
  plt.subplots_adjust(hspace=0.5, wspace=1.0)
  plt.show()

def non_eq_var_ttest(prot1, threshold, LFQ):
  '''Welch test = (mA-mB)/sqrt(sA^2/nA+sB^2/nB) against both controls, one at a time. We only take the variance between our 3 replicates.'''
  df = hd.load_df_table_maxQ(prot1, LFQ)
  is_fourth = 0
  if 'CtrA4' in df:
    is_fourth = 1
  df_rep = df[['Rep1_log10', 'Rep2_log10', 'Rep3_log10']]
  if is_fourth == 0:
    df_ctrA = df[['CtrA1_log10', 'CtrA2_log10', 'CtrA3_log10']]
  else:
    df_ctrA = df[['CtrA1_log10', 'CtrA2_log10', 'CtrA3_log10', 'CtrA4_log10']]
  df_ctrC = df[['CtrC1_log10', 'CtrC2_log10', 'CtrC3_log10']]
  mRep = df_rep.mean(axis = 1)
  mCtrA = df_ctrA.mean(axis = 1)
  mCtrC = df_ctrC.mean(axis = 1)
  log_threshold = np.log10(threshold)
  # df_varrep = df_rep.loc[(df_rep['Rep1_log10'] != log_threshold) & (df_rep['Rep2_log10'] != log_threshold) & (df_rep['Rep3_log10'] != log_threshold)] # keep only when 3 replicates contain the protein. 
  # df_varctrA = df_ctrA
  # for i in df_ctrA.columns:
    # df_varctrA = df_varctrA.loc[df_varctrA[i] != log_threshold]
  # df_varctrC = df_ctrC
  # for i in df_ctrC.columns:
    # df_varctrC = df_varctrC.loc[df_varctrC[i] != log_threshold]
  # df['var_Rep'] = np.mean(df_varrep.var(axis = 1)) # same value for the whole file
  # df['var_ctrA'] = np.mean(df_varctrA.var(axis = 1))
  # df['var_ctrC'] = np.mean(df_varctrC.var(axis = 1))
  # print('var_rep', df['var_Rep'][0])
  # print('varC', df['var_ctrC'][0])
  # print('varA', df['var_ctrA'][0])
  df['tvalueA'] = (mRep - mCtrA)/np.sqrt((df_rep.var(axis = 1)/len(df_rep.columns)+df_ctrA.var(axis = 1)/len(df_ctrA.columns)))
  df['tvalueC'] = (mRep - mCtrC)/np.sqrt((df_rep.var(axis=1)/len(df_rep.columns)+df_ctrC.var(axis = 1)/len(df_ctrC.columns))) # df.var is with unbiaised N-1 by default.
#  print(df['tvalueA'])
  varrep = df_rep.var(axis=1)
  varA = df_ctrA.var(axis = 1)
  nA = len(df_ctrA.columns)
  varC = df_ctrC.var(axis = 1)
  # dfA = np.power(varA/nA+varrep/3, 2)/(np.power(varrep,4)/(3*3*2)+np.power(varA,4)/(nA*nA*(nA-1)))
  dfA = np.power(varA/nA+varrep/3, 2)/(varrep/(3*2)+varA/(nA*(nA-1)))
  print(dfA)
  # df['pvalA'] = stats.t.sf(df['tvalueA'], df= len(df_rep.columns)+len(df_ctrA.columns)-2)
  df['pvalA'] = stats.t.sf(df['tvalueA'], df= dfA)
#  one sided test with the survival function of student = 1-Fct de rep°.
  df['pvalC'] = stats.t.sf(df['tvalueC'], df= len(df_rep.columns)+len(df_ctrC.columns)-2)
  print(len(df_rep.columns)+len(df_ctrC.columns)-2)
  print(stats.t.sf(1.20727, 4))
  #  print(df['pvalA'])

# correction is only on pval used (not where absent of controls).
  testAdone = df_ctrA
  for i in df_ctrA.columns:
    testAdone = testAdone.loc[testAdone[i] == log_threshold]
  testCdone = df_ctrC
  for i in df_ctrC.columns:
    testCdone = testCdone.loc[testCdone[i] == log_threshold]
#  print(testCdone.shape[0] == len(df.loc[df.index.isin(list(testCdone.index.values))]['pvalC']))
  list_pval_C = df.loc[df.index.isin(list(testCdone.index.values))]['pvalC']
  list_pval_A = df.loc[df.index.isin(list(testAdone.index.values))]['pvalA']
  pval_C_corr = statsmodels.stats.multitest.multipletests(df['pvalC'], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
  pval_A_corr = statsmodels.stats.multitest.multipletests(df['pvalA'], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
#  print(testCdone)
  df['pvalC_corr'] = pval_C_corr[1]
  df['C_is'] = pval_C_corr[0]
  df['pvalA_corr'] = pval_A_corr[1]
  df['A_is'] = pval_A_corr[0]
  def pval_is (row, ctrl):
    if row[ctrl] <=0.05 :
      return True
    else: return False 
#  df['C_is'] = df.apply (lambda row: pval_is(row, 'pvalC'), axis=1)
#  df['A_is'] = df.apply (lambda row: pval_is(row, 'pvalA'), axis=1)
  path_batch = "maxQ/SEGREGATED-20200619T092017Z-001/Protein_table/"
  if LFQ == False:
    df.to_csv(path_batch+ prot1[0]+"_"+prot1[1].replace('/', '_')+'_Int.csv')
  else:
    df.to_csv(path_batch+ prot1[0]+"_"+prot1[1].replace('/', '_')+'_LFQ.csv')    


def plot_log10_abundance_comparison(prot1, threshold, LFQ, contaminant_genes, control = 'AC'):
  '''Plot log10(emPAI) value for each gene for controls and test.'''
  fig,ax = plt.subplots()
  width = 12.5 ; height = 6.35 # taille finale de ta figure png
#  width = 14.4 ; height = 7.15 # taille finale de ta figure svg
  fig.set_size_inches(width, height)
#  var_test = hd.load_df_equal_test()
  df = hd.load_df_table_maxQ(prot1, LFQ)
  indexes = df.index
  df['max_intensity'] = df[['Rep1', 'Rep2', 'Rep3']].max(axis=1)
  minval = df[['Rep1', 'Rep2', 'Rep3']].min().min() # min value of the table
  maxval = df[['Rep1', 'Rep2', 'Rep3']].max().max() # max value of the table
  df = df.sort_values(by = 'max_intensity', ascending = False)
  if LFQ == True:
    title_text1 = 'LFQ'
    title_text2 = 'raw'
  else:
    title_text1 = 'raw'
    title_text2 = 'LFQ'
  for i,rep in enumerate(['Rep1_log10', 'Rep2_log10', 'Rep3_log10']):
    ax.scatter(df.index, df[rep], label="Protein test with "+title_text1+' intensity' if i == 0 else "", color='royalblue', alpha=0.5, marker = 'o', s=40)
  plt.title(prot1[0]+' in '+prot1[1]+' with intensity = '+str(threshold)+' as threshold ')
  if 'A' in control:
    for i,rep in enumerate(['CtrA1_log10', 'CtrA2_log10', 'CtrA3_log10']):
      ax.scatter(df.index, df[rep], label="CtrA (without SPA tag) with "+title_text1+ ' intensity' if i == 0 else "", color='red', alpha=0.6, marker = 0, s=40)
    if 'CtrA4' in df :
      ax.scatter(df.index, np.log10(df.CtrA4), color='red', alpha=0.6, marker = 0, s=40)
  if 'C' in control:
    for i,rep in enumerate(['CtrC1_log10', 'CtrC2_log10', 'CtrC3_log10']):
      ax.scatter(df.index, df[rep], label="CtrC (with SPA tag) with "+title_text1+ ' intensity' if i == 0 else "", color='yellowgreen', alpha=0.6, marker = 1, s=40)

  df_all = hd.load_df_table_maxQ(prot1, not(LFQ))
  df = df_all[df_all.index.isin(indexes)]
  df_remove = df_all[~df_all.index.isin(indexes)]
  df_remove['max_intensity'] = df_remove[['Rep1', 'Rep2', 'Rep3']].max(axis=1)
  df_remove = df_remove.sort_values(by = 'max_intensity', ascending = False)
  for i, rep in enumerate(['Rep1_log10', 'Rep2_log10', 'Rep3_log10']):
    ax.scatter(df.index, df[rep], label="Protein test with "+title_text2+" intensity" if i == 0 else "", color='magenta', alpha=0.3, marker = 'o', s=40)
  list_others = ['others1', 'others2', 'others3']
  for i, rep in enumerate(['Rep1_log10', 'Rep2_log10', 'Rep3_log10']):
    ax.scatter(list_others*(df_remove.shape[0]//3)+list_others[0:df_remove.shape[0]%3], df_remove[rep], color='green', label="Protein only identified with "+title_text2+ " intensity for "+str(df_remove.shape[0])+' other genes' if i == 0 else "", alpha=0.3, marker = 'o', s=40)
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

def plot_log10_abundance(prot1, threshold, LFQ):
  '''Plot log10(emPAI) value for each gene for controls and test.'''
  fig,ax = plt.subplots()
  width = 12.5 ; height = 6.35 # taille finale de ta figure png
#  width = 14.4 ; height = 7.15 # taille finale de ta figure svg
  fig.set_size_inches(width, height)
#  var_test = hd.load_df_equal_test()
  df = hd.load_df_table_maxQ(prot1, LFQ)
  indexes = df.index
  df['max_intensity'] = df[['Rep1', 'Rep2', 'Rep3']].max(axis=1)
  minval = df[['Rep1', 'Rep2', 'Rep3']].min().min() # min value of the table
  maxval = df[['Rep1', 'Rep2', 'Rep3']].max().max() # max value of the table
  log_minval = np.log10(minval)
  df = df.sort_values(by = 'max_intensity', ascending = False)
  if LFQ == True:
    title_text1 = 'LFQ'
  else:
    title_text1 = 'raw'
  for i,rep in enumerate(['Rep1_log10', 'Rep2_log10', 'Rep3_log10']):
    ax.scatter(df.index, df[rep], label="Protein test with "+title_text1+' intensity' if i == 0 else "", color='royalblue', alpha=0.5, marker = 'o', s=40)
  plt.title(prot1[0]+' in '+prot1[1]+' with intensity = '+str(threshold)+' as threshold ')
  for i,rep in enumerate(['CtrA1_log10', 'CtrA2_log10', 'CtrA3_log10']):
    ax.scatter(df.index, df[rep], label="CtrA (without SPA tag) with "+title_text1+ ' intensity' if i == 0 else "", color='red', alpha=0.6, marker = 0, s=40)
  if 'CtrA4' in df :
    ax.scatter(df.index, np.log10(df.CtrA4), color='red', alpha=0.6, marker = 0, s=40)
  for i,rep in enumerate(['CtrC1_log10', 'CtrC2_log10', 'CtrC3_log10']):
    ax.scatter(df.index, df[rep], label="CtrC (with SPA tag) with "+title_text1+ ' intensity' if i == 0 else "", color='yellowgreen', alpha=0.6, marker = 1, s=40)

  dftrue = df[df.C_is == True] 
  abs_ctrC =[] # triangle if absent in control
  abs_ctrA = []
  sigA = []
  sigC = [] # get index list for stars in plot
  prot_sig_C = []
  prot_sig_A = [] # get list of proteins where test is significant. 
  for i, my_index in enumerate(df.index):
    prot_sig_C.append(my_index)
    if (df.CtrC1.iloc[i] == threshold) and (df.CtrC2.iloc[i] == threshold) and (df.CtrC3.iloc[i] == threshold):
      abs_ctrC.append(i)
    elif df.C_is.iloc[i] == True:
      sigC.append(i)
    else: 
      prot_sig_C.pop()

    prot_sig_A.append(my_index)
    if 'CtrA4' in df :
      if (df.CtrA1.iloc[i] == threshold and df.CtrA2.iloc[i] == threshold and df.CtrA3.iloc[i] == threshold and df.CtrA4.iloc[i] == threshold):
        abs_ctrA.append(i)
      elif df.A_is.iloc[i] == True:
        sigA.append(i)
      else: 
        prot_sig_A.pop()
    else:
      if (df.CtrA1.iloc[i] == threshold and df.CtrA2.iloc[i] == threshold and df.CtrA3.iloc[i] == threshold):
        abs_ctrA.append(i)
      elif df.A_is.iloc[i] == True:
        sigA.append(i)
      else: 
        prot_sig_A.pop()

  abs_ctrC = [x+0.2 for x in abs_ctrC];  abs_ctrA = [x-0.2 for x in abs_ctrA]
  sigC = [x+0.2 for x in sigC];  sigA = [x-0.2 for x in sigA]
  prot_sig = list(set(prot_sig_A) | set(prot_sig_C)) # list of names of proteins.
  prot_sig.sort()
  print('cont ?', prot_sig)
  print('abs ctr ?', abs_ctrA, abs_ctrC)
  print('all signif :', len(prot_sig))
  c = 0; nc = 0
  for i in prot_sig:
    if i in contaminant_genes:
      c += 1
#      print(i, end = ', ')
  print()
  print('contaminant :', c)
  for i in prot_sig:
    if i not in contaminant_genes:
      nc += 1
#      print(i, end =', ')
  print()
  print('not contaminant :', nc)

  ax.scatter(sigA, [log_minval-0.2]*len(sigA),c='red', marker=(5, 2), label = 'Significant test with CtrA', s=30) # add stars for Significant controls.
  ax.scatter(sigC, [log_minval-0.2]*len(sigC),c='yellowgreen', marker=(5, 2), label = 'Significant test with CtrC', s=30) # add stars for Significant controls.

  ax.scatter(abs_ctrA, [log_minval-0.2]*len(abs_ctrA),c='red', marker='^', label = 'Protein absent of CtrA', s = 20) # add triangles for proteins absent of each replicate of control A.
  ax.scatter(abs_ctrC, [log_minval-0.2]*len(abs_ctrC),c='yellowgreen', marker='^', label = 'Protein absent of CtrC', s= 20) # add triangles for proteins absent of each replicate of control C.

# Confidence interval
  # df_rep = np.log2(df[['Rep1', 'Rep2', 'Rep3']])
  # mean_conf_int = st.mean_confidence_interval(df_rep, 0.95, st.get_global_variance(prot1, threshold))
  # mean_conf_int = mean_conf_int.reindex(index = df.index)
  # print(mean_conf_int)  
  # ax.plot( mean_conf_int['mean'], '-', linewidth=1, color = 'royalblue', alpha = 0.5)
  # ax.fill_between(mean_conf_int.index, mean_conf_int['conf_inf'], mean_conf_int['conf_sup'], color='royalblue', alpha=.15)

  # df_ctr = np.log2(df[['CtrC1', 'CtrC2', 'CtrC3']])
  # mean_conf_int = st.mean_confidence_interval(df_ctr, 0.95, st.get_global_variance(prot1, threshold))
  # mean_conf_int = mean_conf_int.reindex(index = df.index)
  # print(mean_conf_int)
  # ax.plot( mean_conf_int['mean'], '-', linewidth=1, color = 'yellowgreen', alpha = 0.5)
  # ax.fill_between(mean_conf_int.index, mean_conf_int['conf_inf'], mean_conf_int['conf_sup'], color='yellowgreen', alpha=.2)

  # if 'CtrA4' in df:
  #   df_ctr = np.log2(df[['CtrA1', 'CtrA2', 'CtrA3', 'CtrA4']])
  # else:
  #   df_ctr = np.log2(df[['CtrA1', 'CtrA2', 'CtrA3']])
  # mean_conf_int = st.mean_confidence_interval(df_ctr, 0.95, st.get_global_variance(prot1, threshold))
  # mean_conf_int = mean_conf_int.reindex(index = df.index)
  # print(mean_conf_int)
  # ax.plot( mean_conf_int['mean'], '-', linewidth=1, color = 'red', alpha = 0.3)
  # ax.fill_between(mean_conf_int.index, mean_conf_int['conf_inf'], mean_conf_int['conf_sup'], color='red', alpha=.1)

  fig.canvas.draw()
  plt.ylim(log_minval-0.3, np.log10(maxval)+0.3)
  plt.xlabel('Gene name')
  plt.ylabel('log10(intensity) value')
  plt.grid(axis = 'x') # vertical lines
  plt.xticks(rotation=90)

  for ticklabel in ax.get_xticklabels(): # adjust legend.
    if ticklabel.get_text() in contaminant_genes: # contaminant_genes est la liste de genes contaminants
      ticklabel.set_color('orange')
    if ticklabel.get_text() in prey[prot1[0]]:
      ticklabel.set_color('blue')
    if ticklabel.get_text() in interesting_prey[prot1[0]]:
      ticklabel.set_color('cyan')
  handles, labels = ax.get_legend_handles_labels()
  cont_patch = mpatches.Patch(color='orange', label='potential contaminant proteins')
  certain_interactor_patch = mpatches.Patch(color='blue', label='confirmed interactor proteins')
  interesting_interactor_patch = mpatches.Patch(color='cyan', label='interesting interactor proteins')
  handles.extend([cont_patch, certain_interactor_patch, interesting_interactor_patch]) # add to legend
  plt.legend(handles=handles)

  path_batch = "maxQ/Images/my_plots/"
# get an appropriate plot and saved image.
  manager = plt.get_current_fig_manager() # get full screen
  manager.window.showMaximized() # get full screen
  fig.tight_layout()
  fig.subplots_adjust(left=.05, bottom=.2, right=.96, top=.93) # marges
  filename = path_batch+prot1[0]+'_'+prot1[1][:6].replace('/', '_')+'_'+title_text1+'_'+str(threshold)+'_log10values.png'
#  plt.savefig(path_batch+'test.svg') # image vectorisée
  plt.savefig(filename, transparent=False, dpi = 300) # image pixelisée, dpi = résolution
#  plt.show()

df = load_df_maxQ('F2')
df_int = select_intensity(df, False)
df_LFQ = select_intensity(df, True)

print('intensity', df_int.shape)
print('LFQ', df_LFQ.shape)
print('identif type', df[df['Identification_type_sample'].notnull()].shape)
#print(df.iloc[11])
#bn = dt.get_batches(pd_samples, 'DnaA', 'LB log') #batch names
#replicates = [load_df_maxQ(bn[0]), load_df_maxQ(bn[1]), load_df_maxQ(bn[2])] # pandas replicates
used_prot_tuples = dp.good_proteins()
contaminant_genes = dt.load_contaminant_list()
threshold = 1000000
#for prot in [used_prot_tuples[2], used_prot_tuples[3], used_prot_tuples[5]]:
prot_3_rep = prot_three_rep()
# print(prot_two_rep())
for prot in prot_3_rep[0:3]:
  dt.header(prot[0]+ ' in '+prot[1])
  # create_table(prot, threshold, False)
  # create_table(prot, threshold, True)
  # log10_calculation(prot, False)
  # log10_calculation(prot, True)
  non_eq_var_ttest(prot, threshold, False)
  plot_log10_abundance(prot, threshold, LFQ = False)

# plot_var_repartition(prot_3_rep, threshold)
# test_bart_per_all_prots(prot_3_rep, threshold)
