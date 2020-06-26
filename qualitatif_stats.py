import Load_PPI_screen as dt
import Data_processing as dp
import csv
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import handle as hd
import os
import glob
from scipy import stats
import statsmodels.stats.multitest
from matplotlib.colors import LinearSegmentedColormap

pd_samples = dt.load_based_screen_samples()
pd_controls = dt.load_based_screen_controls()
used_prot_tuples = dt.good_proteins() 
controls_typeC = {'LB log':['S1','S2','S3'], 'LB O/N':['O9', 'R4', 'R5'], 'M9 0.2% ac O/N':['R1', 'R2', 'R3']}
controls_typeA = {'LB log':['L1', 'T7', 'T8', 'T9'], 'LB O/N':['C13', 'P5', 'U10'], 'M9 0.2% ac O/N':['A10', 'T5', 'T6']}

def log_calculation(prot1, data):
  '''Calculation of emPAI logs and for a specific protein/condition.'''
  path_batch = "MS data csv reportbuilder-20200408T212019Z-001/Used_proteins_csv/"
  df = hd.load_df_table(prot1, True)
  if 'CtrC4' in df:
    df = df[['Rep1', 'Rep2', 'Rep3', 'CtrC1', 'CtrC2', 'CtrC3', 'CtrA1','CtrA2', 'CtrA3', 'CtrA4']]
  else:
    df = df[['Rep1', 'Rep2', 'Rep3', 'CtrC1', 'CtrC2', 'CtrC3', 'CtrA1','CtrA2', 'CtrA3']]
  for i in df.columns:
    newname = i+'_log2'
    df[newname] = np.log2(df[i])
  df.to_csv(path_batch+ prot1[0]+"_"+prot1[1].replace('/', '_')+'_multipleRep.csv')

def wrong_bartlett_test(all_prots, data, threshold):
  '''Statistical test that checks if variance is similar or not between predatory proteins of a similar bait protein. H0 is equal variance.
  data can be 0 (test) or 1(ctrA) or 2(ctrC).'''
  fig, ax = plt.subplots(3,3)
  for ax1 in range(3):
    for ax2 in range(3):
      prot = all_prots[3*ax1+ax2]
      if data == 0: #test
        data_type = 'test'
        df = hd.load_df_table(prot, True)[['Rep1_log2', 'Rep2_log2', 'Rep3_log2']]
      elif data == 1: #ctrA
        data_type = 'ctrlA'
        df = hd.load_df_table(prot, False)
        if 'CtrA4' in df.columns:
          df = df[['CtrA1', 'CtrA2', 'CtrA3', 'CtrA4']]
        else:
          df = df[['CtrA1', 'CtrA2', 'CtrA3']]
        df = df.apply(lambda x: np.where(x < threshold ,threshold,x))
        df = np.log2(df)
      elif data == 2: #ctrC
        data_type = 'ctrlC'
        df = hd.load_df_table(prot, True)[['CtrC1_log2', 'CtrC2_log2', 'CtrC3_log2']]
      else : print('error in data variable')
      b_test = []
      print(df.shape)
      for i in df.columns:
        df = df.loc[df[i] != np.log2(threshold)]
      print('Number of different proteins :', df.shape[0])
      m = np.zeros((df.shape[0], df.shape[0]))
      for i,idx1 in enumerate(df.index):
        for j,idx2 in enumerate(df.index):
          test = stats.bartlett(df.iloc[i], df.iloc[j])
          b_test.append(test[1])
          m[i, j] = test[1]
      thresh = 0.1
      nodes = [0, 0.00001,thresh, thresh, 1.0]
      colors = ['grey',"yellow",'yellow', 'palegreen', 'blue']
      cmap = LinearSegmentedColormap.from_list("", list(zip(nodes, colors)))
    #  cmap.set_under("gray")
      im = ax[ax1,ax2].imshow(m, cmap=cmap, vmin=0, vmax=1)

      m[m == 0] = 2
      m = np.triu(m,k=1)
    # Counts how many tests are significant or not.
      print('Number of different proteins :', df.shape[0])
      print("Not equal var :", np.sum((m <= thresh) & (m>0)))
      print("Equal var :", np.sum((m > thresh) & (m<2)))
      print("Error test :", np.sum(m == 2))
      ax[ax1,ax2].set_title(prot[0]+' in '+prot[1])
  fig.subplots_adjust(right=0.8)
  cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
  cbar = fig.colorbar(im, cax = cbar_ax, extend="min")
  cbar.set_label('pvalue', rotation=270, labelpad=15)
  fig.suptitle('Variance between prey proteins for '+data_type+ ' replicates')
  plt.show()
  return b_test

#bartlett_test(used_prot_tuples,0, 0.25)


def get_samples_bartlett_test(prot, data, threshold):
  '''Statistical test that checks if variance is similar or not between predatory proteins of a similar bait protein. H0 is equal variance.
  data can be 0 (test) or 1(ctrA) or 2(ctrC).
  Returns a list that contains a list of log2(emPAI) values for each protein.'''
  if data == 0: #test
    data_type = 'test'
    df = hd.load_df_table(prot, False)
    df = df[['Rep1', 'Rep2', 'Rep3']]
  elif data == 1: #ctrA
    data_type = 'ctrlA'
    df = hd.load_df_table(prot, False)
    if 'CtrA4' in df.columns:
      df = df[['CtrA1', 'CtrA2', 'CtrA3', 'CtrA4']]
    else:
      df = df[['CtrA1', 'CtrA2', 'CtrA3']]
  elif data == 2: #ctrC
    data_type = 'ctrlC'
    df = hd.load_df_table(prot, False)
    df = df[['CtrC1', 'CtrC2', 'CtrC3']]
  else : print('error in data variable')
  df = df.apply(lambda x: np.where(x < threshold ,threshold,x))
  for i in df.columns:
    df = df.loc[df[i] != threshold]
  df = np.log2(df)
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
  my_array = [list(df.loc[i]) for i in df.index]
  return my_array

def test_bart(my_array):
  '''Realization of the test.'''
  return stats.bartlett(*my_array)

def per_all_prots():
  '''All prots from all condition with test and controls.'''
  for i in used_prot_tuples:
    dt.header(i[0]+' in '+i[1])
    array = []
    array += get_samples_bartlett_test(i,0,0.25)
    array += get_samples_bartlett_test(i,1,0.25)
    array += get_samples_bartlett_test(i,2,0.25)
    print('pval25 :', test_bart(array)[1])
    array = []
    array += get_samples_bartlett_test(i,0,0.04)
    array += get_samples_bartlett_test(i,1,0.04)
    array += get_samples_bartlett_test(i,2,0.04)
    print('pval04 :', test_bart(array)[1])

#per_all_prots()

def per_cond():
  log = []
  ON = []
  M9 = []
  for i in used_prot_tuples:
    if i[1] == 'LB log':
      log += get_samples_bartlett_test(i, 0, 0.25)
    elif i[1] == 'LB O/N':
      ON += get_samples_bartlett_test(i, 0, 0.25)
    elif i[1] == 'M9 0.2% ac O/N':
      M9 += get_samples_bartlett_test(i, 0, 0.25)
  print('log', test_bart(log)[0])
  print('ON', test_bart(ON)[1])
  print('M9', test_bart(M9)[1])

def per_prot():
  '''Proteins that are present in different growth conditions.'''
  dnaA = []
  seqa = []
  nrdb = []
  for i in used_prot_tuples:
    if i[0] == 'DnaA':
      dnaA += get_samples_bartlett_test(i, 0, 0.25)
    elif i[0] == 'SeqA':
      seqa += get_samples_bartlett_test(i, 0, 0.25)
    elif i[0] == 'NrdB':
      nrdb += get_samples_bartlett_test(i, 0, 0.25)
  print('DnaA', test_bart(dnaA)[1])
  print('SeqA', test_bart(seqa)[1])
  print('NrdB', test_bart(nrdb)[1])

def per_type():
  rep = []
  ctra = []
  ctrc = []
  for i in used_prot_tuples:
  #  dt.header(i[0]+i[1])
    rep += get_samples_bartlett_test(i,0,0.25)
    ctra += get_samples_bartlett_test(i,1,0.25)
    ctrc += get_samples_bartlett_test(i,2,0.25)
  print('rep', test_bart(rep)[1])
  print('ctrA', test_bart(ctra)[1])
  print('ctrC', test_bart(ctrc)[1])
#  samp1 = get_samples_bartlett_test(i, 0, 0.25)
#  samp2 = get_samples_bartlett_test(i, 1, 0.25)
#  samp3 = get_samples_bartlett_test(i, 2, 0.25)
#  print(test_bart(samp3)[1])
#  print('a', test_bart(a)[1])

#per_cond()
#per_prot()
#per_type()

def compare_tests():
  rep = []
  rep += get_samples_bartlett_test(used_prot_tuples[1],0,0.25)
  print(len(rep))
  print(*rep)
  print('bartlett', stats.bartlett(*rep)[1])
  print('levene', stats.levene(*rep)[1])
  print('fligner', stats.fligner(*rep)[1])

#compare_tests()
common_var = hd.load_df_equal_test()
#print(common_var)

df = hd.load_df_table(used_prot_tuples[0], True)[['Rep1_log2', 'Rep2_log2', 'Rep3_log2']]
#print(df)

def mean_confidence_interval(data, confidence, variance = 0):
    n = data.shape[1]
    print(n)
    m = df.mean(axis = 1) # mean
    if variance == 0:
      se = data.sem(axis = 1) # standard error(erreur type = s/racine(n))      
    else:
      se = np.sqrt(variance/n)
    h = se * stats.t.ppf((1 + confidence) / 2., n-1) # t value such that (1-alpha = 0.95)
    print(stats.t.ppf((0.8), n-1))
#    return stats.t.interval(0.95, len(a)-1, loc=np.mean(a), scale=stats.sem(a))
    d = {'mean': m,'conf_inf': m-h, 'conf_sup': m+h}
    return pd.DataFrame(d)

#res = mean_confidence_interval(df, 0.95)
#print(res)
#print(res.shape)

def get_global_variance(prot, threshold):
  '''Returns the global variance for test and ctrC and ctrA.'''
  all_vars = []
  for data in [0,1,2]:
    if data == 0: #test
      data_type = 'test'
      df = hd.load_df_table(prot, False)
      df = df[['Rep1', 'Rep2', 'Rep3']]
    elif data == 1: #ctrA
      data_type = 'ctrlA'
      df = hd.load_df_table(prot, False)
      if 'CtrA4' in df.columns:
        df = df[['CtrA1', 'CtrA2', 'CtrA3', 'CtrA4']]
      else:
        df = df[['CtrA1', 'CtrA2', 'CtrA3']]
    elif data == 2: #ctrC
      data_type = 'ctrlC'
      df = hd.load_df_table(prot, False)
      df = df[['CtrC1', 'CtrC2', 'CtrC3']]
    else : print('error in data variable')
    df = df.apply(lambda x: np.where(x < threshold ,threshold,x))
    for i in df.columns:
      df = df.loc[df[i] != threshold]
    df = np.log2(df)
    nul_var_idx = []
    for i, my_index in enumerate(df.index): # remove row when var = 0.
      nul_var = [df.iloc[i,j] == df.iloc[i,j+1] for j in range(len(df.columns)-1)]
      if nul_var.count(True) == len(nul_var):
        nul_var_idx.append(my_index)
    df = df.drop(nul_var_idx)
    print('Number of different proteins :', df.shape[0])
    df_norm = df.sub(df.mean(axis=1), axis=0)
    all_vars += list(df_norm.var(axis=1))
  return np.mean(all_vars)

def test_normal_equal_var(prot1, threshold):
  '''Test de z'''
  df = hd.load_df_table(prot1, True)
  glob_var = get_global_variance(prot1, threshold)
  print('global_var :', glob_var)
  df['glob_var'] = glob_var
#  print(df)
  df_rep = df[['Rep1_log2', 'Rep2_log2', 'Rep3_log2']]
  if 'CtrA4' in df:
    df_ctrA = df[['CtrA1_log2', 'CtrA2_log2', 'CtrA3_log2', 'CtrA4_log2']]
  else:
    df_ctrA = df[['CtrA1_log2', 'CtrA2_log2', 'CtrA3_log2']]
  df_ctrC = df[['CtrC1_log2', 'CtrC2_log2', 'CtrC3_log2']]
  mRep = df_rep.mean(axis = 1)
  mCtrA = df_ctrA.mean(axis = 1)
  mCtrC = df_ctrC.mean(axis = 1)
  
  df['zvalueA'] = (mRep - mCtrA)/np.sqrt((df['glob_var']/len(df_rep.columns)+df['glob_var']/len(df_ctrA.columns)))
  df['zvalueC'] = (mRep - mCtrC)/np.sqrt((df['glob_var']/len(df_rep.columns)+df['glob_var']/len(df_ctrC.columns)))
  df['pvalA'] = stats.norm.sf(df['zvalueA'])
#  one sided test with the survival function of student = 1-Fct de rep°.
  df['pvalC'] = stats.norm.sf(df['zvalueC'])
  sigA = df['pvalA'][df['pvalA']<0.05]
#  print( df['zvalueA'])
#  print(len(sigA))
#  print(len(df['zvalueA']))
  pval_C_corr = statsmodels.stats.multitest.multipletests(df['pvalC'], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
  pval_A_corr = statsmodels.stats.multitest.multipletests(df['pvalA'], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
#  print(testCdone)
  df['pvalA_corr'] = pval_A_corr[1]
  df['A_is'] = pval_A_corr[0]
  df['pvalC_corr'] = pval_C_corr[1]
  df['C_is'] = pval_C_corr[0]
#  print(df['pvalC'])
  path_batch = "MS data csv reportbuilder-20200408T212019Z-001/Used_proteins_csv/"
  df.to_csv(path_batch+ prot1[0]+"_"+prot1[1].replace('/', '_')+'_multipleRep.csv')


#for i in used_prot_tuples[0:1]:
#  dt.header(i[0]+ 'in '+i[1])
#  test_normal_equal_var(i, 0.25)


def corrected_ttest(prot1, threshold):
  '''Welch test = (mA-mB)/sqrt(sA^2/nA+nB^2/nB) against both controls, one at a time. Variance is estimated for the whole file, it is a global variance.'''
  df = hd.load_df_table(prot1, True)
  is_fourth = 0
  if 'CtrA4' in df:
    is_fourth = 1
  df_rep = df[['Rep1_log2', 'Rep2_log2', 'Rep3_log2']]
  if is_fourth == 0:
    df_ctrA = df[['CtrA1_log2', 'CtrA2_log2', 'CtrA3_log2']]
  else:
    df_ctrA = df[['CtrA1_log2', 'CtrA2_log2', 'CtrA3_log2', 'CtrA4_log2']]
  df_ctrC = df[['CtrC1_log2', 'CtrC2_log2', 'CtrC3_log2']]
  mRep = df_rep.mean(axis = 1)
  mCtrA = df_ctrA.mean(axis = 1)
  mCtrC = df_ctrC.mean(axis = 1)
  log_threshold = np.log2(threshold)
  df_varrep = df_rep.loc[(df_rep['Rep1_log2'] != log_threshold) & (df_rep['Rep2_log2'] != log_threshold) & (df_rep['Rep3_log2'] != log_threshold)] # keep only when 3 replicates contain the protein. 
  df_varctrA = df_ctrA
  for i in df_ctrA.columns:
    df_varctrA = df_varctrA.loc[df_varctrA[i] != log_threshold]
  df_varctrC = df_ctrC
  for i in df_ctrC.columns:
    df_varctrC = df_varctrC.loc[df_varctrC[i] != log_threshold]
  df['var_Rep'] = np.mean(df_varrep.var(axis = 1)) # same value for the whole file
  df['var_ctrA'] = np.mean(df_varctrA.var(axis = 1))
  df['var_ctrC'] = np.mean(df_varctrC.var(axis = 1))
  print('var_rep', df['var_Rep'][0])
  print('varC', df['var_ctrC'][0])
  print('varA', df['var_ctrA'][0])
  df['tvalueA'] = (mRep - mCtrA)/np.sqrt((df['var_Rep']/len(df_rep.columns)+df['var_ctrA']/len(df_ctrA.columns)))
  df['tvalueC'] = (mRep - mCtrC)/np.sqrt((df['var_Rep']/len(df_rep.columns)+df['var_ctrC']/len(df_ctrC.columns)))
#  print(df['tvalueA'])
  df['pvalA'] = stats.t.sf(df['tvalueA'], df= len(df_rep.columns)+len(df_ctrA.columns)-2)
#  one sided test with the survival function of student = 1-Fct de rep°.
  df['pvalC'] = stats.t.sf(df['tvalueC'], df= len(df_rep.columns)+len(df_ctrC.columns)-2)
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
  pval_C_corr = statsmodels.stats.multitest.multipletests(df['pvalC'], alpha=0.1, method='fdr_bh', is_sorted=False, returnsorted=False)
  pval_A_corr = statsmodels.stats.multitest.multipletests(df['pvalA'], alpha=0.1, method='fdr_bh', is_sorted=False, returnsorted=False)
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
  path_batch = "MS data csv reportbuilder-20200408T212019Z-001/Used_proteins_csv/"
  df.to_csv(path_batch+ prot1[0]+"_"+prot1[1].replace('/', '_')+'_multipleRep.csv')

#for i in used_prot_tuples:
#  print(i[0], i[1])
#  corrected_ttest(i, 0.25)

def set_pval(used_prot_tuples):
  '''Basic t-tests with independant samples, not equal var and estimation of variance at each test. '''
  for prot1 in used_prot_tuples:
    print(prot1[0], prot1[1])
    df = hd.load_df_table(prot1, True)
    if 'Rep1_log2' not in df.columns:
      for i in df.columns:
        newname = i+'_log2'
        df[newname] = np.log2(df[i])
    else : print('ok')
    testC = stats.ttest_ind(df[['Rep1_log2', 'Rep2_log2', 'Rep3_log2']], df[['CtrC1_log2', 'CtrC2_log2', 'CtrC3_log2']], equal_var = False, axis = 1)
    real_pval_testC = testC[1]/2 # two-tailed test so you have to divide the obtained pvalue by 2.
    if 'CtrA4' in df:
      testA = stats.ttest_ind(df[['Rep1_log2', 'Rep2_log2', 'Rep3_log2']], df[['CtrA1_log2', 'CtrA2_log2', 'CtrA3_log2', 'CtrA4_log2']], equal_var = False, axis = 1) # two-tailed test so you have to divide the obtained pvalue by 2. 
    else :
      testA = stats.ttest_ind(df[['Rep1_log2', 'Rep2_log2', 'Rep3_log2']], df[['CtrA1_log2', 'CtrA2_log2', 'CtrA3_log2']], equal_var = False, axis = 1) 
    for i in range(len(testA[0])): # here, it is a two-tailed test but we want one tailed so we have to divide the obtained pvalue by 2 when statistic > 0.
      if testA[0][i] < 0:
        testA[1][i] = 1-testA[1][i]/2
      else : testA[1][i] = testA[1][i]/2
      if testC[0][i] < 0:
        testC[1][i] = 1-testC[1][i]/2
      else : testC[1][i] = testC[1][i]/2

    pval_C_corr = statsmodels.stats.multitest.multipletests(testC[1], alpha=0.1, method='fdr_bh', is_sorted=False, returnsorted=False)
    pval_A_corr = statsmodels.stats.multitest.multipletests(testA[1], alpha=0.1, method='fdr_bh', is_sorted=False, returnsorted=False)
#    print(pval_A_corr)
#    print(pval_C_corr)
    df['tstat_C'] = testC[0]
    df['pval_C'] = testC[1]
    df['pval_C_corr'] = pval_C_corr[1]
    df['C_is'] = pval_C_corr[0]
    df['tstat_A'] = testA[0]
    df['pval_A'] = testA[1]
    df['pval_A_corr'] = pval_A_corr[1]
    df['A_is'] = pval_A_corr[0]
    def pval_is (row, ctrl):
      if row[ctrl] <=0.05 :
        return True
      else: return False 
    df['C_is'] = df.apply (lambda row: pval_is(row, 'pval_C'), axis=1)
    df['A_is'] = df.apply (lambda row: pval_is(row, 'pval_A'), axis=1)
    path_batch = "MS data csv reportbuilder-20200408T212019Z-001/Used_proteins_csv/"
    df.to_csv(path_batch+ prot1[0]+"_"+prot1[1].replace('/', '_')+'_multipleRep.csv')

