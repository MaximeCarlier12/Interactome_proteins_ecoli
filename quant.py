from globvar import *
import Load_PPI_screen as dt
import Data_processing as dp
import qualitatif_stats as st
import handle as hd
import glob
import statsmodels.stats.multitest
import matplotlib.patches as mpatches

plt.rcParams.update(params)

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

def load_df_maxQ(bfNumber, norm):
  '''From a batch file number, it gets the right corresponding dataframe. 
  norm = 0 :without norm, 1 : med norm and 2 : bait value norm, 3 : Q1 value, 4: Q3 value.'''
  path_batch = "maxQ/SEGREGATED-20200619T092017Z-001/SEGREGATED/single files/"
  # full_path get the only good file and can be in two different ways (containing 'A_5' or 'A5').
  if len(bfNumber) == 1:
    letter = ''
    number = bfNumber
    full_path = path_batch+"1/"+number
  else :
    letter = bfNumber[0]
    number = bfNumber[1:]
    full_path = path_batch+letter+'/'+letter+number
  if norm == 1:
      df = pd.read_csv(full_path+'_Norm_Med.csv', header = 0)
  elif norm ==2:
      df = pd.read_csv(full_path+'_Norm_Bait.csv', header = 0)
  elif norm ==3:
    df = pd.read_csv(full_path+'_Norm_Q1.csv', header = 0)
  elif norm ==4:
    df = pd.read_csv(full_path+'_Norm_Q3.csv', header = 0)
  else:
    df = pd.read_excel(full_path+".xlsx", sheet_name = 0, header = 0)
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

def normalize_med_file(bn):
  '''Create a new batch file with normalized raw data'''
  df = load_df_maxQ(bn, 0)
  df = select_intensity(df, False)
  med = df['Intensity_sample'].median(axis = 0)
  print('med', med)
  df['Normalized_intensity'] = df['Intensity_sample']/med
  path_batch = "maxQ/SEGREGATED-20200619T092017Z-001/SEGREGATED/single files/"
  if len(bn) == 1:
    letter = ''
    number = bn
    full_path = path_batch+"1/"
  else :
    letter = bn[0]
    number = bn[1:]
    full_path = path_batch+letter+'/'
  df.to_csv(full_path+bn+'_Norm_Med.csv')
  return med

def normalize_q1_file(bn):
  '''Create a new batch file with normalized raw data'''
  df = load_df_maxQ(bn, 0)
  df = select_intensity(df, False)
  q1 = df['Intensity_sample'].quantile(0.25)
  print('q1', q1)
  df['Normalized_intensity'] = df['Intensity_sample']/q1
  path_batch = "maxQ/SEGREGATED-20200619T092017Z-001/SEGREGATED/single files/"
  if len(bn) == 1:
    letter = ''
    number = bn
    full_path = path_batch+"1/"
  else :
    letter = bn[0]
    number = bn[1:]
    full_path = path_batch+letter+'/'
  df.to_csv(full_path+bn+'_Norm_Q1.csv')
  return q1

def normalize_q3_file(bn):
  '''Create a new batch file with normalized raw data'''
  df = load_df_maxQ(bn, 0)
  df = select_intensity(df, False)
  q3 = df['Intensity_sample'].quantile(0.75)
  print('q3', q3)
  df['Normalized_intensity'] = df['Intensity_sample']/q3
  path_batch = "maxQ/SEGREGATED-20200619T092017Z-001/SEGREGATED/single files/"
  if len(bn) == 1:
    letter = ''
    number = bn
    full_path = path_batch+"1/"
  else :
    letter = bn[0]
    number = bn[1:]
    full_path = path_batch+letter+'/'
  df.to_csv(full_path+bn+'_Norm_Q3.csv')
  return q3

def normalize_bait_batch(bn, bait_name):
  '''Bait_name is a gene name that is from PROTEINS variable.'''
  df = load_df_maxQ(bn, 0)
  df = select_intensity(df, False)
  bait_name = bait_name[0].lower()+bait_name[1:]
  if bait_name in df.index:
    bait_value = df.loc[bait_name, 'Intensity_sample']
    print('prot name : '+ bait_name+', prot_value : '+ str(bait_value))
    df['Normalized_intensity'] = df['Intensity_sample']/bait_value
    path_batch = "maxQ/SEGREGATED-20200619T092017Z-001/SEGREGATED/single files/"
    if len(bn) == 1:
      letter = ''
      number = bn
      full_path = path_batch+"1/"
    else :
      letter = bn[0]
      number = bn[1:]
      full_path = path_batch+letter+'/'
    df.to_csv(full_path+bn+'_Norm_Bait.csv')
  else:
    print('not found for :'+ bn)

def create_all_norm_files(prot_3_rep):
  '''For each condition [protein, growth cond], it creates 3 files for the 3 different normalizations (median, q1, q3).'''
  med_values = []
  q3_values = []
  for prot in prot_3_rep:
    dt.header(prot[0]+' in '+prot[1])
    prot_batches = dt.get_batches(pd_samples, prot[0], prot[1])
    prot_batches = get_batches_missing_files(prot_batches)
    for bn in prot_batches :
      print('batch :', bn)
      med_values.append(normalize_med_file(bn))
      normalize_q1_file(bn)
      q3_values.append(normalize_q3_file(bn))
      # normalize_bait_batch(bn, prot[0])
  for cond in CONDITION:
    batchA = dp.get_control_batches(pd_controls, 'MG1655 (TYPE A)' , cond)
    batchA = get_batches_missing_files(batchA)
    for bn in batchA :
      print('batch :', bn)
      med_values.append(normalize_med_file(bn))
      normalize_q1_file(bn)
      q3_values.append(normalize_q3_file(bn))
      # normalize_bait_batch(bn, 'rplB')
    batchC = dp.get_control_batches(pd_controls, 'MG1655 (placI)mVenus-SPA-pUC19 (TYPE C2)' , cond)
    batchC = get_batches_missing_files(batchC)
    for bn in batchC:
      print('batch :', bn)
      med_values.append(normalize_med_file(bn))
      normalize_q1_file(bn)
      q3_values.append(normalize_q3_file(bn))
      # normalize_bait_batch(bn, 'rplB')
  return (med_values, q3_values)

def log10_calculation_new(df):
  '''Calculation of intensity logs10 and for a specific protein/condition, before the creatipn of the table.'''
  for i in df.columns:
    newname = i+'_log10' 
    df[newname] = np.log10(df[i])
  return df

def link_ctr_rep(prot, LFQ, normalize):
  '''Cou t how many proteins are present in controls, among the ones present in the 3 replicates of the tests.'''
  prot_batches = dt.get_batches(pd_samples, prot[0], prot[1])
  prot_batches = get_batches_missing_files(prot_batches)
  rep = [select_intensity(load_df_maxQ(i, normalize), LFQ) for i in prot_batches]
  ctrC = [select_intensity(load_df_maxQ(i, normalize), LFQ) for i in get_batches_missing_files(controls_typeC[prot[1]])]
  ctrA = [select_intensity(load_df_maxQ(i, normalize), LFQ) for i in get_batches_missing_files(controls_typeA[prot[1]])]
  inter_rep = dp.intersect_index(rep)
  list_nb_ctr_per_prot = []
  for i in inter_rep:
    counter = 0
    for batch in ctrC:
      if i in batch.index.values:
        counter += 1
    list_nb_ctr_per_prot.append(counter)
  print('nb prot 3 rep :', len(inter_rep))
  print('nb prot 3 rep & 0 ctrC:', list_nb_ctr_per_prot.count(0))
  print('nb prot 3 rep & 1 ctrC:', list_nb_ctr_per_prot.count(1))
  print('nb prot 3 rep & 2 ctrC:', list_nb_ctr_per_prot.count(2))
  print('nb prot 3 rep & 3 ctrC:', list_nb_ctr_per_prot.count(3))
  list_nb_ctr_per_prot = []
  for i in inter_rep:
    counter = 0
    for batch in ctrA:
      if i in batch.index.values:
        counter += 1
    list_nb_ctr_per_prot.append(counter)
  print('nb prot 3 rep :', len(inter_rep))
  print('nb prot 3 rep & 0 ctrA:', list_nb_ctr_per_prot.count(0))
  print('nb prot 3 rep & 1 ctrA:', list_nb_ctr_per_prot.count(1))
  print('nb prot 3 rep & 2 ctrA:', list_nb_ctr_per_prot.count(2))
  print('nb prot 3 rep & 3 ctrA:', list_nb_ctr_per_prot.count(3))
  print('nb prot 3 rep & 3 ctrA:', list_nb_ctr_per_prot.count(4))

def create_table(prot, LFQ, normalize):
  '''Create a file containing intensity values for each gene present in our 3 test replicates. Absence equals to threshold.'''
  prot_batches = dt.get_batches(pd_samples, prot[0], prot[1])
  prot_batches = get_batches_missing_files(prot_batches)
  rep = [select_intensity(load_df_maxQ(i, normalize), LFQ) for i in prot_batches]
  ctrC = [select_intensity(load_df_maxQ(i, normalize), LFQ) for i in get_batches_missing_files(controls_typeC[prot[1]])]
  ctrA = [select_intensity(load_df_maxQ(i, normalize), LFQ) for i in get_batches_missing_files(controls_typeA[prot[1]])]
  print('nb prot intersect', len(dp.intersect_index(rep)))
  indexes = list(dp.intersect_index(rep))
  print('nb genes:', len(indexes))
  feature_list = ['Rep1', 'Rep2', 'Rep3', 'CtrC1', 'CtrC2', 'CtrC3', 'CtrA1', 'CtrA2', 'CtrA3']
  if len(ctrA) == 4:
    feature_list.append('CtrA4')
  df = pd.DataFrame(0, index=indexes, columns=feature_list)
  path_batch = "maxQ/SEGREGATED-20200619T092017Z-001/Protein_table/"
  for (i,repi) in enumerate(rep):
    for my_index, row in repi.iterrows():
      if my_index in indexes:
        # print(my_index)
        if LFQ == True:
          df.loc[my_index, feature_list[i]] = row.LFQ_intensity_sample
        elif normalize == 0:
          df.loc[my_index, feature_list[i]] = row.Intensity_sample
        else:
          df.loc[my_index, feature_list[i]] = row.Normalized_intensity
  for (i,repi) in enumerate(ctrC):
    repi = select_intensity(repi, LFQ)
    for my_index, row in repi.iterrows():
      if my_index in indexes:
        if LFQ == True:
          df.loc[my_index, feature_list[3+i]] = row.LFQ_intensity_sample
        elif normalize == 0:
          df.loc[my_index, feature_list[3+i]] = row.Intensity_sample
        else:
          df.loc[my_index, feature_list[3+i]] = row.Normalized_intensity
  for (i,repi) in enumerate(ctrA):
    repi = select_intensity(repi, LFQ)
    for my_index, row in repi.iterrows():
      if my_index in indexes:
        if LFQ == True:
          df.loc[my_index, feature_list[6+i]] = row.LFQ_intensity_sample
        elif normalize ==0:
          df.loc[my_index, feature_list[6+i]] = row.Intensity_sample
        else:
          df.loc[my_index, feature_list[6+i]] = row.Normalized_intensity
  if normalize != 0:
    threshold = 0.1
  else:
    threshold = 100000
  print(threshold)
  # threshold = df[df > .0000001].min().min() # min value of the table
  df = df.apply(lambda x: np.where(x < threshold ,threshold,x))
  for col in df[['Rep1', 'Rep2', 'Rep3']].columns:
    df = df.loc[df[col] != threshold] # remove each row where one replicate value equals to threshold
  if normalize != 0:
    df = log10_calculation_new(df)
  else:
    df = log10_calculation_new(df)
  cond = prot[1][:6].replace('/', '_')
  if LFQ == True:
    df.to_csv(path_batch+ prot[0]+"_"+cond+'_LFQ.csv')
  elif normalize == 0:
    df.to_csv(path_batch+ prot[0]+"_"+cond+'_Int.csv')
  elif normalize == 1:
    df.to_csv(path_batch+ prot[0]+"_"+cond+'_Int_Norm_Med.csv')
  elif normalize == 2:
    df.to_csv(path_batch+ prot[0]+"_"+cond+'_Int_Norm_Bait.csv')
  elif normalize == 3:
    df.to_csv(path_batch+ prot[0]+"_"+cond+'_Int_Norm_Q1.csv')
  elif normalize == 4:
    df.to_csv(path_batch+ prot[0]+"_"+cond+'_Int_Norm_Q3.csv')
  return threshold

def save_table_update(df, prot, LFQ, normalize):
  path_batch = "maxQ/SEGREGATED-20200619T092017Z-001/Protein_table/"
  if LFQ == True:
    df.to_csv(path_batch+ prot[0]+"_"+prot[1][:6].replace('/', '_')+'_LFQ.csv')
  elif normalize == 0:
    df.to_csv(path_batch+ prot[0]+"_"+prot[1][:6].replace('/', '_')+'_Int.csv')
  elif normalize == 1:
    df.to_csv(path_batch+ prot[0]+"_"+prot[1][:6].replace('/', '_')+'_Int_Norm_Med.csv')
  elif normalize == 2:
    df.to_csv(path_batch+ prot[0]+"_"+prot[1][:6].replace('/', '_')+'_Int_Norm_Bait.csv')
  elif normalize == 3:
    df.to_csv(path_batch+ prot[0]+"_"+prot[1][:6].replace('/', '_')+'_Int_Norm_Q1.csv')
  elif normalize == 4:
    df.to_csv(path_batch+ prot[0]+"_"+prot[1][:6].replace('/', '_')+'_Int_Norm_Q3.csv')
  else :
    print('error : filename not possible to establish')

def log10_calculation(prot, LFQ, normalize):
  '''Calculation of intensity logs10 and for a specific protein/condition.'''
  df = hd.load_df_table_maxQ(prot, LFQ, normalize)
  if 'CtrC4' in df:
    df = df[['Rep1', 'Rep2', 'Rep3', 'CtrC1', 'CtrC2', 'CtrC3', 'CtrA1','CtrA2', 'CtrA3', 'CtrA4']]
  else:
    df = df[['Rep1', 'Rep2', 'Rep3', 'CtrC1', 'CtrC2', 'CtrC3', 'CtrA1','CtrA2', 'CtrA3']]
  for i in df.columns:
    newname = i+'_log10'
    df[newname] = np.log10(df[i])
  save_table_update(df, prot, LFQ, normalize)

def get_df_to_variance(prot, data, threshold, LFQ, normalize):
  '''Statistical test that checks if variance is similar or not between predatory proteins of a similar bait protein. H0 is equal variance.
  data can be 0 (test) or 1(ctrA) or 2(ctrC).
  Returns a list that contains a list of log10(intensity) values for each protein.'''
  if data == 0: #test
    data_type = 'test'
    df = hd.load_df_table_maxQ(prot, LFQ, normalize)
    df = df[['Rep1', 'Rep2', 'Rep3']]
  elif data == 1: #ctrA
    data_type = 'ctrlA'
    df = hd.load_df_table_maxQ(prot, LFQ, normalize)
    if 'CtrA4' in df.columns:
      df = df[['CtrA1', 'CtrA2', 'CtrA3', 'CtrA4']]
    else:
      df = df[['CtrA1', 'CtrA2', 'CtrA3']]
  elif data == 2: #ctrC
    data_type = 'ctrlC'
    df = hd.load_df_table_maxQ(prot, LFQ, normalize)
    df = df[['CtrC1', 'CtrC2', 'CtrC3']]
  else : print('error in data variable')
  # threshold = df[df > .0001].min().min() # min value of the table
  # df = df.apply(lambda x: np.where(x < threshold ,threshold,x))
  for i in df.columns:
    df = df.loc[df[i] > threshold]
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
  return df

def get_samples_bartlett_test(prot, data, threshold, LFQ, normalize):
  df = get_df_to_variance(prot, data, threshold, LFQ, normalize)
  my_array = [list(df.loc[i]) for i in df.index] #lists of 3 values 
  return my_array

def test_bart(my_array):
  '''Realisation of the test.'''
  return stats.levene(*my_array)

def test_bart_per_all_prots(used_prot_tuples, threshold, LFQ, normalize):
  '''All prots from all condition with test and controls.'''
  for prot in used_prot_tuples:
    dt.header(prot[0]+' in '+prot[1])
    array = []
    array += get_samples_bartlett_test(prot,0,threshold, LFQ, normalize)
    array += get_samples_bartlett_test(prot,1,threshold, LFQ, normalize)
    array += get_samples_bartlett_test(prot,2,threshold, LFQ, normalize)
    # if prot[0] == 'DnaG':
      # print(array)
    # print('Bartlett test for test and ctrA and ctrC :', test_bart(array)[1])
    a = get_samples_bartlett_test(prot,0, threshold, LFQ, normalize)
    print('Bartlett test for norm int in test replicates :', test_bart(a)[1])
    print('Bartlett test for ctrA :', test_bart(get_samples_bartlett_test(prot,1,threshold, LFQ, normalize))[1])
    print('Bartlett test for ctrC :', test_bart(get_samples_bartlett_test(prot,2,threshold, LFQ, normalize))[1])
    print('Global bartlett test :', test_bart(array)[1])

def plot_var_repartition(used_prot_tuples, threshold, LFQ, normalize):
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
    array += get_samples_bartlett_test(prot,0,threshold, LFQ, normalize)
    array += get_samples_bartlett_test(prot,1,threshold, LFQ, normalize)
    array += get_samples_bartlett_test(prot,2,threshold, LFQ, normalize)
    list_var = []
    for j in array:
      list_var.append(np.var(j))
    ax = fig.add_subplot(3,3,1+i-step)
    ax.set_title(prot[0]+' in '+prot[1])
    # ax.scatter(labx, np.sort(list_var))
    ax.hist(list_var, bins = 20)
  plt.subplots_adjust(hspace=0.5, wspace=1.0)
  plt.show()

def plot_rep_proteins_per_batch(used_prot_tuples, LFQ, normalize):
  '''Shows variance repartition for each protein to see if a normal distribution suits well.'''
  data = []
  for prot in used_prot_tuples:
    name = prot[0]+'_'+prot[1][:6]
    print(name)
    prot_batches = dt.get_batches(pd_samples, prot[0], prot[1])
    prot_batches = get_batches_missing_files(prot_batches)
    for bn in prot_batches:
      print(bn)
      df = load_df_maxQ(bn, normalize) # med norm
      med_coeff = normalize_med_file(bn)
      q3_coeff = normalize_q3_file(bn)
      if normalize != 0:
        split_mylist = pd.cut(df['Normalized_intensity'], [0, .5, 1, 5, 10, 100, 500, 2000, 5000, 10000]).value_counts(sort = False)
      else:
        df = select_intensity(df, LFQ)
        split_mylist = pd.cut(df['Intensity_sample'], [1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11]).value_counts(sort = False)
      split_mylist = split_mylist.tolist()
      print([name, bn,len(df.index)]+split_mylist+[med_coeff, q3_coeff])
      data.append([name, bn,len(df.index)]+split_mylist+[med_coeff, q3_coeff])
  CONDITION = ['LB log', 'LB O/N' ,'M9 0.2% ac O/N']
  for cond in CONDITION:
    batchA = dp.get_control_batches(pd_controls, 'MG1655 (TYPE A)' , cond)
    batchA = get_batches_missing_files(batchA)
    batchC = dp.get_control_batches(pd_controls, 'MG1655 (placI)mVenus-SPA-pUC19 (TYPE C2)' , cond)
    batchC = get_batches_missing_files(batchC)
    name = 'ctrA_'+cond[:6]
    for bn in batchA:
      print(bn)
      df = load_df_maxQ(bn, normalize) # med norm
      med_coeff = normalize_med_file(bn)
      q3_coeff = normalize_q3_file(bn)
      if normalize != 0:
        split_mylist = pd.cut(df['Normalized_intensity'], [0, .5, 1, 5, 10, 100, 500, 2000, 5000, 10000]).value_counts(sort = False)
      else:
        df = select_intensity(df, LFQ)
        split_mylist = pd.cut(df['Intensity_sample'], [1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11]).value_counts(sort = False)
      split_mylist = split_mylist.tolist()
      print([name, bn,len(df.index)]+split_mylist+[med_coeff, q3_coeff])
      data.append([name, bn,len(df.index)]+split_mylist+[med_coeff, q3_coeff])
    name = 'ctrC_'+cond[:6]
    for bn in batchC:
      print(bn)
      df = load_df_maxQ(bn, normalize) # med norm or q3 norm
      med_coeff = normalize_med_file(bn)
      q3_coeff = normalize_q3_file(bn)
      if normalize != 0:
        split_mylist = pd.cut(df['Normalized_intensity'], [0, .5, 1, 5, 10, 100, 500, 2000, 5000, 10000]).value_counts(sort = False)
      else:
        df = select_intensity(df, LFQ)
        split_mylist = pd.cut(df['Intensity_sample'], [1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10, 1e11]).value_counts(sort = False)
      split_mylist = split_mylist.tolist()
      print([name, bn,len(df.index)]+split_mylist+[med_coeff, q3_coeff])
      data.append([name, bn,len(df.index)]+split_mylist+[med_coeff, q3_coeff])
  if normalize == 1:
    df = pd.DataFrame(data, columns = ['Prot name', 'Batch name', 'Total', '[0-0.5]', '[0.5-1]', '[1-5]', '[5-10]', '[10-100]', '[100-500]', '[500-2000]', '[2000-5000]', '[5000-10000]', 'Normalized median coeff', 'Normalized q3 coeff']) 
    df.to_csv('maxQ/'+ 'values_distribution_per_batch_med.csv', index = False)
  elif normalize == 4:
    df = pd.DataFrame(data, columns = ['Prot name', 'Batch name', 'Total', '[0-0.5]', '[0.5-1]', '[1-5]', '[5-10]', '[10-100]', '[100-500]', '[500-2000]', '[2000-5000]', '[5000-10000]', 'Normalized median coeff', 'Normalized q3 coeff']) 
    df.to_csv('maxQ/'+ 'values_distribution_per_batch_q3.csv', index = False)
  elif normalize == 0:
    df = pd.DataFrame(data, columns = ['Prot name', 'Batch name', 'Total', '[1e4-1e5]', '[1e5-1e6]', '[1e6-1e7]', '[1e7-1e8]', '[1e8-1e9]', '[1e9-1e10]', '[1e10-1e11]', 'Normalized median coeff', 'Normalized q3 coeff']) 
    df.to_csv('maxQ/'+ 'values_distribution_per_batch_raw.csv', index = False)
  else : 
    print('error with normalize value')

def non_eq_var_ttest(prot, threshold, LFQ, normalize):
  '''Welch test = (mA-mB)/sqrt(sA^2/nA+sB^2/nB) against both controls, one at a time. We only take the variance between our 3 replicates.'''
  df = hd.load_df_table_maxQ(prot, LFQ, normalize)
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
  df['tvalueC'] = (mRep - mCtrC)/np.sqrt((df_rep.var(axis=1)/len(df_rep.columns)+df_ctrC.var(axis = 1)/len(df_ctrC.columns))) # df.var is with unbiased N-1 by default.
  #  print(df['tvalueA'])
  varrep = df_rep.var(axis=1)
  varA = df_ctrA.var(axis = 1)
  nA = len(df_ctrA.columns)
  varC = df_ctrC.var(axis = 1)
  # dfA = np.power(varA/nA+varrep/3, 2)/(np.power(varrep,4)/(3*3*2)+np.power(varA,4)/(nA*nA*(nA-1)))
  dfA = np.power(varA/nA+varrep/3, 2)/(varrep/(3*2)+varA/(nA*(nA-1)))
  print(dfA)
  # df['pvalA'] = stats.t.sf(df['tvalueA'], df= len(df_rep.columns)+len(df_ctrA.columns)-2)
  df['pvalA'] = stats.t.sf(df['tvalueA'], df= len(df_rep.columns)+len(df_ctrA.columns)-2)
  #  one sided test with the survival function of student = 1-Fct de rep°.
  df['pvalC'] = stats.t.sf(df['tvalueC'], df= len(df_rep.columns)+len(df_ctrC.columns)-2)
  print('dof :', len(df_rep.columns)+len(df_ctrC.columns)-2)
  # print('ex test student :', stats.t.sf(1.20727, 4))
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
  adj_pval_C = statsmodels.stats.multitest.multipletests(df['pvalC'], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
  adj_pval_A = statsmodels.stats.multitest.multipletests(df['pvalA'], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
  #  print(testCdone)
  df['adj_pvalC'] = adj_pval_C[1]
  df['C_is'] = adj_pval_C[0]
  df['adj_pvalA'] = adj_pval_A[1]
  df['A_is'] = adj_pval_A[0]
  def pval_is (row, ctrl):
    if row[ctrl] <=0.05 :
      return True
    else: return False 
  #  df['C_is'] = df.apply (lambda row: pval_is(row, 'pvalC'), axis=1)
  #  df['A_is'] = df.apply (lambda row: pval_is(row, 'pvalA'), axis=1)
  save_table_update(df, prot, LFQ, normalize)

def prot_absent_controls(prot, threshold, LFQ, normalize):
  '''On the table file, it adds a column 'Absent_ctrA' and 'Absent_ctrC' that tells if a protein is at least in one replicate of this control or not. '''
  df = hd.load_df_table_maxQ(prot, LFQ, normalize)
  abs_ctrC =[] # triangle if absent in control
  abs_ctrA = [] 
  for i, my_index in enumerate(df.index):
    if (df.CtrC1.iloc[i] == threshold) and (df.CtrC2.iloc[i] == threshold) and (df.CtrC3.iloc[i] == threshold):
      abs_ctrC.append(i)
    if 'CtrA4' in df :
      if (df.CtrA1.iloc[i] == threshold and df.CtrA2.iloc[i] == threshold and df.CtrA3.iloc[i] == threshold and df.CtrA4.iloc[i] == threshold):
        abs_ctrA.append(i)
    else:
      if (df.CtrA1.iloc[i] == threshold and df.CtrA2.iloc[i] == threshold and df.CtrA3.iloc[i] == threshold):
        abs_ctrA.append(i)

  for i, my_index in enumerate(df.index) : 
    if i in abs_ctrA:
      df.loc[my_index, 'Absent_ctrA'] = True
    else : 
      df.loc[my_index, 'Absent_ctrA'] = False
    if i in abs_ctrC:
      df.loc[my_index, 'Absent_ctrC'] = True
    else : 
      df.loc[my_index, 'Absent_ctrC'] = False
  save_table_update(df, prot, LFQ, normalize)

def get_global_variance_per_prot(prot, threshold, LFQ, normalize):
  '''Returns the global variance for test and ctrC and ctrA for each protein/grotwh condition.'''
  all_vars = []
  flat_list = []
  for data in [0,1,2]:
    if data == 0: #test
      data_type = 'test'
      df = hd.load_df_table_maxQ(prot, LFQ, normalize)
      df = df[['Rep1', 'Rep2', 'Rep3']]
    elif data == 1: #ctrA
      data_type = 'ctrlA'
      df = hd.load_df_table_maxQ(prot, LFQ, normalize)
      if 'CtrA4' in df.columns:
        df = df[['CtrA1', 'CtrA2', 'CtrA3', 'CtrA4']]
      else:
        df = df[['CtrA1', 'CtrA2', 'CtrA3']]
    elif data == 2: #ctrC
      data_type = 'ctrlC'
      df = hd.load_df_table_maxQ(prot, LFQ, normalize)
      df = df[['CtrC1', 'CtrC2', 'CtrC3']]
    else : print('error in data variable')
    # df = df.apply(lambda x: np.where(x < threshold ,threshold,x))
    for i in df.columns:
      df = df.loc[df[i] != threshold]
    df = np.log10(df)
    nul_var_idx = []
    for i, my_index in enumerate(df.index): # remove row when var = 0.
      nul_var = [df.iloc[i,j] == df.iloc[i,j+1] for j in range(len(df.columns)-1)]
      if nul_var.count(True) == len(nul_var):
        nul_var_idx.append(my_index)
    df = df.drop(nul_var_idx)
    print('Number of different proteins :', df.shape[0])
    df_norm = df.sub(df.mean(axis=1), axis=0)
    all_vars += list(df_norm.var(axis=1))
    list_pd = df_norm.values.tolist()
    for sublist in list_pd:
      for item in sublist:
        flat_list.append(item)
  print('biaised_estimated var', np.var(flat_list))
  print('unbiased_estimated var', np.var(flat_list, ddof = 1))
  print('over_estimated var', np.mean(all_vars))
  # flat_list = flat_list/np.std(flat_list)
  import statsmodels.api as sm 
  shapiro_test = stats.kstest(flat_list, 'norm')
  print('shapiro test :', shapiro_test)
  # plt.hist(flat_list, bins = 100)
  # plt.draw()
  # plt.show()
  # plt.close()
  return np.var(flat_list, ddof = 1)

def get_common_variances_per_media(used_prot_tuples, threshold, LFQ, normalize):
  '''Calculation of a variance per growth condition (only 3 different variances). This variance will be prefered to '''
  three_common_var = {}
  for gw_cond in CONDITION: # 3 common variances
    flat_list = []
    list_vars = []
    ctrl_added = False
    for prot in used_prot_tuples:
      if prot[1] == gw_cond: # protein to be consiered
        if ctrl_added == False: #other prot
          data_types = [0, 1, 2]
          ctrl_added = True
        else: # my_prot : we add control replicates as well
          data_types = [0]
        for data in data_types:      
          if data == 0: #replicates
            data_type = 'replicates'
            df = hd.load_df_table_maxQ(prot, LFQ, normalize)
            df = df[['Rep1', 'Rep2', 'Rep3']]
          elif data == 1: #ctrA
            data_type = 'ctrlA'
            df = hd.load_df_table_maxQ(prot, LFQ, normalize)
            if 'CtrA4' in df.columns:
              df = df[['CtrA1', 'CtrA2', 'CtrA3', 'CtrA4']]
            else:
              df = df[['CtrA1', 'CtrA2', 'CtrA3']]
          elif data == 2: #ctrC
            data_type = 'ctrlC'
            df = hd.load_df_table_maxQ(prot, LFQ, normalize)
            df = df[['CtrC1', 'CtrC2', 'CtrC3']]
          else : print('error in data variable')
          # df = df.apply(lambda x: np.where(x < threshold ,threshold,x)) # not necessary
          for i in df.columns:
            df = df.loc[df[i] != threshold]
          df = np.log10(df)
          nul_var_idx = []
          for i, my_index in enumerate(df.index): # remove row when var = 0.
            nul_var = [df.iloc[i,j] == df.iloc[i,j+1] for j in range(len(df.columns)-1)]
            if nul_var.count(True) == len(nul_var):
              nul_var_idx.append(my_index)
          df = df.drop(nul_var_idx)
          print('Number of different proteins :', df.shape[0])
          df_norm = df.sub(df.mean(axis=1), axis=0)
          list_vars += list(df_norm.var(axis=1))
          list_pd = df_norm.values.tolist()
          for sublist in list_pd:
            for item in sublist:
              flat_list.append(item)
    three_common_var[gw_cond] = np.var(flat_list, ddof = 1)
  print('shapiro flat list', stats.shapiro(flat_list))
  print('shapiro vars', stats.shapiro(list_vars))
  plt.hist(flat_list, bins = 100)
  plt.show()
  return three_common_var

def test_normal_equal_var(prot, threshold, LFQ, normalize, common_variance):
  '''Test de z'''
  df = hd.load_df_table_maxQ(prot, LFQ, normalize)
  # glob_var = get_global_variance_per_prot(prot, threshold, LFQ, normalize)
  glob_var = common_variance[prot[1]]
  print('global_var :', glob_var)
  df['glob_var'] = glob_var
  #  print(df)
  df_rep = df[['Rep1_log10', 'Rep2_log10', 'Rep3_log10']]
  if 'CtrA4' in df:
    df_ctrA = df[['CtrA1_log10', 'CtrA2_log10', 'CtrA3_log10', 'CtrA4_log10']]
  else:
    df_ctrA = df[['CtrA1_log10', 'CtrA2_log10', 'CtrA3_log10']]
  df_ctrC = df[['CtrC1_log10', 'CtrC2_log10', 'CtrC3_log10']]
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
  adj_pval_C = statsmodels.stats.multitest.multipletests(df['pvalC'][df['Absent_ctrC'] ==  False], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
  adj_pval_A = statsmodels.stats.multitest.multipletests(df['pvalA'][df['Absent_ctrA'] ==  False], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
  #  print(testCdone)
  iterA = 0; iterC = 0
  for i, my_index in enumerate(df.index):
    if df.loc[my_index, 'Absent_ctrA'] == True: # this protein is absent in control
      df.loc[my_index, 'adj_pvalA'] = df.loc[my_index, 'pvalA']
      df.loc[my_index, 'A_is'] = True
    else: # this protein is present in controls
      df.loc[my_index,'adj_pvalA'] = adj_pval_A[1][iterA]
      df.loc[my_index, 'A_is'] = adj_pval_A[0][iterA]
      iterA += 1
    if df.loc[my_index, 'Absent_ctrC'] == True:
      df.loc[my_index, 'adj_pvalC'] = df.loc[my_index, 'pvalC']
      df.loc[my_index, 'C_is'] = True
    else:
      df.loc[my_index,'adj_pvalC'] = adj_pval_C[1][iterC]
      df.loc[my_index, 'C_is'] = adj_pval_C[0][iterC]
      iterC += 1
  #  print(df['pvalC'])
  save_table_update(df, prot, LFQ, normalize)

def fold_change(prot, LFQ, normalize):
  '''Calculate log2 of fold change for each protein. Creates 2 columns fo each control : FC and log2FC. Keep only log2(FC)>1.5 or >0.6 ? It means FC = 2.8 or FC>2'''
  df = hd.load_df_table_maxQ(prot, LFQ, normalize)
  print('mean', df[['Rep1', 'Rep2', 'Rep3']].mean(axis = 1))
  if 'CtrA4' in df.columns:
    print('diff_mean', df[['Rep1', 'Rep2', 'Rep3']].mean(axis = 1)/df[['CtrA1', 'CtrA2', 'CtrA3', 'CtrA4']].mean(axis = 1))
    df['FC_A'] = df[['Rep1', 'Rep2', 'Rep3']].mean(axis = 1)/df[['CtrA1', 'CtrA2', 'CtrA3', 'CtrA4']].mean(axis = 1)
  else:
    df['FC_A'] = df[['Rep1', 'Rep2', 'Rep3']].mean(axis = 1)/df[['CtrA1', 'CtrA2', 'CtrA3']].mean(axis =1)
  df['FC_C'] = df[['Rep1', 'Rep2', 'Rep3']].mean(axis = 1)/df[['CtrC1', 'CtrC2', 'CtrC3']].mean(axis =1)
  df['log2FC_A'] = np.log2(df['FC_A'])
  df['log2FC_C'] = np.log2(df['FC_C'])
  save_table_update(df, prot, LFQ, normalize)

def max_int_replicates(prot, LFQ, normalize):
  '''Know which replicate has highest values.'''
  df = hd.load_df_table_maxQ(prot, LFQ, normalize)
  df = df[['Rep1_log10', 'Rep2_log10', 'Rep3_log10']]
  df['max_intensity'] = df.max(axis=1)
  conditions1 = [
    df['Rep1_log10'] == df['max_intensity'],
    df['Rep2_log10'] == df['max_intensity'],
    df['Rep3_log10'] == df['max_intensity']] # we know the first one
  conditions3 = [
  (df['Rep1_log10'] < df['Rep2_log10']) & (df['Rep1_log10'] < df['Rep3_log10']),
  (df['Rep2_log10'] <= df['Rep1_log10']) & (df['Rep2_log10'] < df['Rep3_log10']),
  (df['Rep3_log10'] <= df['Rep1_log10']) & (df['Rep3_log10'] <= df['Rep2_log10']) ] # we know the last one
  choices = [1,2,3]
  df['First_rep'] = np.select(conditions1, choices, default = 0)
  df['Third_rep'] = np.select(conditions3, choices, default = 0)
  df['Sec_rep'] = 6-df['First_rep']- df['Third_rep']
  print('norm =', normalize)
  print('nb genes', len(df.index))
  dt.header('first')
  print(df['First_rep'].value_counts())
  dt.header('second')
  print(df['Sec_rep'].value_counts())
  dt.header('third')
  print(df['Third_rep'].value_counts())

def plot_log10_abundance_comparison(prot, threshold, LFQ, normalize, contaminant_genes, control = 'AC'):
  '''Plot log10(emPAI) value for each gene for controls and test. Compare LFQ and raw intensity.'''
  fig,ax = plt.subplots()
  width = 12.5 ; height = 6.35 # taille finale de ta figure png
  #  width = 14.4 ; height = 7.15 # taille finale de ta figure svg
  fig.set_size_inches(width, height)
  #  var_test = hd.load_df_equal_test()
  df = hd.load_df_table_maxQ(prot, LFQ, normalize)
  indexes = df.index
  df['max_intensity'] = df[['Rep1', 'Rep2', 'Rep3']].max(axis=1)
  minval = df[['Rep1', 'Rep2', 'Rep3']].min().min() # min value of the table
  maxval = df[['Rep1', 'Rep2', 'Rep3']].max().max() # max value of the table
  df = df.sort_values(by = 'max_intensity', ascending = False)
  if LFQ == True:
    title_text1 = 'LFQ'
    title_text2 = 'raw'
  elif normalize == 0:
    title_text1 = 'raw'
  else:
    title_text1 = 'normalize_raw'
  for i,rep in enumerate(['Rep1_log10', 'Rep2_log10', 'Rep3_log10']):
    ax.scatter(df.index, df[rep], label="Protein test with "+title_text1+' intensity' if i == 0 else "", color='royalblue', alpha=0.5, marker = 'o', s=40)
  plt.title(prot[0]+' in '+prot[1]+' with intensity = '+str(threshold)+' as threshold ')
  if 'A' in control:
    for i,rep in enumerate(['CtrA1_log10', 'CtrA2_log10', 'CtrA3_log10']):
      ax.scatter(df.index, df[rep], label="CtrA (without SPA tag) with "+title_text1+ ' intensity' if i == 0 else "", color='red', alpha=0.6, marker = 0, s=40)
    if 'CtrA4' in df :
      ax.scatter(df.index, np.log10(df.CtrA4), color='red', alpha=0.6, marker = 0, s=40)
  if 'C' in control:
    for i,rep in enumerate(['CtrC1_log10', 'CtrC2_log10', 'CtrC3_log10']):
      ax.scatter(df.index, df[rep], label="CtrC (with SPA tag) with "+title_text1+ ' intensity' if i == 0 else "", color='forestgreen', alpha=0.6, marker = 1, s=40)

  df_all = hd.load_df_table_maxQ(prot, not(LFQ), normalize)
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
  filename = path_batch+prot[0]+'_'+prot[1][:6].replace('/', '_')+'_'+str(threshold)+'_log10values.png'
  #  plt.savefig(path_batch+'test.svg') # image vectorisée
  plt.savefig(filename, transparent=False, dpi = 300) # image pixelisée, dpi = résolution
  #  plt.show()

def plot_log10_abundance(prot, threshold, LFQ, normalize, common_variance):
  '''Plot log10(emPAI) value for each gene for controls and test.'''
  fig,ax = plt.subplots()
  width = 12.5 ; height = 6.35 # taille finale de ta figure png
  #  width = 14.4 ; height = 7.15 # taille finale de ta figure svg
  fig.set_size_inches(width, height)
  #  var_test = hd.load_df_equal_test()
  df = hd.load_df_table_maxQ(prot, LFQ, normalize)
  indexes = df.index
  df['mean_intensity'] = df[['Rep1_log10', 'Rep2_log10', 'Rep3_log10']].mean(axis=1)
  minval = df[['Rep1_log10', 'Rep2_log10', 'Rep3_log10']].min().min() # min value of the replicates. Reduces variability due to the absence of a value. 
  print('log10_minval_replicates', minval)
  minval = np.log10(threshold)
  print('log10_threshold', minval)
  maxval = df[['Rep1_log10', 'Rep2_log10', 'Rep3_log10', 'CtrC1_log10', 'CtrC2_log10', 'CtrC3_log10', 'CtrA1_log10', 'CtrA2_log10', 'CtrA3_log10']].max().max() # max value of the table
  print('maxval', maxval)
  df = df.sort_values(by = 'mean_intensity', ascending = False)
  # df = df.iloc[4:32, ]
  if LFQ == True:
    title_text1 = 'LFQ' # name of the file
  elif normalize == 0:
    title_text1 =  'raw'
  elif normalize == 1:
    title_text1 = 'norm_med_raw'
  elif normalize ==2:
    title_text1 = 'norm_bait_raw'
  elif normalize == 3:
    title_text1 = 'norm_q1_raw'
  elif normalize == 4:
    title_text1 = 'norm_q3_raw'
  # plt.title(prot[0]+' in '+prot[1]+' with intensity = '+str(threshold)+' as threshold ')
  plt.title(prot[0]+' in '+prot[1])
  for i,rep in enumerate(['Rep1_log10', 'Rep2_log10', 'Rep3_log10']):
    ax.scatter(df.index, df[rep], label="Test" if i == 0 else "", color='royalblue', alpha=0.5, marker = 'o', s=40)
  for i,rep in enumerate(['CtrA1_log10', 'CtrA2_log10', 'CtrA3_log10']):
    ax.scatter(df.index, df[rep], label="CtrA (without SPA tag)" if i == 0 else "", color='red', alpha=0.6, marker = 0, s=40)
  if 'CtrA4' in df :
    ax.scatter(df.index, np.log10(df.CtrA4), color='red', alpha=0.6, marker = 0, s=40)
  for i,rep in enumerate(['CtrC1_log10', 'CtrC2_log10', 'CtrC3_log10']):
    ax.scatter(df.index, df[rep], label="CtrC (with SPA tag)" if i == 0 else "", color='forestgreen', alpha=0.6, marker = 1, s=40)

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
  if normalize != 0: # if it has been normalized
    ax.scatter(sigA, [minval-0.1]*len(sigA),c='red', marker=(5, 2), label = 'Significant test with CtrA', s=30) # add stars for Significant controls.
    ax.scatter(sigC, [minval-0.1]*len(sigC),c='forestgreen', marker=(5, 2), label = 'Significant test with CtrC', s=30) # add stars for Significant controls.
    ax.scatter(abs_ctrA, [minval-0.1]*len(abs_ctrA),c='red', marker='^', label = 'Absent of CtrA', s = 20) # add triangles for proteins absent of each replicate of control A.
    ax.scatter(abs_ctrC, [minval-0.1]*len(abs_ctrC),c='forestgreen', marker='^', label = 'Absent of CtrC', s= 20) # add triangles for proteins absent of each replicate of control C.
  else:
    ax.scatter(sigA, [minval-0.2]*len(sigA),c='red', marker=(5, 2), label = 'Significant test with CtrA', s=30) # add stars for Significant controls.
    ax.scatter(sigC, [minval-0.2]*len(sigC),c='forestgreen', marker=(5, 2), label = 'Significant test with CtrC', s=30) # add stars for Significant controls.
    ax.scatter(abs_ctrA, [minval-0.2]*len(abs_ctrA),c='red', marker='^', label = 'Absent of CtrA', s = 20) # add triangles for proteins absent of each replicate of control A.
    ax.scatter(abs_ctrC, [minval-0.2]*len(abs_ctrC),c='forestgreen', marker='^', label = 'Absent of CtrC', s= 20) # add triangles for proteins absent of each replicate of control C.

  # Confidence interval

  df_rep = np.log10(df[['Rep1', 'Rep2', 'Rep3']])
  # mean_conf_int = st.mean_confidence_interval(df_rep, 0.95, get_global_variance_per_prot(prot, threshold, LFQ, normalize))
  mean_conf_int = st.mean_confidence_interval(df_rep, 0.95, common_variance[prot[1]])

  mean_conf_int = mean_conf_int.reindex(index = df.index)
  # print(mean_conf_int)  
  ax.plot( mean_conf_int['mean'], '--', linewidth=0.7, color = 'royalblue', alpha = 0.5)
  ax.fill_between(mean_conf_int.index, mean_conf_int['conf_inf'], mean_conf_int['conf_sup'], color='royalblue', alpha=.08)

  df_ctr = np.log10(df[['CtrC1', 'CtrC2', 'CtrC3']])
  # mean_conf_int = st.mean_confidence_interval(df_ctr, 0.95, get_global_variance_per_prot(prot, threshold, LFQ, normalize))
  mean_conf_int = st.mean_confidence_interval(df_ctr, 0.95, common_variance[prot[1]])

  mean_conf_int = mean_conf_int.reindex(index = df.index)
  # print(mean_conf_int)
  ax.plot( mean_conf_int['mean'], '--', linewidth=0.7, color = 'forestgreen', alpha = 0.5)
  ax.fill_between(mean_conf_int.index, mean_conf_int['conf_inf'], mean_conf_int['conf_sup'], color='forestgreen', alpha=.08)

  if 'CtrA4' in df:
    df_ctr = np.log10(df[['CtrA1', 'CtrA2', 'CtrA3', 'CtrA4']])
  else:
    df_ctr = np.log10(df[['CtrA1', 'CtrA2', 'CtrA3']])
  # mean_conf_int = st.mean_confidence_interval(df_ctr, 0.95, get_global_variance_per_prot(prot, threshold, LFQ, normalize))
  mean_conf_int = st.mean_confidence_interval(df_ctr, 0.95, common_variance[prot[1]])
  mean_conf_int = mean_conf_int.reindex(index = df.index)
  # print(mean_conf_int)
  ax.plot( mean_conf_int['mean'], '--', linewidth=0.7, color = 'red', alpha = 0.3)
  ax.fill_between(mean_conf_int.index, mean_conf_int['conf_inf'], mean_conf_int['conf_sup'], color='red', alpha=.06)

  fig.canvas.draw()
  if normalize != 0:
    plt.ylim(minval-0.2, maxval+0.1)
  else : 
    plt.ylim(minval-0.3, maxval+0.3)
  plt.xlabel('Gene name')
  plt.ylabel('log10(intensity) value')
  plt.grid(axis = 'x') # vertical lines
  plt.xticks(rotation=90)

  for ticklabel in ax.get_xticklabels(): # adjust legend.
    gene_name = ticklabel.get_text()
    if gene_name in contaminant_genes: # contaminant_genes est la liste de genes contaminants
      ticklabel.set_color('orange')
    if gene_name in prey[prot[0]]:
      ticklabel.set_color('blue')
    if gene_name in interesting_prey[prot[0]]:
      ticklabel.set_color('cyan')
  handles, labels = ax.get_legend_handles_labels()
  cont_patch = mpatches.Patch(color='orange', label='Potential contaminant')
  certain_interactor_patch = mpatches.Patch(color='blue', label='Confirmed interactor')
  interesting_interactor_patch = mpatches.Patch(color='cyan', label='Interesting interactor')
  handles.extend([cont_patch, certain_interactor_patch, interesting_interactor_patch]) # add to legend
  plt.legend(handles=handles, loc='upper right')

  path_batch = "maxQ/Images/my_plots/"
  # get an appropriate plot and saved image.
  manager = plt.get_current_fig_manager() # get full screen
  manager.window.showMaximized() # get full screen
  fig.tight_layout()
  fig.subplots_adjust(left=.05, bottom=.2, right=.96, top=.93) # marges
  filename = path_batch+prot[0]+'_'+prot[1][:6].replace('/', '_')+'_'+title_text1+'_log10values.png'
  #  plt.savefig(path_batch+'test.svg') # image vectorisée
  # plt.savefig(filename, transparent=False, dpi = 300) # image pixelisée, dpi = résolution
  plt.savefig('test_holD.png', transparent=False, dpi = 300) # image pixelisée, dpi = résolution
  # plt.show()

def create_entire_final_csv(prot, LFQ, normalize):
  '''Create final files that will be sent to the collaborators. It contains important information only on interesting proteins such as pvalues, FC, putative type, mean between replicates.'''
  df = hd.load_df_table_maxQ(prot, LFQ, normalize)
  df['Mean replicates log10'] = df[['Rep1_log10', 'Rep2_log10', 'Rep3_log10']].mean(axis = 1)
  df['Individual variance replicates log10'] = df[['Rep1_log10', 'Rep2_log10', 'Rep3_log10']].var(axis = 1) # default : ddof = 1
  if 'CtrA4' in df.columns:
    df['Mean ctrA log10'] = df[['CtrA1_log10', 'CtrA2_log10', 'CtrA3_log10', 'CtrA4_log10']].mean(axis = 1)
    df['Individual variance ctrA log10'] = df[['CtrA1_log10', 'CtrA2_log10', 'CtrA3_log10', 'CtrA4_log10']].var(axis = 1)
  else:
    df['Mean ctrA log10'] = df[['CtrA1_log10', 'CtrA2_log10', 'CtrA3_log10']].mean(axis = 1)
    df['Individual variance ctrA log10'] = df[['CtrA1_log10', 'CtrA2_log10', 'CtrA3_log10']].var(axis = 1)
  df['Mean ctrC log10'] = df[['CtrC1_log10', 'CtrC2_log10', 'CtrC3_log10']].mean(axis = 1)
  df['Individual variance ctrC log10'] = df[['CtrC1_log10', 'CtrC2_log10', 'CtrC3_log10']].var(axis = 1)
  df['Gene name'] = df.index
  # keep only interesting proteins (significant)
  print(df.shape)
  bn = dt.get_batches(pd_samples, prot[0], prot[1]) #batch names
  full_data = load_df_maxQ(bn[0], 1)
  for i, my_index in enumerate(df.index):
    df.loc[my_index, 'Protein names'] = full_data.loc[my_index, 'Protein names']
    if my_index in prey[prot[0]]:
      df.loc[my_index, 'Putative protein type'] = 'Well-confirmed prey'
    elif my_index in interesting_prey[prot[0]]:
      df.loc[my_index, 'Putative protein type'] = 'Interesting prey'
    elif my_index in contaminant_genes:
      df.loc[my_index, 'Putative protein type'] = 'Potential contaminant'
    else:
      df.loc[my_index, 'Initial characterization'] = 'None'
  df = df[['Gene name', 'Protein names', 'Mean replicates log10', 'Individual variance replicates log10', 'Mean ctrA log10', 'Individual variance ctrA log10', 'Mean ctrC log10', 'Individual variance ctrC log10', 'glob_var', 'adj_pvalA', 'adj_pvalC', 'log2FC_A', 'log2FC_C', 'Putative protein type']]
  df = df.sort_values(by=['adj_pvalA'])
  df = df.rename(columns={"glob_var": "Global variance", "adj_pvalA": "Adjusted pvalue with ctrA", "adj_pvalC": "Adjusted pvalue with ctrC", "log2FC_A" : "Log2 enrichment with ctrA", "log2FC_C" : "Log2 enrichment with ctrC"})
  comment = "* : a protein is absent from a certain control when the 'Mean ctr log10' value of this control equals -1 (and individual variance of the control equals 0). In that case, the corresponding adjusted pvalue is the pvalue."
  df.loc['last'] = ['INFORMATION']+[comment]+[None]*(len(df.columns)-2) # specify an information at the end of the file.
  # df.loc['INFORMATION'] = [comment]+[None]*len(df.columns-1) # specify an information at the end of the file.
  if LFQ == True:
    title_text1 = 'LFQ' # name of the file
  elif normalize == 0:
    title_text1 =  'raw'
  elif normalize == 1:
    title_text1 = 'med'
  elif normalize ==2:
    title_text1 = 'bait'
  elif normalize == 3:
    title_text1 = 'q1'
  elif normalize == 4:
    title_text1 = 'q3'
  path_batch = "maxQ/Final_results/"+prot[0]+'/'
  df.to_excel(path_batch+prot[0]+'_'+prot[1][:6].replace('/', '_')+'_'+title_text1+'_summary.xlsx', index = False)

def create_significant_prot_final_csv(prot, LFQ, normalize):
  '''Same function than create_entire_final_csv but only for enriched proteins (absent or significantly enriched compared to one control).'''
  path_batch = "maxQ/Final_results/"+prot[0]+'/'
  if LFQ == True:
    title_text1 = 'LFQ' # name of the file
  elif normalize == 0:
    title_text1 =  'raw'
  elif normalize == 1:
    title_text1 = 'med'
  elif normalize ==2:
    title_text1 = 'bait'
  elif normalize == 3:
    title_text1 = 'q1'
  elif normalize == 4:
    title_text1 = 'q3'
  else:
    print('error in normalize value')
  full_path = path_batch+ prot[0]+"_"+prot[1][:6].replace('/', '_')+'_'+title_text1+'_summary.xlsx'

  df_entire = pd.read_excel(full_path, sep=',', header=0, index_col = 0, skipfooter = 1)
  df = hd.load_df_table_maxQ(prot, LFQ, normalize)
  print(df_entire.shape)
  print(df.shape)

  df = df_entire[df['Absent_ctrA'] | df['Absent_ctrC']| ((df['adj_pvalA'] <= 0.05) & df['adj_pvalA']) | ((df['adj_pvalC'] <= 0.05) & df['adj_pvalC'])]
  print(df.shape)
  comment = "* : a protein is absent from a certain control when the 'Mean ctr log10' value of this control equals -1 (and individual variance of the control equals 0). In that case, the corresponding adjusted pvalue is the pvalue."

  df.to_excel(path_batch+prot[0]+'_'+prot[1][:6].replace('/', '_')+'_'+title_text1+'_significant_prots.xlsx', index = True)

def create_putative_proteins_analysis(used_prot_tuples, prot_name, LFQ, normalize):
  path_batch = "maxQ/Final_results/"+prot_name+'/'
  if LFQ == True:
    title_text1 = 'LFQ' # name of the file
  elif normalize == 0:
    title_text1 =  'raw'
  elif normalize == 1:
    title_text1 = 'med'
  elif normalize ==2:
    title_text1 = 'bait'
  elif normalize == 3:
    title_text1 = 'q1'
  elif normalize == 4:
    title_text1 = 'q3'
  else:
    print('error in normalize value')
  # filename = path_batch+prot[0]+'_'+title_text1+'_putative_proteins.csv'
  filename = path_batch+prot_name+'_'+title_text1+'_putative_proteins.xlsx'
  writer = pd.ExcelWriter(filename, engine='xlsxwriter')
  for cond in CONDITION: # for one bait protein, all the conditions are present in a same file
    if (prot_name, cond) in used_prot_tuples: # check if this condition has been studied (absence of replicates ?)
      full_path = path_batch+prot_name+"_"+cond[:6].replace('/', '_')+'_'+title_text1+'_significant_prots.xlsx'
      df_sig = pd.read_excel(full_path, sep=',', header=0, index_col = 0)
      df = hd.load_df_table_maxQ((prot_name, cond[:6].replace('/', '_')), LFQ, normalize)

      cont_NS = []; cont_S = []; interest_S = []; interest_NS = []; unk_S = []; unk_NS = []; conf_S = []; conf_NS = []
      for i in df.index:
        points = 0 
        if i in df_sig.index: # signif means odd, NS means even.
          points += 1
        if i in contaminant_genes:
          points += 10
        if i in interesting_prey[prot_name]:
          points += 20
        if i in prey[prot_name]:
          points += 30
        if points == 0:
          unk_NS.append(i)
        elif points == 1 :
          unk_S.append(i)
        elif points == 10 :
          cont_NS.append(i)
        elif points == 11 :
          cont_S.append(i)
        elif points == 20 :
          interest_NS.append(i)
        elif points == 21 :
          interest_S.append(i)
        elif points == 30 :
          conf_NS.append(i)
        elif points == 31 :
          conf_S.append(i)
        else : 'error : no class gene'
      cont_NS = ', '.join(cont_NS)
      cont_S = ', '.join(cont_S)
      interest_S = ', '.join(interest_S)
      interest_NS = ', '.join(interest_NS)
      unk_S = ', '.join(unk_S)
      unk_NS = ', '.join(unk_NS)
      conf_S = ', '.join(conf_S)
      conf_NS = ', '.join(conf_NS)
      data = {'Enriched': [cont_S, conf_S, interest_S, unk_S], 'Not enriched' : [cont_NS, conf_NS, interest_NS, unk_NS]}
      df_output = pd.DataFrame(data, index = ['Potential contaminant', 'Well-confirmed prey', 'Interesting prey', 'Unknown'])
      df_output.to_excel(writer, sheet_name=cond.replace('/', '_')[:6])
  writer.save()
  # with open(filename, 'w') as f:
  #   all_df[0].to_csv(f)
  # for i in range(1, len(all_df)):
  #   with open(filename, 'a') as f:
  #     all_df[i].to_csv(f)



df = load_df_maxQ('F2', 0)
df_int = select_intensity(df, LFQ = False)
df_LFQ = select_intensity(df, LFQ = True)

# print('intensity', df_int.shape)
# print('LFQ', df_LFQ.shape)
# print('identif type', df[df['Identification_type_sample'].notnull()].shape)
#print(df.iloc[11])
#bn = dt.get_batches(pd_samples, 'DnaA', 'LB log') #batch names
#replicates = [load_df_maxQ(bn[0], False), load_df_maxQ(bn[1], False), load_df_maxQ(bn[2], False)] # pandas replicates
used_prot_tuples = dp.good_proteins()
#for prot in [used_prot_tuples[2], used_prot_tuples[3], used_prot_tuples[5]]:

prot_3_rep = prot_three_rep()
# print(prot_two_rep())
threshold_med = 0.1
threshold_raw = 100000

# for prot in prot_3_rep:
#   threshold_med = create_table(prot, LFQ = False, normalize = 1)
#   prot_absent_controls(prot, threshold_med, False, 1)
#   threshold_q3 = create_table(prot, LFQ = False, normalize = 4)
#   prot_absent_controls(prot, threshold_med, False, 4)

# common_variances_med = get_common_variances_per_media(prot_3_rep, threshold_med, False, 1)
# common_variances_q3 = get_common_variances_per_media(prot_3_rep, threshold_med, False, 4)

# create_all_norm_files(prot_3_rep)
for prot in []:
  dt.header(prot[0]+ ' in '+prot[1])
    # link_ctr_rep(prot, False, 1)
  # threshold_raw = create_table(prot, False, 0)
#   # threshold_bait = create_table(prot, LFQ = False, normalize = 2)
  # threshold_q1 = create_table(prot, False, 3)
  # threshold_q3 = create_table(prot,0 False, 4)
# #   # log10_calculation(prot, False, False)
# #   # log10_calculation(prot, False, False)
#   # for norm in [1]:
#     # print('norm =', norm)
#     # max_int_replicates(prot, False, norm)
#   # non_eq_var_ttest(prot, threshold_med, LFQ = False, normalize = 1)

  # for i in [1,4]: # normalize median and q3
  #   if i==1:
  #     common_var = common_variances_med
  #   else:
  #     common_var = common_variances_q3
  #   threshold_med = create_table(prot, LFQ = False, normalize = i)
  #   fold_change(prot, False, i)
  #   prot_absent_controls(prot, threshold_med, False, i)
  #   test_normal_equal_var(prot, threshold_med, False, i, common_var)
  #   plot_log10_abundance(prot, threshold_med, LFQ = False, normalize = i, common_variance = common_var )
    # create_entire_final_csv(prot, False, i)
    # create_significant_prot_final_csv(prot, False, i)
# for prot_name in PROTEINS:
#   create_putative_proteins_analysis(prot_3_rep, prot_name, False, 1)
#   create_putative_proteins_analysis(prot_3_rep, prot_name, False, 4)
  
  # test_normal_equal_var(prot, threshold_q1, False, 3)
  # plot_log10_abundance(prot, threshold_q1, LFQ = False, normalize = 3)
  # test_normal_equal_var(prot, threshold_q3, False, 4)
  # plot_log10_abundance(prot, threshold_q3, LFQ = False, normalize = 4)
#   # non_eq_var_ttest(prot, threshold_bait, LFQ = False, normalize = 2)
#   # fold_change(prot, False, 2)
#   # plot_log10_abundance(prot, threshold_bait, LFQ = False, normalize = 2)
#   # non_eq_var_ttest(prot, threshold_raw, LFQ = False, normalize = 0)
#   # fold_change(prot, False, 0)
#   # non_eq_var_ttest(prot, threshold_med, LFQ = False, normalize = 1)
#   test_normal_equal_var(prot, threshold_raw, False, 0)
#   # fold_change(prot, False, 1)
#   plot_log10_abundance(prot, threshold_raw, LFQ = False, normalize = 0)

# plot_var_repartition(prot_3_rep, threshold, False, 1)
# plot_rep_proteins_per_batch(prot_3_rep, False, 4)
    # test_bart_per_all_prots([prot], threshold_med, False, 1) #med
# test_bart_per_all_prots(prot_3_rep[0:1], threshold_q3, False, 4) #q3
# test_bart_per_all_prots(prot_3_rep[0:1], threshold_raw, False, 0) #raw
