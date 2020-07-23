from globvar import *
import Load_PPI_screen as dt
import Data_processing as dp
#from matplotlib.ticker import FormatStrFormatter
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
used_prot_tuples = dp.good_proteins() 
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
  path = "MS data csv reportbuilder-20200408T212019Z-001/contaminant_or_not.csv" 
#  df = pd.read_csv(path, header=0, sep =';', usecols = lambda column : column not in ['Potential contaminant','Not contaminant'])
  df = pd.read_csv(path, header=0, sep =';')
  df = df.set_index(['Bait protein', 'Condition'])
  return df

def load_df_table_maxQ(prot1, LFQ, normalize):
  '''Load dataframe from new files with all gene names for raw or LFQ.'''
  path_batch = "maxQ/SEGREGATED-20200619T092017Z-001/Protein_table/"
  if LFQ == True:
    full_path = path_batch+ prot1[0]+"_"+prot1[1].replace('/', '_')+'_LFQ.csv'
  elif normalize == 0:
    full_path = path_batch+ prot1[0]+"_"+prot1[1].replace('/', '_')+'_Int.csv'
  elif normalize == 1:
    full_path = path_batch+ prot1[0]+"_"+prot1[1].replace('/', '_')+'_Int_Norm_Med.csv'
  elif normalize == 2:
    full_path = path_batch+ prot1[0]+"_"+prot1[1].replace('/', '_')+'_Int_Norm_Bait.csv'
  elif normalize == 3:
    full_path = path_batch+ prot1[0]+"_"+prot1[1].replace('/', '_')+'_Int_Norm_Q1.csv'
  elif normalize == 4:
    full_path = path_batch+ prot1[0]+"_"+prot1[1].replace('/', '_')+'_Int_Norm_Q3.csv'
  else:
    print('error in normalize value')
  df = pd.read_csv(full_path, sep=',', header=0, index_col = 0)
  return df
