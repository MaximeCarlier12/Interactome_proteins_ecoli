# conda activate myenv
# conda env export > environment.yml

import numpy as np
import pandas as pd
import os
import glob

os.chdir("/home/carlier/Documents/Stage/Interactome_proteins_ecoli/") # where is the MS data csv reportbuilder-20200408T212019Z-001 folder

CONDITION = ['LB log', 'LB O/N' ,'M9 0.2% ac O/N']
CONTROL = ['MG1655 (TYPE A)', 'MG1655 (placI)mVenus-SPA-pUC19 (TYPE C2)', 'MG1655']
PROTEIN = ['DnaA', 'DiaA', 'Hda', 'SeqA', 'HolD', 'DnaB', 'DnaG', 'NrdB']

def header(msg):
  print('-'*50)
  print('---- '+msg+' :')

def load_df(bfNumber):
  '''From a batch file number, it gets the right corresponding dataframe.'''

  path_batch = "MS data csv reportbuilder-20200408T212019Z-001/MS data csv reportbuilder/"
  # full_path get the only good file and can be in two different ways (containing 'A_5' or 'A5').
  if len(bfNumber) == 1:
    letter = ''
    number = bfNumber
    full_path = glob.glob(path_batch+"batch 1/*"+"_"+number+".reportbuilder.csv")
  else :
    letter = bfNumber[0]
    number = bfNumber[1]
    full_path = glob.glob(path_batch+"batch "+letter+"/*"+letter+"_"+number+".reportbuilder.csv") # A_5
    if full_path == []:
      full_path = glob.glob(path_batch+"batch "+letter+"/*"+letter+number+".reportbuilder.csv")
  df = pd.read_csv(full_path[0], sep=',', header=23) # A5
# We remove the first lines because they are not part of data. 
  if not "Accession" in df.columns: # some files have one additional row in the description part. 
    df = pd.read_csv(full_path[0], sep=',', header=24)
  df = df.set_index("Accession")
  df = df[df.Description.str.contains("Escherichia coli", case=False)] # Not case sensitive
  #print('empai', sum(df['emPAI']), df['Database'].iloc[1])
  return df

def load_sample_code():
  '''get sample code part of MS xlsx file'''

  MS = pd.read_excel("MS PPI screen - sample coding, 200226.xlsx", sheet_name = 0, header = 0)
  print(MS.iloc[0:10]['sample code'])
  return MS

def load_based_screen_controls(): 
  '''get controls part of based screen in MS xlsx file.
Returns a panda containing 3 rows : the resin unspecific control, the tag unspecific control and the untagged whole cell lysates (which is not a control)'''

  df = pd.read_excel("MS PPI screen - sample coding, 200226.xlsx", sheet_name = 1, header = 2, skipfooter = 44-13)
  df = df.rename(columns = {'LB log.1':'repeat LB log','LB O/N.1' : 'repeat LB O/N', 'M9 0.2% ac O/N.1' : 'repeat M9 0.2% ac O/N', 'strain - tagged protein (all in MG1655 genetic background)':'protein'}) # rename repeat number for each condition. 
  df = df.drop(df.columns[[9,10,11,12]], axis=1)
  typeResin = df.loc[df['strain'] == 'MG1655 (TYPE A)']
  typeTag = df.loc[df['strain'] == 'MG1655 (placI)mVenus-SPA-pUC19 (TYPE C2)']
  typeAbundant = df.loc[df['strain'] == 'MG1655']
  return pd.concat([typeResin, typeTag, typeAbundant])

def load_based_screen_samples(): 
  '''get samples part of based screen in MS xlsx file'''

  df = pd.read_excel("MS PPI screen - sample coding, 200226.xlsx", sheet_name = 1, header = 16)
  df = df.rename(columns = {'LB log.1':'repeat LB log','LB O/N.1' : 'repeat LB O/N', 'M9 0.2% ac O/N.1' : 'repeat M9 0.2% ac O/N', 'strain - tagged protein (all in MG1655 genetic background)':'protein'}) # rename repeat number for each condition. 
  return df

def get_control_replicates(pd_control, type_control, condition):
  ''' From the pd.control file, get pd files of the different replicates available for a certain type of control and a certain type of condition.'''
  if not condition in CONDITION : 
    print("Error condition")
    return None
  if not type_control in CONTROL : 
    print("Error control")
    return None
  control_row = pd_control[pd_control.strain== type_control] # get only one row.
  #repeat_condition = 'repeat '+condition
  batches = control_row[condition].iloc[0]  # get batches of the different replicates.
  batches = batches.replace(';',',') # replace ';' by ',' because there are two separator ways in the initial file. 
  batches = batches.replace(" ", "") # remove spaces
  batches = batches.split(',') # separate the replicates
#  header('batches')    
#  print(batches)
  return [load_df(batches[i]) for i in range(len(batches))]

def get_batches(pd_samples, protein, condition):
  ''' MS is the sample file. From one protein and one condition, it gets all batch names.'''
  if not condition in CONDITION : 
    print("Error condition")
    return None
  if not protein in PROTEIN : 
    print("Error protein")
    return None
  prot_row = pd_samples[pd_samples.protein== protein] # get only one row.
  repeat_condition = 'repeat '+condition
  if prot_row[repeat_condition].iloc[0] != 3:
    print("Be carreful, there are not 3 replicates")
  batches = prot_row[condition].iloc[0]  # get batches of the different replicates.
  batches = batches.replace(';',',') # replace ';' by ',' because there are two separator ways in the initial file. 
  batches = batches.replace(" ", "") # remove spaces
  batches = batches.split(',') # separate the replicates
#  header('batches')    
#  print(batches)
  return batches

def get_protein_replicates(pd_samples, protein, condition):
  batches = get_batches(pd_samples, protein, condition)
  return [load_df(batches[i]) for i in range(len(batches))]

def get_replicates(pd, obj, condition, control=False):
  ''' Global function with a boolean indicating if we want replicates of a control or a protein.
Obj : control or protein.'''

  if control == True:
    return get_control_replicates(pd, obj,condition)
  else:
    return get_protein_replicates(pd, obj, condition)

def intersect_index(replicates):
  ''' From a list of replicates, gives a list of accession number proteins present in all replicates.'''
  intersect = replicates[0].index
  for i in replicates : 
    intersect = intersect.intersection(i.index)
  return intersect

#df = load_df("A1")
#load_sample_code()
all_controls = load_based_screen_controls()
one_control = get_control_replicates(all_controls, 'MG1655','LB log')
all_samples = load_based_screen_samples()
#[rep1,rep2,rep3] = get_protein_replicates(all_samples, 'DnaB', 'M9 0.2% ac O/N')
#intersect = intersect_index([rep2,rep3])
#header('intercept')
#print("rep3",len(rep3.index), "rep2", len(rep2.index), "intersect", len(intersect))

def count_databases():
  '''Count intersection for each file'''

  header("test database")
  nb_inter = []
  bd = []
  bd2 = []
  for p in PROTEIN :
    check = False;check_c = False
    for c in CONDITION : 
      check = False
      check_c = False
      rep = get_protein_replicates(all_samples, p, c)
      inter = len(intersect_index(rep))
      if inter != 0:
        check = True
      nb_inter.append(inter)
      if check == False :
        print("Error for", p, "in condition", c) # for each protein, in there at least one control : no
      if check == True:
        bd.append(rep[0]['Database'][0])
        ctr = 'MG1655 (placI)mVenus-SPA-pUC19 (TYPE C2)'
        rep = get_control_replicates(all_controls, ctr, c)
        bd2.append(rep[0]['Database'][0])
        inter = len(intersect_index(rep))
        if inter != 0:
          check_c = True
        else : print("error for control in condition", c)

  print(nb_inter, len(nb_inter))
  nb_inter1 = [nb_inter[i] for i in range(len(nb_inter)) if nb_inter[i] == 0]
  nb_inter2 = []
  print(bd)
  print(bd2)
  print([bd[i] == bd2[i] for i in range(len(bd))])
  for c in CONDITION : 
    check = False
    for ctr in CONTROL :
      rep = get_control_replicates(all_controls, ctr, c)
      inter = len(intersect_index(rep))
      if inter != 0:
        check = True
      nb_inter2.append(inter)
    print(check) # for each condition, in there at least one control
  print(nb_inter2, len(nb_inter2))
  nb_inter3 = [nb_inter2[i] for i in range(len(nb_inter2)) if nb_inter2[i] == 0]
  print("number of protein with different databases : ", len(nb_inter1), "out of", len(nb_inter))
  print("number of controls with different databases :", len(nb_inter3), "out of", len(nb_inter2))

def good_proteins():
  '''List of tuple containing the protein to be used and its condition'''
  bdp = ''; bdc = ''
  good_proteins = []
  for p in PROTEIN :
#    check = False; check_c = False
    for c in CONDITION : 
      check = False
      check_c = False
      rep = get_protein_replicates(all_samples, p, c)
      inter = len(intersect_index(rep))
      if inter != 0:
        check = True
        bdp = rep[0]['Database'][0]
        ctr = 'MG1655 (placI)mVenus-SPA-pUC19 (TYPE C2)'
        rep = get_control_replicates(all_controls, ctr, c)
        bdc = rep[0]['Database'][0]
        if bdp == bdc:
          good_proteins.append((p,c))
  return good_proteins
#print(good_proteins())
#count_databases()
