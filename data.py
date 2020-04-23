# conda activate myenv
# conda env export > environment.yml

import numpy as np
import matplotlib as mt
import pandas as pd
import os
import glob
os.chdir("/home/carlier/Documents/Stage/Interactome_proteins_ecoli/")

CONDITION = ['LB log', 'LB O/N' ,'M9 0.2% ac O/N']
PROTEIN = ['DnaA', 'DiaA', 'Hda', 'SeqA', 'HolD', 'DnaB', 'DnaG', 'NrdB']

def header(msg):
  print('-'*50)
  print('[ '+msg+' ]')

def load_df(bfNumber):
  '''From a batch file number, it gets the right corresponding dataframe.'''
  path_batch = "MS data csv reportbuilder-20200408T212019Z-001/MS data csv reportbuilder/"
  # full_path get the only good file and can be in two different ways (containing 'A_5' or 'A5').
  full_path = glob.glob(path_batch+"batch "+bfNumber[0]+"/*"+bfNumber[0]+"_"+bfNumber[1]+".reportbuilder.csv")
  if full_path == []:
    full_path = glob.glob(path_batch+"batch "+bfNumber[0]+"/*"+bfNumber[0]+bfNumber[1]+".reportbuilder.csv")
  df = pd.read_csv(full_path[0], sep=',', header=23)
# We remove the first lines because they are not part of data. 
  df = df.set_index("Accession")
  #header("df.head")
  #print(df.head) #check it is correct. 
  return df

def load_sample_code():
  '''get sample code part of MS xlsx file'''

  MS = pd.read_excel("MS PPI screen - sample coding, 200226.xlsx", sheet_name = 0, header = 0)
  print(MS.iloc[0:10]['sample code'])
  return MS

def load_based_screen_controls(): 
  '''get controls part of based screen in MS xlsx file.
Returns a list of pandas containing : the resin unspecific control, the tag unspecific control and the untagged whole cell lysates'''

  df = pd.read_excel("MS PPI screen - sample coding, 200226.xlsx", sheet_name = 1, header = 2, skipfooter = 44-13)
  df = df.rename(columns = {'LB log.1':'repeat LB log','LB O/N.1' : 'repeat LB O/N', 'M9 0.2% ac O/N.1' : 'repeat M9 0.2% ac O/N', 'strain - tagged protein (all in MG1655 genetic background)':'protein'}) # rename repeat number for each condition. 
  df = df.drop(df.columns[[9,10,11,12]], axis=1)
  typeResin = df.loc[df['strain'] == 'MG1655 (TYPE A)']
  typeTag = df.loc[df['strain'] == 'MG1655 (placI)mVenus-SPA-pUC19 (TYPE C2)']
  typeAbundant = df.loc[df['strain'] == 'MG1655']
  return [typeResin, typeTag, typeAbundant]

def load_based_screen_samples(): 
  '''get samples part of based screen in MS xlsx file'''

  df = pd.read_excel("MS PPI screen - sample coding, 200226.xlsx", sheet_name = 1, header = 16)
  header("index")
  df = df.rename(columns = {'LB log.1':'repeat LB log','LB O/N.1' : 'repeat LB O/N', 'M9 0.2% ac O/N.1' : 'repeat M9 0.2% ac O/N', 'strain - tagged protein (all in MG1655 genetic background)':'protein'}) # rename repeat number for each condition. 
  return df

def get_replicates(MS, protein, condition):
  ''' From one protein and one condition, it gets all replicates dataframes.'''

  prot_row = MS[MS.protein=='DnaA']
  print(prot_row.dtypes)
  repeat_condition = 'repeat '+condition
  if prot_row[repeat_condition].iloc[0] == 3:
    batches = prot_row[condition].iloc[0]  # get batches of the different replicates.
    batches = batches.replace(';',',') # replace ';' by ','
    batches = batches.replace(" ", "") # remove spaces
    batches = batches.split(',') # separate the replicates
    header('batches')    
    print(batches)
    return [load_df(batches[0]), load_df(batches[1]), load_df(batches[2])]
  else :
    print("Not enough replicates")
    return None

#df = load_df("A1")
#load_sample_code()
controls = load_based_screen_controls()
samples = load_based_screen_samples()
#replicates = get_replicates(MS, 'DnaA', 'M9 0.2% ac O/N')
#print(replicates)


