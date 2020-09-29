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
  '''From a batch file number, it gets the right corresponding dataframe. emPAI data'''

  path_batch = "MS data csv reportbuilder-20200408T212019Z-001/MS data csv reportbuilder/"
  # full_path get the only good file and can be in two different ways (containing 'A_5' or 'A5').
  if len(bfNumber) == 1:
    letter = ''
    number = bfNumber
    full_path = glob.glob(path_batch+"batch 1/*"+"_"+number+".reportbuilder.csv")
  else :
    letter = bfNumber[0]
    number = bfNumber[1:]
    full_path = glob.glob(path_batch+"batch "+letter+"/*"+letter+"_"+number+".reportbuilder.csv") # A_5
    if full_path == []:
      full_path = glob.glob(path_batch+"batch "+letter+"/*"+letter+number+".reportbuilder.csv")
  df = pd.read_csv(full_path[0], sep=',', header=23) # A5
# We remove the first lines because they are not part of data. 
  if not "Accession" in df.columns: # some files have one additional row in the description part. 
    df = pd.read_csv(full_path[0], sep=',', header=24)
  df = df.set_index("Accession")
#  df = df[df.Description.str.contains("Escherichia coli", case=False)] # Not case sensitive
  #print('empai', sum(df['emPAI']), df['Database'].iloc[1])
  return df

def load_sample_code():
  '''get sample code part of MS xlsx file. For each bf number, the protein name and growth conditions are mentionned.'''

  MS = pd.read_excel("MS PPI screen - sample coding, 200226.xlsx", sheet_name = 0, header = 0)
  print(MS.iloc[0:10]['sample code'])
  return MS

def load_contaminant_list():
  '''Returns a list that contains all potential contaminant genes'''
  MS = pd.read_excel("common contaminants PPI data.xlsx", sheet_name = 0, header = 0)
#  print(MS.iloc[0:10]['gene name'])
#  print(MS.loc[1])
  gene_cont = list(MS['gene name'])
  gene_cont = [x for x in gene_cont if str(x) != 'nan']
  for i in gene_cont:
    if len(i.split(',')) > 1:
      gene_cont.remove(i)
      gene_cont += i.split(',')
  return gene_cont

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
  '''get samples part of based screen in MS xlsx file. For each protein, there are batch numbers for the different growth conditions.'''

  df = pd.read_excel("MS PPI screen - sample coding, 200226.xlsx", sheet_name = 1, header = 16)
  df = df.rename(columns = {'LB log.1':'repeat LB log','LB O/N.1' : 'repeat LB O/N', 'M9 0.2% ac O/N.1' : 'repeat M9 0.2% ac O/N', 'strain - tagged protein (all in MG1655 genetic background)':'protein'}) # rename repeat number for each condition. 
  return df

def get_batches(pd_samples, protein, condition):
  ''' MS is the sample file. From one protein and one condition, it gets all batch names.'''
  if not condition in CONDITION : 
    print("Error condition")
    return None
  if not protein in PROTEIN :
    print("Error protein")
    return None
  prot_row = pd_samples[pd_samples.protein== protein] # get only one row.
  if len(prot_row) == 0:
    batches = []
  else:
    batches = prot_row[condition].iloc[0]  # get batches of the different replicates.
    batches = batches.replace(';',',') # replace ';' by ',' because there are two separator ways in the initial file. 
    batches = batches.replace(" ", "") # remove spaces
    batches = batches.split(',') # separate the replicates
#  header('batches')    
#  print(batches)
  return batches

def load_new_batch_codes_controls():
  '''Same function as load_based_screen_controls but for new analysis files.'''
  df = pd.read_excel("maxQ/New_data/ANALYSIS CODES .xlsx", sheet_name = 0, header = 5, skipfooter = 27-15)
  typeResin = df.loc[df['strain'] == 'MG1655 (TYPE 1 CONTROL)']
  typeResin = typeResin.replace('MG1655 (TYPE 1 CONTROL)', 'MG1655 (TYPE A)')
  typeTag = df.loc[df['strain'] == 'MG1655 (placI)mVenus-SPA-pUC19 (TYPE 2 CONTROL)']
  typeTag = typeTag.replace('MG1655 (placI)mVenus-SPA-pUC19 (TYPE 2 CONTROL)', 'MG1655 (placI)mVenus-SPA-pUC19 (TYPE C2)')
  typeAbundant = df.loc[df['strain'] == 'MG1655 (TYPE 3 CONTROL)']
  typeAbundant = typeAbundant.replace('MG1655 (TYPE 3 CONTROL)', 'MG1655')
  return pd.concat([typeResin, typeTag, typeAbundant])

def load_new_batch_codes_samples():
  '''Same function as load_based_screen_samples but for new analysis files.'''
  df = pd.read_excel("maxQ/New_data/ANALYSIS CODES .xlsx", sheet_name = 0, header = 18)
  df = df.rename(columns = {'strain - tagged protein (all in MG1655 genetic background) ':'protein'})
  return df
  