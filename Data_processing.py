import Load_PPI_screen as dt
from Bio import Entrez
import csv
import numpy as np
import pandas as pd
import os
import glob
os.chdir("/home/carlier/Documents/Stage/Interactome_proteins_ecoli/") # where is the MS data

all_controls = dt.load_based_screen_controls()
all_samples = dt.load_based_screen_samples()

def header(msg):
  print('-'*50)
  print('---- '+msg+' :')

def accession_list(df):
  prot = []
  string = ''
  for i in df.index.values:
    string += ','
    string += i
  string = string.replace(',','',1) # remove the first ","
  return string

def test_bioP():
  bfNumber = 'W1'
  df = dt.load_df(bfNumber)
  string = accession_list(df)
  header('bioP')
  Entrez.email = 'maxime.Carlier@insa-lyon.fr'
#  handle = Entrez.efetch(db="protein", id="gi|490117830", rettype="gp", retmode="xml")
  handle = Entrez.efetch(db='protein', id='gi|1525625145', rettype="gp", retmode="xml") # Nprot
  #handle = Entrez.efetch(db="protein", id="P0ABD5", rettype="gb", retmode="xml") #Sprot
  record = Entrez.read(handle)
  header('different keys')
  print(record[0].keys()) # [0] to have the first id given. 
  header('source_db')
  print(record[0]['GBSeq_source'])
  header('accession number')
  print(record[0]['GBSeq_other-seqids'])
#  print(record[0]['GBSeq_organism'])
#  print(record[0]['GBSeq_length'])
  header('features')
  print(record[0]['GBSeq_feature-table'])
  header('feature Protein')
  print(record[0]['GBSeq_feature-table'][2])
#  a = []
#  for i in range(len(record)):
#    b =[]
#    for feature in record[i]['GBSeq_feature-table']:
#      if feature['GBFeature_key'] in ['CDS', 'Protein', 'gene']:
#        b.append(feature['GBFeature_key'])
#    print((b))    
#    a.append(b)
  handle.close()

#test_bioP()

def add_columns_df(df, string):
  access = []
  definition = []
  name = []
  organism = []
  gene = []
  Entrez.email = 'maxime.Carlier@insa-lyon.fr'
  handle = Entrez.efetch(db="protein", id=string, rettype="gb", retmode="xml") # Nprot
  result = Entrez.read(handle)
  for res in result :
    accession = ""
    for i in res['GBSeq_other-seqids']:
      if 'gi' in i:
        accession = i
    access.append(accession)
    definition.append(res['GBSeq_definition'])
    is_gene = False
    is_prot = False
    for feature in res['GBSeq_feature-table']:
      if feature['GBFeature_key'] == 'Protein':
        name.append(feature['GBFeature_quals'][0]['GBQualifier_value'])
        is_prot = True
      elif feature['GBFeature_key'] == 'gene':
        if is_gene == False :
          gene.append(feature['GBFeature_quals'][0]['GBQualifier_value'])        
        is_gene = True
      elif feature['GBFeature_key'] == 'source':
        organism.append(feature['GBFeature_quals'][0]['GBQualifier_value'])
      elif feature['GBFeature_key'] == 'CDS':
        for j in feature['GBFeature_quals']:
          if j['GBQualifier_name']== 'gene':    
            gene.append(j['GBQualifier_value'])
            is_gene = True
    if is_gene == False : 
      gene.append('None')
    if is_prot == False : 
      name.append('None')
  handle.close()
  df['Accession_number'] = access # recupere tous les numéros d'accession
  df['Protein_name'] = name # recupere tous les numéros d'accession
  df['Organism'] = organism # recupere tous les numéros d'accession
  df['Gene_name'] = gene # recupere tous les numéros d'accession
  df = df.reset_index()
  df = df.set_index('Accession_number')
  return df

def create_csv_genes(bfNumber):
# W1, W2 sprot et G2 ncbiprot
# A1, E1 ncbinr I4 ncbiprot
# F2, U6, U7 : ncbiprot
  df = dt.load_df(bfNumber)
  df = add_columns_df(df, accession_list(df))
  df = df[~df.Gene_name.str.contains("None", case=False)] # Remove rows without genes.
  df = df[df.Organism.str.contains("Escherichia coli", case=False)] # Remove rows wrong organism
  for i in pd.unique(df['Family']): # Remove redundant rows, keep max sig seq. 
    if len(df.loc[df['Family'] == i]) > 1 :
      df_fam = df.loc[df['Family'] == i]
      df_fam = df_fam.sort_values(by=['Num. of significant sequences','Num. of significant matches'] , ascending=False)
      indexNames = df_fam.iloc[1:].index
      df.drop(indexNames , inplace=True)

  path_batch = "MS data csv reportbuilder-20200408T212019Z-001/MS data csv reportbuilder/"
  if len(bfNumber) == 1:
      letter = ''
      number = bfNumber
  else :
    letter = bfNumber[0]
    number = bfNumber[1:]
  full_path = glob.glob(path_batch+"batch "+letter+"/*"+letter+number+".reportbuilder.csv")
  df.to_csv(path_batch+"batch "+letter+"/genes_"+letter+number+'.csv')

#create_csv_genes('T10')

def create_csv_genes_good_proteins():
  header('proteins')
  all_samples = dt.load_based_screen_samples()
  for my_tuple in dt.good_proteins():
    batch_names = dt.get_batches(all_samples,my_tuple[0], my_tuple[1])
    for bname in batch_names:
      print(bname)
      create_csv_genes(bname)
  header('controls')
  controls_typeC = ['S1','S2','S3', 'O9', 'R4', 'R5', 'R1', 'R2', 'R3']
  for bname in controls_typeC:
    print(bname)
    create_csv_genes(bname)
  controls_typeA = ['L1', 'T7', 'T8', 'T9', 'C13', 'P5', 'U10', 'A10', 'T5', 'T6']
  for bname in controls_typeA:
    create_csv_genes(bname)

#create_csv_genes_good_proteins()

def load_df_genes(bfNumber):
  '''Load dataframe from new files with all gene names.'''
  os.chdir("/home/carlier/Documents/Stage/Interactome_proteins_ecoli/") # where is the MS data
  path_batch = "MS data csv reportbuilder-20200408T212019Z-001/MS data csv reportbuilder/"
  letter = bfNumber[0]
  number = bfNumber[1:]
  full_path = path_batch+"batch "+letter+"/genes_"+letter+number+'.csv'
  df = pd.read_csv(full_path, sep=',', header=0)
  df = df.set_index("Accession_number")
  return df

def remove_duplicate_genes(rep):
  '''One replicate '''
  unique_genes = rep.Gene_name.unique()
  a = rep[rep.Gene_name.isin(unique_genes)]
  if len(a.index) != len(unique_genes):
    duplicate = a.groupby('Gene_name').count()
    indexNames = duplicate[duplicate.Accession == 1].index
    duplicate = duplicate.drop(indexNames)
    for i in duplicate.index:
      df_duplicate = a.loc[a['Gene_name'] == i]
      indexNames = df_duplicate.iloc[1:].index
      a = a.drop(indexNames)
  return a

#rep = load_df_genes('U3')
#print(len(rep))
#res = remove_duplicate_genes(rep)
#print(len(res))

def create_csv_unique_gene(bfNumber):
  df = load_df_genes(bfNumber)
  df = remove_duplicate_genes(df)
  if len(bfNumber) == 1:
      letter = ''
      number = bfNumber
  else :
    letter = bfNumber[0]
    number = bfNumber[1:]
  path_batch = "MS data csv reportbuilder-20200408T212019Z-001/MS data csv reportbuilder/"
  df.to_csv(path_batch+"batch "+letter+"/unique_gene_"+letter+number+'.csv')

def create_all_csv_unique_gene():
  used_prot_tuples = dt.good_proteins()
  for prot1 in used_prot_tuples:
    prot_batches = dt.get_batches(all_samples, prot1[0], prot1[1])
    for bname in prot_batches:
      print(bname)
      create_csv_unique_gene(bname)
  controls_typeC = ['S1','S2','S3', 'O9', 'R4', 'R5', 'R1', 'R2', 'R3']
  for bname in controls_typeC:
    create_csv_unique_gene(bname)
  controls_typeA = ['L1', 'T7', 'T8', 'T9', 'C13', 'P5', 'U10', 'A10', 'T5', 'T6']
  for bname in controls_typeA:
    create_csv_unique_gene(bname)
  

#create_all_csv_unique_gene()
# be careful, C13 and A10 are in NCBInr database (controls_typeA in LB O/N and M9 0.2% ac O/N)

