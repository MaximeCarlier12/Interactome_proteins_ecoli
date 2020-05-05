import data as dt
from Bio import Entrez
import csv
import numpy as np
import pandas as pd
import os
import glob

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

test_bioP()

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
  header('Genes for NCBIprot')
  print(gene)
  df['Accession_number'] = access # recupere tous les numéros d'accession
  df['Protein_name'] = name # recupere tous les numéros d'accession
  df['Organism'] = organism # recupere tous les numéros d'accession
  df['Gene_name'] = gene # recupere tous les numéros d'accession
  df = df.reset_index()
  df = df.set_index('Accession_number')
  return df

def create_csv(bfNumber):
# W1, W2 sprot et G2 ncbiprot
# A1, E1 ncbinr I4 ncbiprot
# F2, U6, U7 : ncbiprot
  df = dt.load_df(bfNumber)
  df = add_columns_df(df, accession_list(df))
  os.chdir("/home/carlier/Documents/Stage/Interactome_proteins_ecoli/") # where is the MS data
  path_batch = "MS data csv reportbuilder-20200408T212019Z-001/MS data csv reportbuilder/"
  if len(bfNumber) == 1:
      letter = ''
      number = bfNumber
  else :
    letter = bfNumber[0]
    number = bfNumber[1]
  full_path = glob.glob(path_batch+"batch "+letter+"/*"+letter+number+".reportbuilder.csv")
  df.to_csv(path_batch+"batch "+letter+"/"+letter+number+'.csv')

#create_csv('F2')
#create_csv('U6')
#create_csv('U7')

