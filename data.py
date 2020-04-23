# conda activate myenv

import numpy as np
import matplotlib as mt
import pandas as pd
import os

os.chdir("/home/carlier/Documents/Stage/Interactome_proteins_ecoli/")
def header(msg):
  print('-'*50)
  print('[ '+msg+' ]')

def load_df(batchNumber, fileNumber):
  path_batch = "MS data csv reportbuilder-20200408T212019Z-001/MS data csv reportbuilder/"
  header_names = ['Family,Member','Database','Accession','Score','Mass','Num. of matches','Num. of significant matches','Num. of sequences','Num. of significant sequences','emPAI','Description']
  df = pd.read_csv(path_batch+"batch "+batchNumber+"/609201756glin_"+batchNumber+"_"+fileNumber+".reportbuilder.csv", sep=',', names=header_names, skiprows=[i for i in range(0,27)])
# We remove the first lines because they are not part of data. 
  df = df.set_index("Accession")
  header("df.head")
  print(df.dtypes) # check it is correct. 

load_df("A","1")
