import Load_PPI_screen as dt
import Data_processing as dp
import handle as hd
import visualization as vz
import qualitatif_stats as st
import csv
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import glob

#pd_samples = dt.load_based_screen_samples()
#pd_controls = dt.load_based_screen_controls()
used_prot_tuples = dp.good_proteins()
contaminant_genes = dt.load_contaminant_list()
#controls_typeC = {'LB log':['S1','S2','S3'], 'LB O/N':['O9', 'R4', 'R5'], 'M9 0.2% ac O/N':['R1', 'R2', 'R3']}
#controls_typeA = {'LB log':['L1', 'T7', 'T8', 'T9'], 'LB O/N':['C13', 'P5', 'U10'], 'M9 0.2% ac O/N':['A10', 'T5', 'T6']}

hd.create_table(used_prot_tuples)
threshold = hd.create_condensed_table(used_prot_tuples, 0.25)
for prot1 in used_prot_tuples:
  dt.header(prot1[0]+' in '+prot1[1])
  st.log_calculation(prot1)
  st.test_normal_equal_var(prot1, threshold)
  vz.plot_log2_emPAI(prot1, threshold, contaminant_genes)
#for prot1 in used_prot_tuples:
#  df = hd.load_df_table(prot1, True)
#  print(prot1[0], prot1[1][:6]+'\t', 'TestC', df[df.C_is ==True].shape[0],'\t', 'TestA', df[df.A_is ==True].shape[0])

plt.close('all')
