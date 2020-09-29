import Load_PPI_screen as dt
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
from scipy import stats

controls_typeC = {'LB log':['S1','S2','S3'], 'LB O/N':['O9', 'R4', 'R5'], 'M9 0.2% ac O/N':['R1', 'R2', 'R3']}
controls_typeA = {'LB log':['T7', 'T8', 'T9'], 'LB O/N':['P5', 'U10', 'X17'], 'M9 0.2% ac O/N':['T5', 'T6', 'X16']}
missing_files = ['A1', 'A2', 'M4', 'M5', 'A9', 'A10', 'I1', 'I2']
# missing_files = []

prey = {'DnaA':'diaA', 'DiaA':'dnaA', 'DnaB':'dnaC', 'DnaG':'dnaB', 'NrdB':'nrdA', 'HolD':['dnaE', 'dnaN', 'dnaQ', 'dnaX', 'holA', 'holB', 'holC', 'holE'], 'SeqA':'', 'Hda':''} # no preys for hda and seqA.

interesting_prey = {'DnaA':['purR', 'eno'], 'DiaA':['gldA', 'ndh', 'wbbK', 'rfaF', 'rfaB', 'rfaG','rfaP','RfaD','rfaB','gmhA'], 'DnaB':'tdcB', 'DnaG':['nrdB', 'glgB', 'amyA', 'glgA', 'seqA'], 'NrdB':['dnaN', 'skp'], 'HolD':['topB'], 'SeqA':'rfaD', 'Hda': ''} # no prey for Hda.

CONDITION = ['LB log', 'LB O/N' ,'M9 0.2% ac O/N']
CONTROL = ['MG1655 (TYPE A)', 'MG1655 (placI)mVenus-SPA-pUC19 (TYPE C2)', 'MG1655']
PROTEINS = ['DnaA', 'DiaA', 'Hda', 'SeqA', 'HolD','DnaB', 'DnaG', 'NrdB']

mydir = "/home/carlier/Documents/Stage/Interactome_proteins_ecoli/"
os.chdir(mydir)
# folder_to_save_all_results = 'New_data/'
folder_to_save_all_results = 'New_data_4_replicates/'

pd_samples = dt.load_based_screen_samples()
pd_controls = dt.load_based_screen_controls()
contaminant_genes = dt.load_contaminant_list()

# Samples new codes can be for new data (first september or with 4 replicates (25 september)).
# pd_samples_new_codes = dt.load_new_batch_codes_samples()
data_codes = {'protein': ['DiaA', 'SeqA', 'HolD', 'DnaG', 'Hda'], 'LB log': ['A1, E1, I4, X6', 'F4, U4, U5, X13, X14', '', '', ''], 'LB O/N': ['A2, I5, P1, X7', '', '2, F8, U11, X4', '', ''], 'M9 0.2% ac O/N': ['A3, I6, T3', '','', 'B11, F7, H3, X15', 'B12, I3, T4, X8']}
pd_samples_new_codes = pd.DataFrame(data_codes, columns=['protein', 'LB log', 'LB O/N', 'M9 0.2% ac O/N'])

print(pd_samples_new_codes)
pd_controls_new_codes = dt.load_new_batch_codes_controls()

normalization_types = {0:'Raw intensity', 1:'Median normalization', 2:'Bait protein normalization'}
normalize_values = {0:'Not normalized', 1:'Median normalization', 2:'Bait normalization',  3:'Q1 normalization',  4:'Q3 normalization'}

params = {
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
   'axes.labelsize': 12.5, # names of axes
   'font.size': 11,
   'legend.fontsize': 9, # size legend
   'xtick.labelsize': 10, # values in axes
   'ytick.labelsize': 10,
   'text.usetex': False,
   'axes.linewidth':1.5, #0.8
   'axes.titlesize':18, # title of graph
   'axes.spines.top':True,
   'axes.spines.right':True,
   'font.family': "Arial"
   }
plt.rcParams.update(params)