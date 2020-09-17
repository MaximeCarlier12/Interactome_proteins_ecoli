import Load_PPI_screen as dt
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
from scipy import stats

controls_typeC = {'LB log':['S1','S2','S3'], 'LB O/N':['O9', 'R4', 'R5'], 'M9 0.2% ac O/N':['R1', 'R2', 'R3']}
controls_typeA = {'LB log':['T7', 'T8', 'T9'], 'LB O/N':['P5', 'U10', 'X17'], 'M9 0.2% ac O/N':['T5', 'T6', 'X16']}
missing_files = ['A1', 'A2', 'M4', 'M5', 'A9', 'A10', 'I1', 'I2']

prey = {'DnaA':'diaA', 'DiaA':'dnaA', 'DnaB':'dnaC', 'DnaG':'dnaB', 'NrdB':'nrdA', 'HolD':['dnaE', 'dnaN', 'dnaQ', 'dnaX', 'holA', 'holB', 'holC', 'holE'], 'SeqA':'', 'Hda':''} # no preys for hda and seqA.

interesting_prey = {'DnaA':['purR', 'eno'], 'DiaA':['gldA', 'ndh', 'wbbK', 'rfaF', 'rfaB', 'rfaG','rfaP','RfaD','rfaB','gmhA'], 'DnaB':'tdcB', 'DnaG':['nrdB', 'glgB', 'amyA', 'glgA', 'seqA'], 'NrdB':['dnaN', 'skp'], 'HolD':['topB'], 'SeqA':'rfaD', 'Hda': ''} # no prey for Hda.

CONDITION = ['LB log', 'LB O/N' ,'M9 0.2% ac O/N']
CONTROL = ['MG1655 (TYPE A)', 'MG1655 (placI)mVenus-SPA-pUC19 (TYPE C2)', 'MG1655']
PROTEINS = ['DnaA', 'DiaA', 'Hda', 'SeqA', 'HolD','DnaB', 'DnaG', 'NrdB']

mydir = "/home/carlier/Documents/Stage/Interactome_proteins_ecoli/"
os.chdir(mydir)

pd_samples = dt.load_based_screen_samples()
pd_controls = dt.load_based_screen_controls()
contaminant_genes = dt.load_contaminant_list()

pd_samples_new_codes = dt.load_new_batch_codes_samples()
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