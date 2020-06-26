import Load_PPI_screen as dt
import Data_processing as dp
import handle as hd
import csv
import numpy as np
import matplotlib.pyplot as plt
import qualitatif_stats as st
import pandas as pd
import matplotlib.patches as mpatches
import os
import glob

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

prey = {'DnaA':'diaA', 'DiaA':'dnaA', 'DnaB':'dnaC', 'DnaG':'dnaB', 'NrdB':'nrdA', 'HolD':['dnaE', 'dnaN', 'dnaQ', 'dnaX', 'holA', 'holB', 'holC', 'holE'], 'SeqA':''}
interesting_prey = {'DnaA':['purR', 'eno'], 'DiaA':['gldA', 'ndh', 'wbbK', 'rfaF', 'rfaB', 'rfaG','rfaP','RfaD','rfaB','gmhA'], 'DnaB':'tdcB', 'DnaG':['nrdB', 'glgB', 'amyA', 'glgA', 'seqA'], 'NrdB':['dnaN', 'skp'], 'HolD':['topB'], 'SeqA':'rfaD'}

def sum_empai(used_prot_tuples):
  '''Print sum(emPAI) per replicate per protein studied)'''
  prots = []
  name_prots = []
  for prot1 in used_prot_tuples:
    df = load_df_table(prot1, True)
    sums = df[['Rep1', 'Rep2', 'Rep3', 'CtrC1', 'CtrC2', 'CtrC3', 'CtrA1', 'CtrA2', 'CtrA3']].sum(axis=0)
    name = prot1[0]+'_'+prot1[1][:4]
    name_prots.append(name)
    prots.append(sums)
    print(sums[0])
  for i in range(len(prots)):
    plt.scatter([name_prots[i]]*3, prots[i][:3] , label="Protein test" if i == 0 else "", color='royalblue')
    plt.scatter([name_prots[i]]*3, prots[i][3:6] , label="CtrC" if i == 0 else "", color='chartreuse')
    plt.scatter([name_prots[i]]*3, prots[i][6:] , label="CtrA" if i == 0 else "", color='red')
  plt.legend()
  plt.title('Differences of emPAI sums between files')
  plt.xlabel('Protein studied')
  plt.ylabel('sum(emPAI) value')
  plt.xticks(rotation=90)
  plt.grid(axis = 'x') # vertical lines
  manager = plt.get_current_fig_manager()
  manager.window.showMaximized()
  plt.tight_layout()
  plt.show()

def plot_emPAI(prot1, control = 'AC'):
  '''Plot emPAI value for each gene for controls and test in log scale.'''
  df = load_df_table(prot1, True)
  empRep = []; empCtr = []
  df['max_empai'] = df[['Rep1', 'Rep2', 'Rep3']].max(axis=1)
  df = df.sort_values(by = 'max_empai', ascending = False)
  maxval = max(df.max(axis = 1)) # max value of the table
  minval = min(df.max(axis = 1)) # max value of the table
  #  empRep.append(np.mean([row.Rep1, row.Rep2, row.Rep3]))
  #  empCtr.append(np.mean([row.Ctr1, row.Ctr2, row.Ctr3]))
  plt.scatter(df.index, df.Rep1, label="Protein test", color='royalblue')
  plt.scatter(df.index, df.Rep2, color='royalblue')
  plt.scatter(df.index, df.Rep3, color='royalblue')
  plt.title(prot1[0]+' in '+prot1[1])
  if 'A' in control:
    plt.scatter(df.index, df.CtrA1, label = "CtrA (without SPA tag)", color='red')
    plt.scatter(df.index, df.CtrA2, color='red')
    plt.scatter(df.index, df.CtrA3, color='red')
    if 'CtrA4' in df :
      plt.scatter(df.index, df.CtrA4, color='red')
  if 'C' in control:
    plt.scatter(df.index, df.CtrC1, label = "CtrC (with SPA tag)", color='chartreuse')
    plt.scatter(df.index, df.CtrC2, color='chartreuse')
    plt.scatter(df.index, df.CtrC3, color='chartreuse')
  plt.xticks(rotation=90)
  plt.ylim(-0.02, 100)
  plt.xlabel('Gene name')
  plt.ylabel('emPAI value')
  plt.yscale('symlog')
  plt.grid(axis = 'x') # vertical lines
  plt.legend()
  manager = plt.get_current_fig_manager()
  manager.window.showMaximized()
  plt.tight_layout()
  plt.show()

def plot_log2_emPAI(prot1, threshold, contaminant_genes, control = 'AC'):
  '''Plot log2(emPAI) value for each gene for controls and test.'''

  fig,ax = plt.subplots()
  width = 12.5 ; height = 6.35 # taille finale de ta figure png
#  width = 14.4 ; height = 7.15 # taille finale de ta figure svg
  fig.set_size_inches(width, height)
  var_test = hd.load_df_equal_test()
  df = hd.load_df_table(prot1, True)
  df['max_empai'] = df[['Rep1', 'Rep2', 'Rep3']].max(axis=1)
  minval = df[['Rep1', 'Rep2', 'Rep3']].min().min() # min value of the table
  maxval = df[['Rep1', 'Rep2', 'Rep3']].max().max() # max value of the table
  df = df.sort_values(by = 'max_empai', ascending = False)
  for i,rep in enumerate(['Rep1', 'Rep2', 'Rep3']):
    ax.scatter(df.index, np.log2(df[rep]), label="Protein test" if i == 0 else "", color='royalblue', alpha=0.3, marker = 'o', s=40)
  plt.title(prot1[0]+' in '+prot1[1]+' with emPAI = '+str(threshold)+' as threshold ')
  if 'A' in control:
    for i,rep in enumerate(['CtrA1', 'CtrA2', 'CtrA3']):
      ax.scatter(df.index, np.log2(df[rep]), label="CtrA (without SPA tag)" if i == 0 else "", color='red', alpha=0.6, marker = 0, s=40)
    if 'CtrA4' in df :
      ax.scatter(df.index, np.log2(df.CtrA4), color='red', alpha=0.6, marker = 0, s=40)
  if 'C' in control:
    for i,rep in enumerate(['CtrC1', 'CtrC2', 'CtrC3']):
      ax.scatter(df.index, np.log2(df[rep]), label="CtrC (with SPA tag)" if i == 0 else "", color='yellowgreen', alpha=0.6, marker = 1, s=40)

  dftrue = df[df.C_is == True] 
  abs_ctrC =[] # triangle if absent in control
  abs_ctrA = []
  sigA = []
  sigC = [] # get index list for stars in plot
  prot_sig_C = []
  prot_sig_A = [] # get list of proteins where test is significant. 
  for i, my_index in enumerate(df.index):
    prot_sig_C.append(my_index)
    if (df.CtrC1.iloc[i] == threshold) and (df.CtrC2.iloc[i] == threshold) and (df.CtrC3.iloc[i] == threshold):
      abs_ctrC.append(i)
    elif df.C_is.iloc[i] == True:
      sigC.append(i)
    else: 
      prot_sig_C.pop()

    prot_sig_A.append(my_index)
    if 'CtrA4' in df :
      if (df.CtrA1.iloc[i] == threshold and df.CtrA2.iloc[i] == threshold and df.CtrA3.iloc[i] == threshold and df.CtrA4.iloc[i] == threshold):
        abs_ctrA.append(i)
      elif df.A_is.iloc[i] == True:
        sigA.append(i)
      else: 
        prot_sig_A.pop()
    else:
      if (df.CtrA1.iloc[i] == threshold and df.CtrA2.iloc[i] == threshold and df.CtrA3.iloc[i] == threshold):
        abs_ctrA.append(i)
      elif df.A_is.iloc[i] == True:
        sigA.append(i)
      else: 
        prot_sig_A.pop()

  abs_ctrC = [x+0.2 for x in abs_ctrC];  abs_ctrA = [x-0.2 for x in abs_ctrA]
  sigC = [x+0.2 for x in sigC];  sigA = [x-0.2 for x in sigA]
  prot_sig = list(set(prot_sig_A) | set(prot_sig_C)) # list of names of proteins.
  prot_sig.sort()
  print('all signif :', len(prot_sig))
  c = 0; nc = 1
  for i in prot_sig:
    if i in contaminant_genes:
      c += 1
#      print(i, end = ', ')
  print()
  print('contaminant :', c)
  for i in prot_sig:
    if i not in contaminant_genes:
      nc += 1
#      print(i, end =', ')
  print()
  print('not contaminant :', nc)

  ax.scatter(sigA, [np.log2(minval)-0.4]*len(sigA),c='red', marker=(5, 2), label = 'Significant test with CtrA', s=30) # add stars for Significant controls.
  ax.scatter(sigC, [np.log2(minval)-0.4]*len(sigC),c='yellowgreen', marker=(5, 2), label = 'Significant test with CtrC', s=30) # add stars for Significant controls.

  ax.scatter(abs_ctrA, [np.log2(minval)-0.4]*len(abs_ctrA),c='red', marker='^', label = 'Protein absent of CtrA', s = 20) # add triangles for proteins absent of each replicate of control A.
  ax.scatter(abs_ctrC, [np.log2(minval)-0.4]*len(abs_ctrC),c='yellowgreen', marker='^', label = 'Protein absent of CtrC', s= 20) # add triangles for proteins absent of each replicate of control C.

#  df_rep = np.log2(df[['Rep1', 'Rep2', 'Rep3']])
#  mean_conf_int = st.mean_confidence_interval(df_rep, 0.95, st.get_global_variance(prot1, threshold))
#  mean_conf_int = mean_conf_int.reindex(index = df.index)
#  ax.fill_between(mean_conf_int.index, mean_conf_int['conf_inf'], mean_conf_int['conf_sup'], color='b', alpha=.1)

  fig.canvas.draw()
  plt.ylim(np.log2(minval)-0.6, np.log2(maxval)+0.5)
  plt.xlabel('Gene name')
  plt.ylabel('log2(emPAI) value')
  plt.grid(axis = 'x') # vertical lines
  plt.xticks(rotation=90)

  for ticklabel in ax.get_xticklabels(): # adjust legend.
    if ticklabel.get_text() in contaminant_genes: # contaminant_genes est la liste de genes contaminants
      ticklabel.set_color('orange')
    if ticklabel.get_text() in prey[prot1[0]]:
      ticklabel.set_color('blue')
    if ticklabel.get_text() in interesting_prey[prot1[0]]:
      ticklabel.set_color('cyan')
  handles, labels = ax.get_legend_handles_labels()
  cont_patch = mpatches.Patch(color='orange', label='potential contaminant proteins')
  certain_interactor_patch = mpatches.Patch(color='blue', label='confirmed interactor proteins')
  interesting_interactor_patch = mpatches.Patch(color='cyan', label='interesting interactor proteins')
  handles.extend([cont_patch, certain_interactor_patch, interesting_interactor_patch]) # add to legend
  plt.legend(handles=handles)
  path_batch = "../Images/emPAI/log2values/"
# get an appropriate plot and saved image.
  manager = plt.get_current_fig_manager() # get full screen
  manager.window.showMaximized() # get full screen
  fig.tight_layout()
  fig.subplots_adjust(left=.05, bottom=.2, right=.96, top=.93) # marges
  filename = path_batch+prot1[0]+'_'+prot1[1][:6].replace('/', '_')+'_'+str(threshold)+'_pval_.05_log2values.png'
#  plt.savefig(path_batch+'test.svg') # image vectorisée
  plt.savefig(filename, transparent=False, dpi = 300) # image pixelisée, dpi = résolution
#  plt.show()
