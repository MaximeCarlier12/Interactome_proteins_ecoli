import Load_PPI_screen as dt
import Data_processing as dp
import handle as hd
import csv
import numpy as np
import matplotlib.pyplot as plt
import qualitatif_stats as st
import quant as qt
import pandas as pd
import matplotlib.patches as mpatches
from matplotlib_venn import venn3, venn2
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
pd_samples = dt.load_based_screen_samples()
used_prot_tuples = dp.good_proteins()
pd_controls = dt.load_based_screen_controls()
controls_typeC = {'LB log':['S1','S2','S3'], 'LB O/N':['O9', 'R4', 'R5'], 'M9 0.2% ac O/N':['R1', 'R2', 'R3']}
controls_typeA = {'LB log':['L1', 'T7', 'T8', 'T9'], 'LB O/N':['C13', 'P5', 'U10'], 'M9 0.2% ac O/N':['A10', 'T5', 'T6']}
missing_files = ['A1', 'A2', 'M4', 'M5', 'A9', 'A10', 'I1', 'I2']

def venn_diagram(data, names):
  '''Plot a venn2 or venn3 diagram. data is a list of replicates.'''
  set_array = []
  # print(data[0])
  # print(data[0].columns)
  # print(data[0].index)
  for rep in data:
    if 'emPAI' in rep.columns:
      rep = rep.set_index(['Gene_name'])
    set_array.append(set(rep.index))
  if len(data) == 3:
    venn3(set_array, names)   # venn3 works for three sets
  elif len(data) == 2:
    venn2(set_array, names)   # venn3 works for three sets
  elif len(data) == 4:
    venn3(set_array[:3], names[:3])   # venn3 works for three sets
  else : print('error, please change data length')

def venn_three_rep(used_prot_tuples, data_type):
  '''Venn diagram of replicate of protein test for all used proteins.'''
  plt.suptitle('Venn diagrams in replicates of a protein test')
  j = 0
  for (i,prot1) in enumerate(used_prot_tuples):
    prot_batches = dt.get_batches(pd_samples, prot1[0], prot1[1])
    if data_type == 0:
      rep = [hd.load_df_unique_gene(i) for i in prot_batches]
    elif data_type == 1:
      prot_batches = qt.get_batches_missing_files(prot_batches)
      rep = [qt.select_intensity(qt.load_df(i), LFQ = False) for i in prot_batches]
    elif data_type == 2:
      prot_batches = qt.get_batches_missing_files(prot_batches)
      rep = [qt.select_intensity(qt.load_df(i), LFQ = True) for i in prot_batches]
    if i == 9 or i == 18:
        j += 9
        manager = plt.get_current_fig_manager() # get full screen
        manager.window.showMaximized() # get full screen
        plt.show()
        plt.suptitle('Venn diagrams in replicates of a protein test')
    plt.subplot(3,3, i+1-j)
    plt.title(prot1[0]+' in '+prot1[1])
    venn_diagram(rep, prot_batches)
  manager = plt.get_current_fig_manager() # get full screen
  manager.window.showMaximized() # get full screen
  plt.show()

def venn_two_rep(used_prot_tuples, data_type):
  '''Venn diagram of replicate of protein test for all used proteins.'''
  plt.suptitle('Venn diagrams with only two replicates of a protein test')
  j = 0
  for (i,prot1) in enumerate(used_prot_tuples):
    prot_batches = dt.get_batches(pd_samples, prot1[0], prot1[1])
    if data_type == 0:
      rep = [hd.load_df_unique_gene(i) for i in prot_batches]
    elif data_type == 1:
      prot_batches = qt.get_batches_missing_files(prot_batches)
      rep = [qt.select_intensity(qt.load_df(i), LFQ = False) for i in prot_batches]
    elif data_type == 2:
      prot_batches = qt.get_batches_missing_files(prot_batches)
      rep = [qt.select_intensity(qt.load_df(i), LFQ = True) for i in prot_batches]
    if i == 9 or i == 18:
        j += 9
        plt.show()
        plt.suptitle('Venn diagrams in replicates of a protein test')
    plt.subplot(3,3, i+1-j)
    plt.title(prot1[0]+' in '+prot1[1])
    venn_diagram(rep, prot_batches)
  plt.show()

def venn_ctr(used_prot_tuples, controls_dico, control_type, data_type):
  '''Venn diagram of replicates of a given control for all used proteins. We check if some files can't be used (e.g. A10)'''
  plt.suptitle('Venn diagrams for replicates in controls type'+control_type)
  for i, cond in enumerate(controls_dico.keys()):
    if data_type == 0:
      ctr = [hd.load_df_unique_gene(bname) for bname in controls_dico[cond]]
    elif data_type == 1:
      ctr = [qt.select_intensity(qt.load_df(bname), LFQ = False) for bname in qt.get_batches_missing_files(controls_dico[cond])]
    elif data_type == 2:
      ctr = [qt.select_intensity(qt.load_df(bname), LFQ = True) for bname in qt.get_batches_missing_files(controls_dico[cond])]
    plt.subplot(1,3, i+1)
    plt.title('Condition '+cond)
    venn_diagram(ctr, controls_dico[cond])
  plt.show()

def venn_inter(used_prot_tuples, controls_dico):
  '''Venn diagram between intersect control and test.'''
  plt.suptitle('Venn diagrams of intersection of controls and replicates')
  for (i,prot1) in enumerate(used_prot_tuples):
    prot_batches = dt.get_batches(pd_samples, prot1[0], prot1[1])
    rep = [hd.load_df_unique_gene(i) for i in prot_batches]
    interR = df_intersect(rep)
    ctr = [hd.load_df_unique_gene(i) for i in controls_dico[prot1[1]]]
    interC = df_intersect(ctr)
    plt.subplot(3,3, i+1)
    plt.title(prot1[0]+' in '+prot1[1])
    venn_diagram([interR, interC])
  plt.show()

prot_3_reps = qt.prot_three_rep()
# venn_three_rep(qt.prot_three_rep(), 2)
# venn_two_rep(qt.prot_two_rep(), 2)
# venn_ctr(prot_3_reps, controls_typeA, 'A', 2)
# venn_ctr(prot_3_reps, controls_typeC, 'C', 2)

def sum_log_abundance(used_prot_tuples, data):
  '''Print sum(log(abundance)) per replicate and protein studied. 
  Different types of abundance : data = 0 : emPAI, data = 1 : raw intensity, data = 2 : LFQ intensity.'''
  fig,ax = plt.subplots()
  prots = []
  name_prots = []
  for prot1 in used_prot_tuples:
    if data == 0 :
      df = hd.load_df_table(prot1, True)
      sums = np.log2(df[['Rep1', 'Rep2', 'Rep3', 'CtrC1', 'CtrC2', 'CtrC3', 'CtrA1', 'CtrA2', 'CtrA3']]).sum(axis=0)
      plt.ylabel('sum(emPAI)')
    elif data ==1 :
      df = qt.load_df_table(prot1, False)
      plt.ylabel('sum(raw intensity)')
      sums = np.log10(df[['Rep1', 'Rep2', 'Rep3', 'CtrC1', 'CtrC2', 'CtrC3', 'CtrA1', 'CtrA2', 'CtrA3']]).sum(axis=0)
    elif data == 2:
      df = qt.load_df_table(prot1, True)
      plt.ylabel('sum(LFQ intensity)')
      sums = np.log10(df[['Rep1', 'Rep2', 'Rep3', 'CtrC1', 'CtrC2', 'CtrC3', 'CtrA1', 'CtrA2', 'CtrA3']]).sum(axis=0)
    else : print('error in data variable')
    name = prot1[0]+'_'+prot1[1][:4]
    name_prots.append(name)
    prots.append(sums)
#    print(sums[0])
  for i in range(len(prots)):
    ax.scatter([name_prots[i]]*3, prots[i][:3] , label="Protein test" if i == 0 else "", color='royalblue')
    ax.scatter([name_prots[i]]*3, prots[i][3:6] , label="CtrC" if i == 0 else "", color='chartreuse')
    ax.scatter([name_prots[i]]*3, prots[i][6:] , label="CtrA" if i == 0 else "", color='red')
  plt.legend()
  plt.title('Differences of abundance between files')
  plt.xlabel('Protein studied')
  plt.xticks(rotation=90)
  plt.grid(axis = 'x') # vertical lines
  manager = plt.get_current_fig_manager()
  manager.window.showMaximized()
  plt.tight_layout()
  plt.show()

# sum_log_abundance(used_prot_tuples, 0)
sum_log_abundance(prot_3_reps, 1)
# sum_log_abundance(prot_3_reps, 2)

def plot_emPAI(prot1, control = 'AC'):
  '''Plot emPAI value for each gene for controls and test in log scale. It is not very visual.'''
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

def plot_log2_emPAI(prot1, threshold, contaminant_genes):
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
  for i,rep in enumerate(['CtrA1', 'CtrA2', 'CtrA3']):
    ax.scatter(df.index, np.log2(df[rep]), label="CtrA (without SPA tag)" if i == 0 else "", color='red', alpha=0.6, marker = 0, s=40)
  if 'CtrA4' in df :
    ax.scatter(df.index, np.log2(df.CtrA4), color='red', alpha=0.6, marker = 0, s=40)
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

  df_rep = np.log2(df[['Rep1', 'Rep2', 'Rep3']])
  mean_conf_int = st.mean_confidence_interval(df_rep, 0.95, st.get_global_variance(prot1, threshold))
  mean_conf_int = mean_conf_int.reindex(index = df.index)
  print(mean_conf_int)  
  ax.plot( mean_conf_int['mean'], '-', linewidth=1, color = 'royalblue', alpha = 0.5)
  ax.fill_between(mean_conf_int.index, mean_conf_int['conf_inf'], mean_conf_int['conf_sup'], color='royalblue', alpha=.15)

  df_ctr = np.log2(df[['CtrC1', 'CtrC2', 'CtrC3']])
  mean_conf_int = st.mean_confidence_interval(df_ctr, 0.95, st.get_global_variance(prot1, threshold))
  mean_conf_int = mean_conf_int.reindex(index = df.index)
  print(mean_conf_int)
  ax.plot( mean_conf_int['mean'], '-', linewidth=1, color = 'yellowgreen', alpha = 0.5)
  ax.fill_between(mean_conf_int.index, mean_conf_int['conf_inf'], mean_conf_int['conf_sup'], color='yellowgreen', alpha=.2)

  if 'CtrA4' in df:
    df_ctr = np.log2(df[['CtrA1', 'CtrA2', 'CtrA3', 'CtrA4']])
  else:
    df_ctr = np.log2(df[['CtrA1', 'CtrA2', 'CtrA3']])
  mean_conf_int = st.mean_confidence_interval(df_ctr, 0.95, st.get_global_variance(prot1, threshold))
  mean_conf_int = mean_conf_int.reindex(index = df.index)
  print(mean_conf_int)
  ax.plot( mean_conf_int['mean'], '-', linewidth=1, color = 'red', alpha = 0.3)
  ax.fill_between(mean_conf_int.index, mean_conf_int['conf_inf'], mean_conf_int['conf_sup'], color='red', alpha=.1)

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
