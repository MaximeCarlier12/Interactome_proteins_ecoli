from globvar import *
import Load_PPI_screen as dt
import Data_processing as dp
import handle as hd
import qualitatif_stats as st
import quant as qt
import matplotlib.patches as mpatches
from matplotlib_venn import venn3, venn2
import glob

plt.rcParams.update(params)
used_prot_tuples = dp.good_proteins()

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
  '''Venn diagram of replicate of protein test for all used proteins. '''
  venn_data_type = {0:'emPAI', 1:'raw int', 2:'LFQ int'}
  plt.suptitle('Venn diagrams in replicates of a protein test')
  j = 0
  for (i,prot) in enumerate(used_prot_tuples):
    prot_batches = dt.get_batches(pd_samples, prot[0], prot[1])
    if data_type == 0:
      rep = [hd.load_df_unique_gene(i) for i in prot_batches]
    elif data_type == 1:
      prot_batches = qt.get_batches_missing_files(prot_batches)
      rep = [qt.select_intensity(qt.load_df_maxQ(i, 0), LFQ = False) for i in prot_batches]
    elif data_type == 2:
      prot_batches = qt.get_batches_missing_files(prot_batches)
      rep = [qt.select_intensity(qt.load_df_maxQ(i, 0), LFQ = True) for i in prot_batches]
    if i == 9 or i == 18:
        j += 9
        manager = plt.get_current_fig_manager() # get full screen
        manager.window.showMaximized() # get full screen
        plt.show()
        plt.suptitle('Venn diagrams in replicates of a protein test')
    plt.subplot(3,3, i+1-j)
    plt.title(prot[0]+' in '+prot[1])
    venn_diagram(rep, prot_batches)
  manager = plt.get_current_fig_manager() # get full screen
  manager.window.showMaximized() # get full screen
  plt.show()

def venn_two_rep(used_prot_tuples, data_type):
  '''Venn diagram of replicate of protein test for all used proteins.'''
  venn_data_type = {0:'emPAI', 1:'raw int', 2:'LFQ int'}
  plt.suptitle('Venn diagrams with only two replicates of a protein test')
  j = 0
  for (i,prot) in enumerate(used_prot_tuples):
    prot_batches = dt.get_batches(pd_samples, prot[0], prot[1])
    if data_type == 0:
      rep = [hd.load_df_unique_gene(i) for i in prot_batches]
    elif data_type == 1:
      prot_batches = qt.get_batches_missing_files(prot_batches)
      rep = [qt.select_intensity(qt.load_df_maxQ(i, 0), LFQ = False) for i in prot_batches]
    elif data_type == 2:
      prot_batches = qt.get_batches_missing_files(prot_batches)
      rep = [qt.select_intensity(qt.load_df_maxQ(i, 0), LFQ = True) for i in prot_batches]
    if i == 9 or i == 18:
        j += 9
        plt.show()
        plt.suptitle('Venn diagrams in replicates of a protein test')
    plt.subplot(3,3, i+1-j)
    plt.title(prot[0]+' in '+prot[1])
    venn_diagram(rep, prot_batches)
  plt.show()

def venn_ctr(used_prot_tuples, controls_dico, control_type, data_type):
  '''Venn diagram of replicates of a given control for all used proteins. We check if some files can't be used (e.g. A10)'''
  venn_data_type = {0:'emPAI', 1:'raw int', 2:'LFQ int'}
  plt.suptitle('Venn diagrams for replicates in controls type'+control_type)
  for i, cond in enumerate(controls_dico.keys()):
    if data_type == 0:
      ctr = [hd.load_df_unique_gene(bname) for bname in controls_dico[cond]]
    elif data_type == 1:
      ctr = [qt.select_intensity(qt.load_df_maxQ(bname, 0), LFQ = False) for bname in qt.get_batches_missing_files(controls_dico[cond])]
    elif data_type == 2:
      ctr = [qt.select_intensity(qt.load_df_maxQ(bname, 0), LFQ = True) for bname in qt.get_batches_missing_files(controls_dico[cond])]
    plt.subplot(1,3, i+1)
    plt.title('Condition '+cond)
    venn_diagram(ctr, controls_dico[cond])
  plt.show()

def venn_inter(used_prot_tuples, controls_dico):
  '''Venn diagram between intersect control and test.'''
  plt.suptitle('Venn diagrams of intersection of controls and replicates')
  for (i,prot) in enumerate(used_prot_tuples):
    prot_batches = dt.get_batches(pd_samples, prot[0], prot[1])
    rep = [hd.load_df_unique_gene(i) for i in prot_batches]
    interR = df_intersect(rep)
    ctr = [hd.load_df_unique_gene(i) for i in controls_dico[prot[1]]]
    interC = df_intersect(ctr)
    plt.subplot(3,3, i+1)
    plt.title(prot[0]+' in '+prot[1])
    venn_diagram([interR, interC])
  plt.show()

prot_3_reps = qt.prot_three_rep()
# venn_three_rep(qt.prot_three_rep(), 1)
# venn_two_rep(qt.prot_two_rep(), 2)
# venn_ctr(prot_3_reps, controls_typeA, 'A', 2)
# venn_ctr(prot_3_reps, controls_typeC, 'C', 2)

def first_version_sum_log_abundance(used_prot_tuples, data):
  '''Print sum(log(abundance)) per replicate and protein studied. 
  Different types of abundance : data = 0 : emPAI, data = 1 : raw intensity, data = 2 : LFQ intensity, data = 3 : normalized raw intensity.'''
  venn_data_type = {0:'emPAI', 1:'raw int', 2:'LFQ int'}
  fig,ax = plt.subplots()
  prots = []
  name_prots = []
  batch_names = []
  batch_namesA = []
  batch_namesC = []
  for prot in used_prot_tuples:
    prot_batches = dt.get_batches(pd_samples, prot[0], prot[1])
    batchA = dp.get_control_batches(pd_controls, 'MG1655 (TYPE A)' , prot[1])
    batchC = dp.get_control_batches(pd_controls, 'MG1655 (placI)mVenus-SPA-pUC19 (TYPE C2)' , prot[1])
    if data == 0 :
      df = hd.load_df_table(prot, True)
      sums = df[['Rep1', 'Rep2', 'Rep3', 'CtrC1', 'CtrC2', 'CtrC3', 'CtrA1', 'CtrA2', 'CtrA3']].sum(axis=0)
      plt.ylabel('sum(emPAI)')
      plt.yscale('symlog')
    elif data ==1 :
      prot_batches = qt.get_batches_missing_files(prot_batches)
      batch_names.append(prot_batches)
      batchA = qt.get_batches_missing_files(batchA)
      batch_namesA.append(batchA)
      batchC = qt.get_batches_missing_files(batchC)
      batch_namesC.append(batchC)
      df = hd.load_df_table_maxQ(prot, False, 0)
      plt.ylabel('log10(sum(raw intensity))')
      if 'CtrA4' in df.columns:
        sums = np.log10(df[['Rep1', 'Rep2', 'Rep3', 'CtrC1', 'CtrC2', 'CtrC3', 'CtrA1', 'CtrA2', 'CtrA3', 'CtrA4']].sum(axis=0))
      else:
        sums = np.log10(df[['Rep1', 'Rep2', 'Rep3', 'CtrC1', 'CtrC2', 'CtrC3', 'CtrA1', 'CtrA2', 'CtrA3']].sum(axis=0))
    elif data == 2:
      df = hd.load_df_table_maxQ(prot, True, 0)
      plt.ylabel('log10(sum(LFQ intensity))')
      sums = np.log10(df[['Rep1', 'Rep2', 'Rep3', 'CtrC1', 'CtrC2', 'CtrC3', 'CtrA1', 'CtrA2', 'CtrA3']].sum(axis=0))
    elif data == 3 :
      df = hd.load_df_table_maxQ(prot, False, 1)
      plt.ylabel('sum(normalized median raw intensity))')
      sums = np.log10(df[['Rep1', 'Rep2', 'Rep3', 'CtrC1', 'CtrC2', 'CtrC3', 'CtrA1', 'CtrA2', 'CtrA3']].sum(axis=0))
      # sums = df[['Rep1_log10', 'Rep2_log10', 'Rep3_log10', 'CtrC1_log10', 'CtrC2_log10', 'CtrC3_log10', 'CtrA1_log10', 'CtrA2_log10', 'CtrA3_log10']].sum(axis=0)
    else : print('error in data variable')
    name = prot[0]+'_'+prot[1][:6]
    name_prots.append(name)
    prots.append(sums)
#    print(sums[0])
  for i in range(len(prots)):
    ax.scatter([name_prots[i]]*3, prots[i][:3] , label="Protein test" if i == 0 else "", color='royalblue')
    ax.scatter([name_prots[i]]*3, prots[i][3:6] , label="CtrC" if i == 0 else "", color='chartreuse')
    print('here', prots[i][6:])
    print('df', df.columns)
    print(name_prots[i])
    if len(batch_namesA[i]) == 4:
      ax.scatter([name_prots[i]]*4, prots[i][6:] , label="CtrA" if i == 0 else "", color='red')
    elif len(batch_namesA[i]) == 3:
      ax.scatter([name_prots[i]]*3, prots[i][6:] , label="CtrA" if i == 0 else "", color='red')
    else : print('error in length batch_namesA') # simple check
    print(batch_names[i])
    print(prots[i][:3])
    for j,bn in enumerate(batch_names[i]):
      ax.annotate(bn, (name_prots[i], prots[i][j]), color = 'royalblue')
    for j,bn in enumerate(batch_namesC[i]):
      ax.annotate(bn, (name_prots[i], prots[i][3+j]), color = 'chartreuse')
    for j,bn in enumerate(batch_namesA[i]):
      print(bn)
      print('j', j)
      print(name_prots[i])
      print(len(prots[i]))
      print(prots[i][j+6])
      ax.annotate(bn, (name_prots[i], prots[i][6+j]), color = 'red')
  plt.legend()
  plt.title('Differences of abundance between files')
  plt.xlabel('Protein studied')
  plt.xticks(rotation=90)
  plt.grid(axis = 'x') # vertical lines
  manager = plt.get_current_fig_manager()
  manager.window.showMaximized()
  plt.tight_layout()
  plt.show()

def sum_log_abundance(used_prot_tuples, data):
  '''Print sum(log(abundance)) per replicate and protein studied. 
  Different types of abundance : data = 0 : emPAI, data = 1 : raw intensity, data = 2 : LFQ intensity, data = 3 : normalized median raw intensity, data = 4 : normalized bait raw intensity.'''
  fig,ax = plt.subplots()
  prots = []
  name_prots = []
  batch_names = []
  batch_namesA = []
  batch_namesC = []
  for prot in used_prot_tuples:
    name = prot[0]+'_'+prot[1][:6]
    name_prots.append(name)
    prot_batches = dt.get_batches(pd_samples, prot[0], prot[1])
    prot_batches = qt.get_batches_missing_files(prot_batches)
    if (data == 4) and ('U5' in prot_batches): # SeqA is absent from U5 file so we can't normalize this file.
      prot_batches.remove('U5')
    batch_names.append(prot_batches)
    # all_batches = prot_batches+batchC+batchA
    sums = []
    for i in prot_batches:
      if data == 1: #raw
        df = qt.select_intensity(qt.load_df_maxQ(i, 0), LFQ = False)
        nb_genes = len(df.index)
        avg_int = df[['Intensity_sample']].sum()/nb_genes
      elif data ==2: # lfq
        df = qt.select_intensity(qt.load_df_maxQ(i, 0), LFQ = True)
        nb_genes = len(df.index)
        avg_int = df[['LFQ_intensity_sample']].sum()/nb_genes
      elif data==3: # median
        print('number of proteins :', end = ' ')
        df = qt.select_intensity(qt.load_df_maxQ(i, 1), LFQ = False)
        nb_genes = len(df.index)
        print(i+' : '+ str(nb_genes), end = '\t')
        avg_int = df[['Normalized_intensity']].sum()/nb_genes
      elif data==4: # bait
        print('number of proteins :', end = ' ')
        df = qt.select_intensity(qt.load_df_maxQ(i, 2), LFQ = False)
        nb_genes = len(df.index)
        print(i+' : '+ str(nb_genes), end = '\t')
        avg_int = df[['Normalized_intensity']].sum()/nb_genes
      elif data==5: # q1
        print('number of proteins :', end = ' ')
        df = qt.select_intensity(qt.load_df_maxQ(i, 3), LFQ = False)
        nb_genes = len(df.index)
        print(i+' : '+ str(nb_genes), end = '\t')
        avg_int = df[['Normalized_intensity']].sum()/nb_genes
      elif data==6: # q3
        print('number of proteins :', end = ' ')
        df = qt.select_intensity(qt.load_df_maxQ(i, 4), LFQ = False)
        nb_genes = len(df.index)
        print(i+' : '+ str(nb_genes), end = '\t')
        avg_int = df[['Normalized_intensity']].sum()/nb_genes
      sums.append(np.log10(avg_int))
    print()
    prots.append(sums)
  CONDITION = ['LB log', 'LB O/N' ,'M9 0.2% ac O/N']
  for cond in CONDITION:
    name_prots.append('ctr_'+cond[:6])
    batchA = dp.get_control_batches(pd_controls, 'MG1655 (TYPE A)' , cond)
    batchA = qt.get_batches_missing_files(batchA)
    batchC = dp.get_control_batches(pd_controls, 'MG1655 (placI)mVenus-SPA-pUC19 (TYPE C2)' , cond)
    batchC = qt.get_batches_missing_files(batchC)
    batch_names.append(batchC+batchA)
    sums = []
    print('number of proteins :', end = ' ')
    for i in batchC+batchA :
      if data ==1:
        df = qt.select_intensity(qt.load_df_maxQ(i, 0), LFQ = False)
        sums.append(np.log10(df[['Intensity_sample']].sum()/len(df.index)))
      elif data ==2:
        df = qt.select_intensity(qt.load_df_maxQ(i, 0), LFQ = True)
        sums.append(np.log10(df[['LFQ_intensity_sample']].sum()/len(df.index)))
      elif data ==3:
        df = qt.select_intensity(qt.load_df_maxQ(i, 1), LFQ = False)
        print(i+' : '+ str(len(df['Normalized_intensity'])), end = '\t')
        sums.append(np.log10(df[['Normalized_intensity']].sum()/len(df.index)))
      elif data ==4:
        df = qt.select_intensity(qt.load_df_maxQ(i, 2), LFQ = False)
        print(i+' : '+ str(len(df['Normalized_intensity'])), end = '\t')
        sums.append(np.log10(df[['Normalized_intensity']].sum()/len(df.index)))
      elif data ==5:
        df = qt.select_intensity(qt.load_df_maxQ(i, 3), LFQ = False)
        print(i+' : '+ str(len(df['Normalized_intensity'])), end = '\t')
        sums.append(np.log10(df[['Normalized_intensity']].sum()/len(df.index)))
      elif data ==6:
        df = qt.select_intensity(qt.load_df_maxQ(i, 4), LFQ = False)
        print(i+' : '+ str(len(df['Normalized_intensity'])), end = '\t')
        sums.append(np.log10(df[['Normalized_intensity']].sum()/len(df.index)))
    print()
    prots.append(sums)
  # print('prots', prots)

  for i in range(len(prots)): # scatter and anotate
    ax.scatter([name_prots[i]]*len(prots[i]), prots[i][:] , label="batch used" if i == 0 else "", color='royalblue')
    # ax.scatter([name_prots[i]]*3, prots[i][3:6] , label="CtrC" if i == 0 else "", color='chartreuse')
    if 'ctr_' in name_prots[i]:
      for j,bn in enumerate(batch_names[i]):
        if j <3:
          if j%2 == 0:
            ax.annotate(bn, (name_prots[i], prots[i][j]), color = 'green', xytext=(-4, 0), textcoords='offset points', horizontalalignment='right', verticalalignment='center')
          else:
            ax.annotate(bn, (name_prots[i], prots[i][j]), color = 'green', xytext=(4, 0), textcoords='offset points', horizontalalignment='left', verticalalignment='center')
        else:
          if j%2 == 0:
            ax.annotate(bn, (name_prots[i], prots[i][j]), color = 'red', xytext=(-4, 0), textcoords='offset points', horizontalalignment='right', verticalalignment='center')
          else:
            ax.annotate(bn, (name_prots[i], prots[i][j]), color = 'red', xytext=(4, 0), textcoords='offset points', horizontalalignment='left', verticalalignment='center')
    else:
      for j,bn in enumerate(batch_names[i]):
          if j%2 == 0:
            ax.annotate(bn, (name_prots[i], prots[i][j]), color = 'royalblue', xytext=(-4, 0), textcoords='offset points', horizontalalignment='right', verticalalignment='center')
          else:
            ax.annotate(bn, (name_prots[i], prots[i][j]), color = 'royalblue', xytext=(4, 0), textcoords='offset points', horizontalalignment='left', verticalalignment='center')
  plt.legend()
  plt.ylabel('log10(mean(intensity))')
  plt.title('Differences of abundance between files')
  plt.xlabel('Protein studied')
  plt.xticks(rotation=90)
  plt.grid(axis = 'x') # vertical lines
  fig.tight_layout()
  manager = plt.get_current_fig_manager()
  manager.window.showMaximized()
  plt.show()

# sum_log_abundance(used_prot_tuples, 0)
# sum_log_abundance(prot_3_reps, 1) #raw
# sum_log_abundance(prot_3_reps, 2) #lfq
# sum_log_abundance(prot_3_reps, 3) #med
# sum_log_abundance(prot_3_reps, 4) #bait
# sum_log_abundance(prot_3_reps, 5) #q1
# sum_log_abundance(prot_3_reps, 6) #q3

def plot_emPAI(prot, control = 'AC'):
  '''Plot emPAI value for each gene for controls and test in log scale. It is not very visual.'''
  df = load_df_table(prot, True)
  empRep = []; empCtr = []
  df['max_empai'] = df[['Rep1', 'Rep2', 'Rep3']].max(axis=1)
  df = df.sort_values(by = 'max_empai', ascending = False)
  maxval = max(df.max(axis = 1)) # max value of the table
  minval = min(df.min(axis = 1)) # max value of the table
  #  empRep.append(np.mean([row.Rep1, row.Rep2, row.Rep3]))
  #  empCtr.append(np.mean([row.Ctr1, row.Ctr2, row.Ctr3]))
  plt.scatter(df.index, df.Rep1, label="Protein test", color='royalblue')
  plt.scatter(df.index, df.Rep2, color='royalblue')
  plt.scatter(df.index, df.Rep3, color='royalblue')
  plt.title(prot[0]+' in '+prot[1])
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

def plot_log2_emPAI(prot, threshold, contaminant_genes):
  '''Plot log2(emPAI) value for each gene for controls and test.'''

  fig,ax = plt.subplots()
  width = 12.5 ; height = 6.35 # taille finale de ta figure png
#  width = 14.4 ; height = 7.15 # taille finale de ta figure svg
  fig.set_size_inches(width, height)
  var_test = hd.load_df_equal_test()
  df = hd.load_df_table(prot, True)
  df['max_empai'] = df[['Rep1', 'Rep2', 'Rep3']].max(axis=1)
  minval = df[['Rep1', 'Rep2', 'Rep3']].min().min() # min value of the table
  maxval = df[['Rep1', 'Rep2', 'Rep3']].max().max() # max value of the table
  df = df.sort_values(by = 'max_empai', ascending = False)
  for i,rep in enumerate(['Rep1', 'Rep2', 'Rep3']):
    ax.scatter(df.index, np.log2(df[rep]), label="Protein test" if i == 0 else "", color='royalblue', alpha=0.3, marker = 'o', s=40)
  plt.title(prot[0]+' in '+prot[1]+' with emPAI = '+str(threshold)+' as threshold ')
  for i,rep in enumerate(['CtrA1', 'CtrA2', 'CtrA3']):
    ax.scatter(df.index, np.log2(df[rep]), label="CtrA (without SPA tag)" if i == 0 else "", color='red', alpha=0.6, marker = 0, s=40)
  if 'CtrA4' in df :
    ax.scatter(df.index, np.log2(df.CtrA4), color='red', alpha=0.6, marker = 0, s=40)
  for i,rep in enumerate(['CtrC1', 'CtrC2', 'CtrC3']):
    ax.scatter(df.index, np.log2(df[rep]), label="CtrC (with SPA tag)" if i == 0 else "", color='forestgreen', alpha=0.6, marker = 1, s=40)

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
  ax.scatter(sigC, [np.log2(minval)-0.4]*len(sigC),c='forestgreen', marker=(5, 2), label = 'Significant test with CtrC', s=30) # add stars for Significant controls.

  ax.scatter(abs_ctrA, [np.log2(minval)-0.4]*len(abs_ctrA),c='red', marker='^', label = 'Protein absent of CtrA', s = 20) # add triangles for proteins absent of each replicate of control A.
  ax.scatter(abs_ctrC, [np.log2(minval)-0.4]*len(abs_ctrC),c='forestgreen', marker='^', label = 'Protein absent of CtrC', s= 20) # add triangles for proteins absent of each replicate of control C.

  df_rep = np.log2(df[['Rep1', 'Rep2', 'Rep3']])
  mean_conf_int = st.mean_confidence_interval(df_rep, 0.95, st.get_global_variance(prot, threshold))
  mean_conf_int = mean_conf_int.reindex(index = df.index)
  print(mean_conf_int)  
  ax.plot( mean_conf_int['mean'], '--', linewidth=0.7, color = 'royalblue', alpha = 0.5)
  ax.fill_between(mean_conf_int.index, mean_conf_int['conf_inf'], mean_conf_int['conf_sup'], color='royalblue', alpha=.08)

  df_ctr = np.log2(df[['CtrC1', 'CtrC2', 'CtrC3']])
  mean_conf_int = st.mean_confidence_interval(df_ctr, 0.95, st.get_global_variance(prot, threshold))
  mean_conf_int = mean_conf_int.reindex(index = df.index)
  print(mean_conf_int)
  ax.plot( mean_conf_int['mean'], '--', linewidth=0.7, color = 'forestgreen', alpha = 0.5)
  ax.fill_between(mean_conf_int.index, mean_conf_int['conf_inf'], mean_conf_int['conf_sup'], color='forestgreen', alpha=.08)

  if 'CtrA4' in df:
    df_ctr = np.log2(df[['CtrA1', 'CtrA2', 'CtrA3', 'CtrA4']])
  else:
    df_ctr = np.log2(df[['CtrA1', 'CtrA2', 'CtrA3']])
  mean_conf_int = st.mean_confidence_interval(df_ctr, 0.95, st.get_global_variance(prot, threshold))
  mean_conf_int = mean_conf_int.reindex(index = df.index)
  print(mean_conf_int)
  ax.plot( mean_conf_int['mean'], '--', linewidth=0.7, color = 'red', alpha = 0.3)
  ax.fill_between(mean_conf_int.index, mean_conf_int['conf_inf'], mean_conf_int['conf_sup'], color='red', alpha=.06)

  fig.canvas.draw()
  plt.ylim(np.log2(minval)-0.6, np.log2(maxval)+0.5)
  plt.xlabel('Gene name')
  plt.ylabel('log2(emPAI) value')
  plt.grid(axis = 'x') # vertical lines
  plt.xticks(rotation=90)

  for ticklabel in ax.get_xticklabels(): # adjust legend.
    if ticklabel.get_text() in contaminant_genes: # contaminant_genes est la liste de genes contaminants
      ticklabel.set_color('orange')
    if ticklabel.get_text() in prey[prot[0]]:
      ticklabel.set_color('blue')
    if ticklabel.get_text() in interesting_prey[prot[0]]:
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
  filename = path_batch+prot[0]+'_'+prot[1][:6].replace('/', '_')+'_'+str(threshold)+'_pval_.05_log2values.png'
#  plt.savefig(path_batch+'test.svg') # image vectorisée
  plt.savefig(filename, transparent=False, dpi = 300) # image pixelisée, dpi = résolution
  # plt.show()
