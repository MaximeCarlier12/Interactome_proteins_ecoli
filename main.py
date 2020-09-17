from globvar import *
import Load_PPI_screen as dt
import Data_processing as dp
import handle as hd
import visualization as vz
import qualitatif_stats as st
import quant as qt
import glob

# data 1 : emPAI
def emPAI_analysis(): # saved tables and plots
  create_csv_genes_good_proteins() # create files 'genes_batchname.csv' type. It has to be done only once.
  create_all_csv_unique_gene() #create files 'unique_gene_batchname.csv' type. It has to be done only once.
  # Be careful, C13 and A10 are in NCBInr database (controls_typeA in LB O/N and M9 0.2% ac O/N)
  used_prot_tuples = dp.good_proteins() # load proteins that have been used for emPAI (database problems)
  venn_three_rep(qt.prot_three_rep(), 0) # venn diagram for condition with 3 replicates
  hd.create_all_tables(used_prot_tuples) # creation of all tables that contain proteins values in replicates and controls.
  threshold = hd.create_all_condensed_tables(used_prot_tuples, 0.25) # creation of all condensed tables 
  # (i.e those with at least 2 replicates). 0.25 is the threshold for log2(emPAI) values.
  for prot in used_prot_tuples:
    dt.header(prot[0]+' in '+prot[1]) # simple display
    st.log_calculation(prot) # log2 calculations
    st.test_normal_equal_var(prot, threshold) # calculation of pvalues with global variance. Every result is saved in xlsx file.
    vz.plot_log2_emPAI(prot, threshold, contaminant_genes) # plot log2(emPAI) values = huge graph

def quant_analysis(LFQ_bool, normalization_method):
  '''Normalization method can be 0 (without norm), 1 (median norm), 2 (bait value norm), 3 (Q1 value) or 4 (Q3 value).
  LFQ_bool can be True (we use LFQ intensity) or False (we used raw intensity).'''
  prot_3_rep = qt.prot_three_rep() # get all proteins that contain 3 replicates that can be used with the current data. Returns tuple with bait protein and growth condition.
  # qt.create_all_norm_files(prot_3_rep) # Create all normalized files (q1, q3, median) for given proteins. It has to be done only once.
  vz.venn_three_rep(prot_3_rep, 1) # represents venn diagram between the 3 replicates, for each protein + growth condition.
  for prot in prot_3_rep:
    dt.header(prot[0]+ ' in '+prot[1]) # simple display
    threshold = qt.create_table(prot, LFQ = LFQ_bool, normalize = normalization_method) # creation of a table for each condition (bait+growth+media)
    #used i.e those with at least 2 replicates). 0.25 is the threshold for log2(emPAI) values. It saves the table in a new file.
    qt.prot_absent_controls(prot, threshold, LFQ_bool, normalization_method) # modify the table file and 
    # add columns 'Absent_ctrA' and 'Absent_ctrC' that are booleans and tell if the protein is present in at least one replicate of the control or not.
  common_variances = qt.get_common_variances_per_media(prot_3_rep, threshold, LFQ_bool, normalization_method) # calculation of the common variance (from every single protein+growth condition) for each media condition (list of 3 values).
  for prot in prot_3_rep: # loop for each condition used.
    qt.fold_change(prot, LFQ_bool, normalization_method) # calculate fold change for each protein and add it to file table.
    qt.test_normal_equal_var(prot, threshold, LFQ_bool, normalization_method, common_variances) # calculate tests for each protein and add it to file table.
    qt.plot_log10_abundance(prot, threshold, LFQ_bool, normalization_method, common_variances) # plot log10(normalized values)
    qt.create_entire_final_csv(prot, LFQ_bool, normalization_method) # create files that summarize results. It contains important information for all prey proteins in the 3 replicates such as pvalues, FC, putative type, mean between replicates.
    qt.create_significant_prot_final_csv(prot, LFQ_bool, normalization_method) # create files that summarize results, only for enriched proteins (present in 3 replicates AND absent or significantly enriched compared to one control). Same information as create_entire_final_csv
  for prot_name in PROTEINS:
    qt.create_putative_proteins_analysis(prot_3_rep, prot_name, LFQ_bool, normalization_method) # create 1 file per bait protein that contains the name of all enriched prey proteins. There is one sheet per growth condition (i.e 3 sheets).

# quant_analysis(False, 1)

plt.close('all')
