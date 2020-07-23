from globvar import *
import Load_PPI_screen as dt
import Data_processing as dp
import handle as hd
import visualization as vz
import qualitatif_stats as st
import glob

used_prot_tuples = dp.good_proteins()

hd.create_table(used_prot_tuples)
threshold = hd.create_condensed_table(used_prot_tuples, 0.25)
for prot1 in used_prot_tuples[0:1]:
  dt.header(prot1[0]+' in '+prot1[1])
  st.log_calculation(prot1)
  st.test_normal_equal_var(prot1, threshold)
  vz.plot_log2_emPAI(prot1, threshold, contaminant_genes)
#for prot1 in used_prot_tuples:
#  df = hd.load_df_table(prot1, True)
#  print(prot1[0], prot1[1][:6]+'\t', 'TestC', df[df.C_is ==True].shape[0],'\t', 'TestA', df[df.A_is ==True].shape[0])

plt.close('all')
