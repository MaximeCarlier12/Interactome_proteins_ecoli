rm(list=ls())
setwd("~/Documents/Stage/Interactome_proteins_ecoli")
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')
BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)
library(readr)
filename = 'maxQ/SEGREGATED-20200619T092017Z-001/Protein_table/HolD_LB log_Int_Norm_Med.csv'
data = read.csv(file = filename, row.names = 1)
EnhancedVolcano(data, lab = rownames(data), 
                x = 'log2FC_C', y = 'adj_pvalC', 
                pCutoff = 0.05, FCcutoff = 1.5, 
                xlim = c(-6,10), ylim = c(-1,15), 
                legendLabels = c('NS', expression(Log[2]~FC), 'adjusted p-value', expression(p-value~and~log[2]~FC)),
                title = 'Replicates versus green control CtrC (tagged)') #lab, ylab ? legend ? 
EnhancedVolcano(data, lab = rownames(data), 
                x = 'log2FC_A', y = 'adj_pvalA', 
                pCutoff = 0.05, FCcutoff = 1.5, 
                xlim = c(-6,8), ylim = c(-1,17), 
                legendLabels = c('NS', expression(Log[2]~FC), 'adjusted p-value', expression(p-value~and~log[2]~FC)),
                title = 'Replicates versus red control CtrA (untagged)')
                
