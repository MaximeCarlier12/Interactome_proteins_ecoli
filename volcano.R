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
                x = 'log2FC_C', y = 'pvalC_corr', 
                pCutoff = 0.05, FCcutoff = 1.5, 
                xlim = c(-6,10), ylim = c(-1,15), 
                title = 'Replicates versus green control CtrC (tagged)') %xlab, ylab ? legend ? 
EnhancedVolcano(data, lab = rownames(data), 
                x = 'log2FC_A', y = 'pvalA_corr', 
                pCutoff = 0.05, FCcutoff = 1.5, 
                xlim = c(-6,8), ylim = c(-1,17), 
                title = 'Replicates versus red control CtrA (untagged)')
