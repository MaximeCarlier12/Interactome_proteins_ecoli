rm(list=ls())
setwd("~/Documents/Stage/Interactome_proteins_ecoli")
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')
BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)
library(readr)
PROTEIN = c('DnaA', 'DiaA', 'Hda', 'SeqA', 'HolD', 'DnaB', 'DnaG', 'NrdB')
CONDITION = c('LB log', 'LB O_N' ,'M9 0.2')
path = 'maxQ/New_data/Protein_table/'
for (p in PROTEIN){
  for (c in CONDITION){
    filename = paste0(path,p,'_',c,'_Int_Norm_Med.csv')
    data = read.csv(file = filename, row.names = 1)
    new_file1 = paste0('maxQ/New_data/Figures/Volcanos/Volcano_', p,'_',c, '_ctrC.svg')
    Volcano1 = EnhancedVolcano(data, lab = rownames(data), 
                        x = 'log2FC_C', y = 'adj_pvalC', 
                        pCutoff = 0.05, FCcutoff = 1.5, 
                        xlim = c(min(data$log2FC_C)-1,max(data$log2FC_C)+1), ylim = c(-1, -log10(min(data$adj_pvalC))+1), 
                        legendLabels = c('NS', expression(Log[2]~FC), 'adjusted p-value', expression(p-value~and~log[2]~FC)),
                        title = 'Replicates versus green control CtrC (tagged)') #lab, ylab ? legend ? 
    svg(filename=new_file1, 
        width=9.4, 
        height=6.69, 
        pointsize=12)
    print(Volcano1)
    dev.off()
    new_file2 = paste0('maxQ/New_data/Figures/Volcanos/Volcano_', p,'_',c, '_ctrA.svg')
    Volcano2 = EnhancedVolcano(data, lab = rownames(data), 
                               x = 'log2FC_A', y = 'adj_pvalA', 
                               pCutoff = 0.05, FCcutoff = 1.5, 
                               xlim = c(min(data$log2FC_A)-1,max(data$log2FC_A)+1), ylim = c(-1, -log10(min(data$adj_pvalA))+1), 
                               legendLabels = c('NS', expression(Log[2]~FC), 'adjusted p-value', expression(p-value~and~log[2]~FC)),
                               title = 'Replicates versus red control CtrA (untagged)') #lab, ylab ? legend ? 
    svg(filename=new_file2, 
        width=9.4, 
        height=6.69, 
        pointsize=12)
    print(Volcano2)
    dev.off()
  }
}







filename = "maxQ/SEGREGATED-20200619T092017Z-001/Protein_table/HolD_LB log_Int_Norm_Med.csv"
data = read.csv(file = filename, row.names = 1)
#png(filename="maxQ/New_data/Figures/test.png", 
#    units="in", 
#    width=9.4, 
#    height=6.69, 
#    pointsize=12, 
#    res=72)
svg(filename="maxQ/New_data/Figures/Std_SVG.svg", 
    width=9.4, 
    height=6.69, 
    pointsize=12)
EnhancedVolcano(data, lab = rownames(data), 
                x = 'log2FC_C', y = 'adj_pvalC', 
                pCutoff = 0.05, FCcutoff = 1.5, 
                xlim = c(-6,10), ylim = c(-1,15), 
                legendLabels = c('NS', expression(Log[2]~FC), 'adjusted p-value', expression(p-value~and~log[2]~FC)),
                title = 'Replicates versus green control CtrC (tagged)') #lab, ylab ? legend ? 
dev.off()
EnhancedVolcano(data, lab = rownames(data), 
                x = 'log2FC_A', y = 'adj_pvalA', 
                pCutoff = 0.05, FCcutoff = 1.5, 
                xlim = c(-6,8), ylim = c(-1,17), 
                legendLabels = c('NS', expression(Log[2]~FC), 'adjusted p-value', expression(p-value~and~log[2]~FC)),
                title = 'Replicates versus red control CtrA (untagged)')
                
