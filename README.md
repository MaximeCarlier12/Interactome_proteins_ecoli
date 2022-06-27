Python code developped by MAXIME CARLIER for the article "Protein interaction network analysis reveals growth conditions-specific crosstalk between chromosomal DNA replication and other cellular processes in E. coli". The objective is to identify significant protein-protein interactions from mass spectrometry data.

Input files were obtained using either Mascot or MaxQuant softwares.

FILES :

environment.yml : requirement file (Python packages and versions)

globvar.py : global variables

load_ppi_screen.py : functions that load data from raw files

data_processing.py : functions for the treatment of emPAI data

qualitative_stats.py : functions for statistical analyses

visualisation.py : functions for visualisation (venn diagrams, plots...)

quant.py : functions for protein quantification

main.py : main functions with examples (emPAI and maxQuant data)

volcano.R : R code for creating volcano plots
