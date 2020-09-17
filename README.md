# Interactome_proteins_ecoli
Characterization of the interactome of proteins involved in DNA replication in E.coli, in collaboration with the poland Monica Glinkowska's lab.
This code is intended to make an bioinformatic analysis of MS/MS data. The goal is to identify the interaction between proteins in several condition growths and for different bait proteins (that have a role in DNA replication).

FILES : 

environment.yml : file to use in order to have exactly the same package as me so that you are sure that the codes can be used without version package problems. 

globvar.py : variables that are similar to files such as constants.

Load_ppi_screen.py : functions that load data from the main files sent by the poland lab (sample codes, controls, batches...)

Data_processing.py : treatment on emPAI data that allows analysis (get gene name, appropriate replicates to use...)

qualitatif_stats.py : functions specific for statistical analysis ()

visualization.py : functions specific for plots (all venn diagrams, big plots for emPAI data and first visual plots to know our data better)

quant.py : 

main.py : example to follow if you want to use this code. One example with emPAI (not the best, quite limited) and another with maxQuant (with the new data that is to say that each condition has 3 replicates). 

volcano.R : code to implemente volcano plots for each protein-condition

INSTRUCTIONS : 

To create an virtual environment from 'environment.yml', do : 
- conda env create -f environment.yml
Then, activate your environment with : 
- conda activate myenv
