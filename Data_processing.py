from globvar import *
import Load_PPI_screen as dt
from Bio import Entrez
import glob

os.chdir("/home/carlier/Documents/Stage/Interactome_proteins_ecoli/") # where is the MS data

def intersect_index(replicates):
  ''' From a list of pd file replicates, gives a list of indexes present in all replicates.'''
  intersect = replicates[0].index
  for i in replicates : 
    intersect = intersect.intersection(i.index)
  return intersect

def get_control_batches(pd_control, type_control, condition):
  ''' From the pd.control file, get pd files of the different replicates available for a certain type of control and a certain type of condition.'''
  if not condition in CONDITION : 
    print("Error condition")
    return None
  if not type_control in CONTROL : 
    print("Error control")
    return None
  control_row = pd_control[pd_control.strain== type_control] # get only one row.
  #repeat_condition = 'repeat '+condition
  batches = control_row[condition].iloc[0]  # get batches of the different replicates.
  batches = batches.replace(';',',') # replace ';' by ',' because there are two separator ways in the initial file. 
  batches = batches.replace(" ", "") # remove spaces
  batches = batches.split(',') # separate the replicates
  # dt.header('batches')    
  # print(batches)
  return batches

def get_control_replicates(pd_control, type_control, condition):
  ''' From the pd.control file, get pd files of the different replicates available for a certain type of control and a certain type of condition.'''
  batches = get_control_batches(pd_control, type_control, condition)
  return [dt.load_df(batches[i]) for i in range(len(batches))]

def get_protein_replicates(pd_samples, protein, condition):
  '''Loads replicate files for a specific protein, with emPAI files. '''
  batches = dt.get_batches(pd_samples, protein, condition)
  return [dt.load_df(batches[i]) for i in range(len(batches))]

def get_replicates(pd, obj, condition, control=False):
  ''' Global function with a boolean indicating if we want replicates of a control or a protein.
Obj : control or protein.'''

  if control == True:
    return get_control_replicates(pd, obj,condition)
  else:
    return get_protein_replicates(pd, obj, condition)

def count_databases():
  '''For each file, print an error if databases are different, with a summary of usable files.'''

  dt.header("test database")
  nb_inter = []
  bd = []
  bd2 = []
  for p in PROTEINS :
    check = False;check_c = False
    for c in CONDITION : 
      check = False
      check_c = False
      rep = get_protein_replicates(pd_samples, p, c)
      inter = len(intersect_index(rep))
      if inter != 0:
        check = True
      nb_inter.append(inter)
      if check == False :
        print("Error for", p, "in condition", c) # for each protein, in there at least one control : no
      if check == True:
        bd.append(rep[0]['Database'][0])
        ctr = 'MG1655 (placI)mVenus-SPA-pUC19 (TYPE C2)'
        rep = get_control_replicates(pd_controls, ctr, c)
        bd2.append(rep[0]['Database'][0])
        inter = len(intersect_index(rep))
        if inter != 0:
          check_c = True
        else : print("error for control in condition", c)

  print(nb_inter, len(nb_inter))
  nb_inter1 = [nb_inter[i] for i in range(len(nb_inter)) if nb_inter[i] == 0]
  nb_inter2 = []
  print(bd)
  print(bd2)
  print([bd[i] == bd2[i] for i in range(len(bd))])
  for c in CONDITION : 
    check = False
    for ctr in CONTROL :
      rep = get_control_replicates(pd_controls, ctr, c)
      inter = len(intersect_index(rep))
      if inter != 0:
        check = True
      nb_inter2.append(inter)
    print(check) # for each condition, True if there is at least one control with 3 replicates from the same database.
  print('number of different genes for interesting proteins :', nb_inter2)
  nb_inter3 = [nb_inter2[i] for i in range(len(nb_inter2)) if nb_inter2[i] == 0]
  print("number of protein with different databases :", len(nb_inter1), "out of", len(nb_inter))
  print("number of controls with different databases :", len(nb_inter3), "out of", len(nb_inter2))

#count_databases()

def good_proteins():
  '''List of tuple containing the protein to be used and its condition (because of data problems we cannot use all files). 
  To do so, we check if the intersection of index replicates is empty or not.'''
  bdp = ''; bdc = ''
  good_proteins = []
  for p in PROTEINS :
#    check = False; check_c = False
    for c in CONDITION : 
      check = False
      check_c = False
      rep = get_protein_replicates(pd_samples, p, c)
      inter = len(intersect_index(rep))
      if inter != 0:
        check = True
        bdp = rep[0]['Database'][0]
        ctr = 'MG1655 (placI)mVenus-SPA-pUC19 (TYPE C2)'
        rep = get_control_replicates(pd_controls, ctr, c)
        bdc = rep[0]['Database'][0]
        if bdp == bdc:
          good_proteins.append((p,c))
  return good_proteins

#dt.header('aa')
#print(good_proteins())

def accession_list(df):
  '''Returns a list of indexes that is to say a list of accession numbers.'''
  prot = []
  string = ''
  for i in df.index.values:
    string += ','
    string += i
  string = string.replace(',','',1) # remove the first ","
  return string

def test_bioP():
#  bfNumber = 'W1'
#  df = dt.load_df(bfNumber)
#  string = accession_list(df)
  dt.header('bioP')
  Entrez.email = 'maxime.Carlier@insa-lyon.fr'
  handle = Entrez.efetch(db="protein", id="KFI00769.1", rettype="gp", retmode="xml")
#  handle = Entrez.efetch(db='protein', id='WP_074458172.1', rettype="gp", retmode="xml") # Nprot
  #handle = Entrez.efetch(db="protein", id="P0ABD5", rettype="gb", retmode="xml") #Sprot
  record = Entrez.read(handle)
  dt.header('different keys')
  print(record[0].keys()) # [0] to have the first id given. 
  dt.header('source_db')
  print(record[0]['GBSeq_source'])
  dt.header('accession number')
  print(record[0]['GBSeq_other-seqids'])
#  print(record[0]['GBSeq_organism'])
#  print(record[0]['GBSeq_length'])
  dt.header('features')
  print(record[0]['GBSeq_feature-table'])
  dt.header('feature Protein')
  print(record[0]['GBSeq_feature-table'][2])
#  a = []
#  for i in range(len(record)):
#    b =[]
#    for feature in record[i]['GBSeq_feature-table']:
#      if feature['GBFeature_key'] in ['CDS', 'Protein', 'gene']:
#        b.append(feature['GBFeature_key'])
#    print((b))    
#    a.append(b)
  handle.close()

#test_bioP()

def add_columns_df(df, string):
  '''Find the appropriate gene_name for each protein, and also protein name, organism, ID. Return a df that contains these new columns.'''
  access = []
  definition = []
  name = []
  organism = []
  gene = []
  Entrez.email = 'maxime.Carlier@insa-lyon.fr'
  handle = Entrez.efetch(db="protein", id=string, rettype="gb", retmode="xml") # Nprot
  result = Entrez.read(handle)
  handle.close()
  for res in result :
    accession = ""
    for i in res['GBSeq_other-seqids']:
      if 'gi' in i:
        accession = i
    access.append(accession)
    definition.append(res['GBSeq_definition'])
    is_gene = False
    is_prot = False
    is_org = False
    for feature in res['GBSeq_feature-table']:
      if feature['GBFeature_key'] == 'Protein':
        name.append(feature['GBFeature_quals'][0]['GBQualifier_value'])
        is_prot = True
      elif feature['GBFeature_key'] == 'gene':
        if is_gene == False :
          gene.append(feature['GBFeature_quals'][0]['GBQualifier_value'])        
          is_gene = True
      elif feature['GBFeature_key'] == 'source':
        if is_org == False:
          organism.append(feature['GBFeature_quals'][0]['GBQualifier_value'])
          is_org = True
      elif feature['GBFeature_key'] == 'CDS':
        if is_gene == False:
          for j in feature['GBFeature_quals']:
            if is_gene == False:
              if j['GBQualifier_name']== 'gene':
                gene.append(j['GBQualifier_value'])
                is_gene = True
#    if is_gene == False:
#      for feature in res['GBSeq_feature-table']:
#        if feature['GBFeature_key'] == 'Region':
#          if is_gene == False:
#            gene.append(feature['GBFeature_quals'][0]['GBQualifier_value'])
#            is_gene = True
    if is_gene == False : 
      gene.append('None')
    if is_prot == False : 
      name.append('None')
    if is_org == False : 
      organism.append('None')
  print(gene)
  print(len(gene))
  print(len(df.index))
  df['Accession_number'] = access # recupere tous les numéros d'accession
  df['Protein_name'] = name # recupere tous les noms de proteines
  df['Organism'] = organism # recupere tous les numéros d'organisms
  df['Gene_name'] = gene # recupere tous les numéros de gènes
  df = df.reset_index()
  df = df.set_index('Accession_number')
  return df

def create_csv_genes(bfNumber):
  '''Create new csv files with information such as gene_name...'''
# W1, W2 sprot et G2 ncbiprot
# A1, E1 ncbinr I4 ncbiprot
# F2, U6, U7 : ncbiprot
  df = dt.load_df(bfNumber)
  print('add_col')
  df = add_columns_df(df, accession_list(df))
  df = df[~df.Gene_name.str.contains("None", case=False)] # Remove rows without genes.
#  df = df[df.Organism.str.contains("Escherichia coli", case=False)]
#  df = df[df['Num. of significant sequences']>1]
  for i in pd.unique(df['Family']): # Remove redundant rows, keep max sig seq. 
    if len(df.loc[df['Family'] == i]) > 1 :
      df_fam = df.loc[df['Family'] == i]
      if len(df_fam[df_fam.Organism.str.contains("Escherichia coli", case=False)].index) >= 1:
        df_fam = df_fam[df_fam.Organism.str.contains("Escherichia coli", case=False)] # Remove rows wrong organism if at least one from e.coli
      df_fam = df_fam.sort_values(by=['Num. of significant sequences','Num. of significant matches'] , ascending=False)
      indexNames = df_fam.iloc[1:].index
      df.drop(indexNames , inplace=True)

  path_batch = "MS data csv reportbuilder-20200408T212019Z-001/MS data csv reportbuilder/"
  if len(bfNumber) == 1:
      letter = ''
      number = bfNumber
  else :
    letter = bfNumber[0]
    number = bfNumber[1:]
  full_path = glob.glob(path_batch+"batch "+letter+"/*"+letter+number+".reportbuilder.csv")
  df.to_csv(path_batch+"batch "+letter+"/genes_"+letter+number+'.csv')

#create_csv_genes('T10')

def create_csv_genes_good_proteins():
  '''Create all the new files that will be used.'''
  dt.header('proteins')
  for my_tuple in dt.good_proteins():
    batch_names = dt.get_batches(pd_samples,my_tuple[0], my_tuple[1])
    for bname in batch_names:
      print(bname)
#      create_csv_genes(bname)
  dt.header('controls')
  controls_typeC = ['S1','S2','S3', 'O9', 'R4', 'R5', 'R1', 'R2', 'R3']
  for bname in controls_typeC:
    print(bname)
#    create_csv_genes(bname)
  controls_typeA = ['L1', 'T7', 'T8', 'T9', 'C13', 'P5', 'U10', 'A10', 'T5', 'T6']
  for bname in controls_typeA:
    print(bname)
#    create_csv_genes(bname)

def load_df_genes(bfNumber):
  '''Load dataframe from new files with all gene names.'''
  os.chdir("/home/carlier/Documents/Stage/Interactome_proteins_ecoli/") # where is the MS data
  path_batch = "MS data csv reportbuilder-20200408T212019Z-001/MS data csv reportbuilder/"
  letter = bfNumber[0]
  number = bfNumber[1:]
  full_path = path_batch+"batch "+letter+"/genes_"+letter+number+'.csv'
  df = pd.read_csv(full_path, sep=',', header=0)
  df = df.set_index("Accession_number")
  return df

def remove_duplicate_genes(rep):
  '''Some genes may be present in several copies and we only want them once. We keep the first value identified because they are supposed to be identical.'''
  rep["Gene_name"] = rep['Gene_name'].apply(lambda x: x.split('_')[0])
  unique_genes = rep.Gene_name.unique()
  a = rep[rep.Gene_name.isin(unique_genes)]
  if len(a.index) != len(unique_genes):
    duplicate = a.groupby('Gene_name').count()
    indexNames = duplicate[duplicate.Accession == 1].index # get index when only one row is present per gene
    duplicate = duplicate.drop(indexNames)
    for i in duplicate.index:
      df_duplicate = a.loc[a['Gene_name'] == i]
      indexNames = df_duplicate.iloc[1:].index
      a = a.drop(indexNames) # remove when a gene is repeated. 
  return a

#rep = load_df_genes('U3')
#print(len(rep))
#res = remove_duplicate_genes(rep)
#print(len(res))

def create_csv_unique_gene(bfNumber):
  '''Remove all duplicate genes and create new csv files for a given bfNUmber.'''
  df = load_df_genes(bfNumber)
  df = remove_duplicate_genes(df)
  if len(bfNumber) == 1:
      letter = ''
      number = bfNumber
  else :
    letter = bfNumber[0]
    number = bfNumber[1:]
  path_batch = "MS data csv reportbuilder-20200408T212019Z-001/MS data csv reportbuilder/"
  df.to_csv(path_batch+"batch "+letter+"/unique_gene_"+letter+number+'.csv')

def create_all_csv_unique_gene():
  '''Create new csv that contain each gene only once for each bfNumber used.'''
  used_prot_tuples = dt.good_proteins()
  for prot1 in used_prot_tuples:
    prot_batches = dt.get_batches(pd_samples, prot1[0], prot1[1])
    for bname in prot_batches:
      print(bname)
      create_csv_unique_gene(bname)
  controls_typeC = ['S1','S2','S3', 'O9', 'R4', 'R5', 'R1', 'R2', 'R3']
  for bname in controls_typeC:
    create_csv_unique_gene(bname)
  controls_typeA = ['L1', 'T7', 'T8', 'T9', 'C13', 'P5', 'U10', 'A10', 'T5', 'T6']
  for bname in controls_typeA:
    create_csv_unique_gene(bname)

#create_csv_genes_good_proteins()
#create_all_csv_unique_gene()
# be careful, C13 and A10 are in NCBInr database (controls_typeA in LB O/N and M9 0.2% ac O/N)

