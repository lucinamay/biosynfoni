import sys,os

import pandas as pd
from rdkit import Chem

sys.path.append(os.path.abspath(os.path.join(sys.path[0], os.pardir, "src")))
from biosynfoni.inoutput import readr, outfile_namer, csv_writr


init_wd = os.getcwd()

synthetic_zinc_numfile = 'my_synzincnum.tmp' #all synthetic zinc numbers
all_zinc_numfile = 'my_allzincnum.tmp' #all zinc numbers
all_smilesfile = 'my_allsmiles.tmp' #indexes corresp. to all zincnum?????

#final file with synthetic smiles:
outfilename='my_syns.smi'

if not os.path.exists(outfilename):
    #making a df with all synthetic zinc ID's
    syn_num = pd.read_csv(synthetic_zinc_numfile)
    syn_num.columns = ['synzinc_num']
    syn_num.astype(int)

    #making a df with all zinc -- smiles pairs
    all_zincsmiles = pd.read_csv(all_zinc_numfile,)
    all_zincsmiles.columns = ['zinc_id_num']
    all_zincsmiles.astype(int)

    all_smiles = pd.read_csv(all_smilesfile)
    all_smiles.columns = ['smiles']

    all_zincsmiles['smiles'] = all_smiles['smiles']

    #making a df with only synthetic zinc -- smiles pairs
    synsmiles = all_zincsmiles[all_zincsmiles['zinc_id_num'].isin(
                                syn_num['synzinc_num'])]

    synsmiles['smiles'].to_csv('my_synsmiles.smi',header=False,index=False)
    synsmiles.to_csv('my_syns.smi',sep='\t',header=False,index=False)
else:
    synsmiles = pd.read_csv(outfilename,sep='\t',header=None)
    synsmiles.columns = ['zinc_id_num','smiles']
smiles_list=synsmiles['smiles'].tolist()
index_list = synsmiles['zinc_id_num'].tolist()
mols_list = []
mols_ids = []
error_list = []
for s_i in range(len(smiles_list)):
    mol = Chem.MolFromSmiles(smiles_list[s_i])
    if mol:
        mols_list.append(mol)
        mols_ids.append(index_list[s_i])
    else:
        error_list.append([index_list[s_i],smiles_list[s_i]])
            

if error_list:
    rootname = outfile_namer("my_molfromsmi_errors")
    csv_writr(error_list,f"{rootname}.tsv",sep='\t')

rootname = outfile_namer("my_synmols")
writer = Chem.SDWriter(f"{rootname}.sdf")
for ind_mol in mols_list:
    writer.write(ind_mol)

