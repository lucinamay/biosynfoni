import sys,os

import pandas as pd
from rdkit import Chem

sys.path.append(os.path.abspath(os.path.join(sys.path[0], os.pardir, "src")))
from biosynfoni.inoutput import readr


init_wd = os.getcwd()

synthetic_zinc_numfile = 'my_synzincnum.tmp'
all_zinc_numfile = 'my_allzincnum.tmp'
all_smilesfile = 'my_allsmiles.tmp' #indexes corr to all zincnum

syn_num = pd.read_csv(synthetic_zinc_numfile)
syn_num.columns = ['synzinc_num']
syn_num.astype(int)

all_zincsmiles = pd.read_csv(all_zinc_numfile,)
all_zincsmiles.columns = ['zinc_id_num']
all_zincsmiles.astype(int)

all_smiles = pd.read_csv(all_smilesfile)
all_smiles.columns = ['smiles']

all_zincsmiles['smiles'] = all_smiles['smiles']

synsmiles = all_zincsmiles[all_zincsmiles['zinc_id_num'].isin(
                            syn_num['synzinc_num'])]

synsmiles['smiles'].to_csv('my_synsmiles.smi',header=False,index=False)
synsmiles.to_csv('my_syns.smi',sep='\t',header=False,index=False)

smiles_list=list(synsmiles['smiles'])
mols_list = []
error_list = []
for s_i in range(len(smiles_list)):
    mol = Chem.MolFromSmiles(smiles_list[s_i])
    if mol:
        mols_list.append(mol)
    else:
        error_list.append(smiles_list[s_i])
        try:
            synsmiles.loc[synsmiles[s_i], 'error'] = 1
        except:
            synsmiles.loc[s_i, 'error'] = 1
            

if error_list:
    errors_df = synsmiles[synsmiles['error']==1]
    errors_df.to_csv('my_molfromsmi_errors.smi',sep='\t',header=False,index=False)

writer = Chem.SDWriter('my_synmols.sdf')
for ind_mol in mols_list:
    writer.write(mol)

