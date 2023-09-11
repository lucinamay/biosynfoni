"""

||||||||||||  ㅇㅅㅇ  |||||||||||||||
____________________________________

title: cOnVERTure (cdk_rdk_parser)  |
language: python                    |
author: Lucina-May Nollen           | 
institute: WUR Bioinformatics       |
student no: 1197142                 |
____________________________________

||||||||||||  ()()()  |||||||||||||||

description:    parses (COCONUT) CDK-style sdf's for its CNP-ID, SMILEs,
                InChi and molecular formula, returning 
                [0] a sdf of mol objects
                [1] a tsv with related info for further use in fingerprinter 
                    program
                    [0] CNP-ID
                    [1] InChi
                    [2] SMILES
                    [3] molecular formula
                    [4] NPLS
                [2] a csv with stats on errors
                    [0] (str) -- CNP-ID of entry causing InChi2mol errors
                    [1] (str) -- CNP-ID of entry causing SMILES2mol errors
                    [2] (list) -- CNP-ID of entries causing errors in both
                    [3] (int)  -- total number of entries

style:          attempting to follow PEP8 styleguide

"""
#--------------------------------- IMPORTS-------------------------------------

from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit import RDLogger  #for muting warnings
from sys import argv
from datetime import date
import os

#================================ FUNCTIONS ===================================
#================================== input =====================================

def readr(filename): #takes 9 sec for full file(24764158 line)
    """
    reads in file in ls of lines (stripped)
    input: (str) filename -- file to be read
           (None / int) len_stop -- stop value for shorter reading
    """
    i = 0
    #dev function:
    #len_stop = None
    all_lines = []
    with open (filename,'r') as sdf:
        for line in sdf:
            all_lines.append(line.strip())
    return all_lines

def entry_parser(lines): #takes 19 secs for 59k entries
    entries = []
    current_entry = []
    for line in lines:
        current_entry.append(line)
        if line == '$$$$':
            entries.append(current_entry)
            current_entry = []  #new entry
    return entries

#================================ converter ====================================

def propertifier(sdf_entries): #22 sec
    """lst of properties in entry as follows:
    every entry is [[coconut_id],[inchi],[smiles],[molecular_formula],[NPLS]]
    #for later:     [molecular weight], [name],]
    """
    entries_properties = []
    for entry in sdf_entries:
        #print(entry)
        i = 0
        coconut_id, inchi,smiles, molecular_formula = '','','',''
        NPLS = ''
        for i in range(len(entry)):
            if entry[i] == '> <coconut_id>':
                coconut_id = entry[i+1]
            elif entry[i] == '> <inchi>':
                inchi = entry[i+1]
            elif entry[i]== '> <SMILES>':
                smiles = entry[i+1]
            elif entry[i]== '> <molecular_formula>':
                molecular_formula = entry[i+1]
            elif entry[i] == '> <NPL_score>':
                NPLS = entry[i+1]
        entries_properties.append([coconut_id, inchi,\
                                   smiles, molecular_formula, NPLS])

    return entries_properties

def molifier(NP_property_list, representation_info = False, stats_out = True):
    """ makes mol class out of the NPs from database, using InChi formula,
    if InChi formula is not correct, it tries the SMILES. Keeps track of errors
    in InChi and SMILES, by adding the error-causing entries to the respective
    error-collection
    input:  (list) NP_property_list -- list of entries where
                    [0] CNP-ID
                    [1] InChi
                    [2] SMILES
                    [3] molecular formula
                    [4] NPLS (NaPLeS)
            (bool) stats_out = True -- yes or no return stats

    returns:    (list) all_mols -- mol classes of all successful inchi/smiles
                (list) prop_list -- 
                    [0] CNP-ID
                    [1] InChi
                    [2] SMILES
                    [3] molecular formula
                    [4] NPLS (NaPLeS)
                (list) stats -- if set to true, returns
                    [0] (str) -- entry CNP causing InChi2mol errors
                    [1] (str) -- entry CNP causing SMILES2mol errors
                    [2] (list) -- entries causing errors in both
                    [3] (int)  -- total number of entries
    """
    inchiNone = []
    smilesNone = []
    bothNone = []
    all_mols = []
    prop_list = []

    for entry in NP_property_list:
        smiles_mol = Chem.MolFromSmiles(entry[2])
        inchi_mol = Chem.MolFromInchi(entry[1])
        if inchi_mol is not None:
            all_mols.append(inchi_mol)
            if not representation_info:
                prop_list.append([entry[0]]+entry[3:]) 
                #avoiding faulty represent.
            elif representation_info:
                prop_list.append(entry)
        elif smiles_mol is not None:
            all_mols.append(smiles_mol)
            if not representation_info:
                prop_list.append([entry[0]]+entry[3:]) 
                #avoiding faulty represent.
            elif representation_info:
                prop_list.append(entry)
        else: #these do not end up in the mol collection
            print('{} not successfully molified'.format(entry[0]))
        #stats
        if stats_out:
            if inchi_mol is None:
                inchiNone.append(entry[0])
            if smiles_mol is None:
                smilesNone.append(entry[0])
            if (smiles_mol is None) and (inchi_mol is None):
                bothNone.append(entry[0])
    print('done with molifying  all entries')
    #output (yes or no stats)
    if stats_out:
        stats = [inchiNone, smilesNone, bothNone, len(NP_property_list)]
        return all_mols, prop_list, stats 
    else: 
        return all_mols, prop_list 

#================================= output =====================================

def sdf_writr(mols, outfile):
    """writes sdf of mols"""
    writer = Chem.SDWriter(outfile)
    for mol in mols:
        writer.write(mol)

def csv_writr(lst, outfile, sep = ","):
    """writes csv of list of lists
    """
    with open(outfile, 'w') as of:
        for np in lst:
            if isinstance(np, list):
                of.write(sep.join(np))
            elif isinstance(np,str):
                of.write(np)
            of.write('\n')
    return None

def outfile_namer(filename_root: str, addition = '', extraroot=True) -> str:
    """gives outfile names in the right formats (with the date 'MMDD_' prefix)
    input: (str) filename_root -- main name (e.g. the cdk sdf  filename)
    requires datetime -- 'from datetime import date' 
    """
    today = date.today()
    my_date = today.strftime("%m%d")
    
    if extraroot:
        #in case filename_root is not rooty enough:
        realroot = filename_root.split('/')[-1].split('.')[:-1]

    if addition:
        outfile_name = "{}_{}_{}".format(my_date, realroot, addition)
    else: 
        outfile_name = "{}_{}".format(my_date, real_root)

    return outfile_name

def output_direr(dirname='./output'):
    """moves to the right directory, if dir doesn't exist yet, will make dir
    (default uses relative location!: './output')
    """
    if not os.path.exists(dirname):
        os.mkdir(dirname)
    os.chdir(dirname)
    return None
#++++++++++++++++++++++++++++++++++ main ++++++++++++++++++++++++++++++++++++++

def main():
    #input
    cdk_sdf = argv[1]

    print('\n', 10*'=','\n', 'cOnVERTURE', '\n',10*'=','\n')
    #handling
    lines = readr(cdk_sdf)
    sdf_entries = entry_parser(lines)
    properties = propertifier(sdf_entries)
    all_mols, ids_props, stats = molifier(
        properties,representation_info = True) 
    
    #output
    actual_name = cdk_sdf.strip('.sdf').split('/')[-1]
    outfile_name = outfile_namer(actual_name) #str

    out_rdk_sdf_name = ''.join(outfile_name + '_rdk.sdf')
    out_info_name = ''.join(outfile_name + '_info.tsv')
    out_stats_name = ''.join(outfile_name + '_stats.csv')
    
    output_direr('./output')
    
    #sdf_writr(all_mols, out_rdk_sdf_name)
    csv_writr(ids_props, out_info_name, sep= "\t")
    csv_writr(stats, out_stats_name)
    
    return None

if __name__ == "__main__":
    main()
