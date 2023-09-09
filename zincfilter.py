"""

||||||||||||  ㅇㅅㅇ  |||||||||||||||
____________________________________

title: zincfilter                   |
language: python                    |
license: MIT                        |
author: Lucina-May Nollen           | 
institute: WUR Bioinformatics       |
student no: 1197142                 |
____________________________________

||||||||||||  ()()()  |||||||||||||||

description:    for collection (sdf) of natural products + zinc (NP + non-NP)
                removes all natural products from the collection, to (ideally)
                leave only non-natural products.
                does this through SMILES
                returns:
                - sdf of 'non-NP'

style:          attempting to follow PEP8 styleguide

"""
#--------------------------------- IMPORTS-------------------------------------

from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit import RDLogger  #for muting warnings
from sys import argv
from datetime import date
import os

#================================= CLASSES ====================================
class MyFile:
    def __init__(self, path):
        self.path = path
        self.name = path.split("/")[-1].split(".")[0]

    def readr(self):
        with open(self.path, 'r') as filetoread:
            lines = [line.strip() for line in filetoread]
        return lines

class MolSupplier(MyFile):
    """molecule supplier file (.sdf)
    """
    def __init__(self, path, info_titles: list):
        #further clean up the name-------------
        MyFile.__init__(self, path) #or use super()? --> does not work
        self.name = self.name.strip("_rdk") #does not strip if not exist
        self.supplier = Chem.SDMolSupplier(self.path)

        #extra info of interest: (list) of titles of the extra info in sdf
        self.info_titles = info_titles
    
    def infopath(self) -> str:
        """returns path for file to save additional info to of sdf
        """
        return "{}{}".format(path.strip("_rdk.sdf"),"_info.tsv")
        
    def statspath(self) -> str:
        """returns path for file containing error-producing molfiles"""
        return "{}{}".format(path.strip("_rdk.sdf"),"_stats.csv")
    
    def molslist(self) -> list:
        return [x for x in self.supplier]
    
    def entries(self) -> list:
        lines = self.readr()
        entries, current_entry = [],[] #initialise
        for line in lines:
            current_entry.append(line)
            if line == '$$$$':
                entries.append(current_entry)
                current_entry = []  #new entry
        return entries
    
    def _empty_infodict(self)-> dict:
        """make dictionary of the info-type-descriptors so the order of them
        in the sdf does not matter for parsing and always returns the right
        info in the right 'column' of the info array"""
        info_dict = {}
        for title in self.info_titles():
            info_dict[title]=''
        return info_dict

    def _info_announcers(self) -> list:
        """makes the titles represent sdf format of extra info"""
        return [">  <{}>".format(name) for name in self.info_titles]
    
    def _info_per_entry(self, entry) -> list:
        #empty dictionary to hold un-ordered info
        info_dict = self._empty_infodict()
        info_announcers = self._info_announcers()

        #go through entry looking for the announcers
        for e_i in range(len(entry)): #e_i==line index for entry
            if not entry[e_i].startswith(">"): #prevent unnecces. checking
                continue
            else:
                for a_i in info_announcers:
                    if entry[e_i].startswith(info_announcer[a_i]):
                        info_dict[titles[a_i]] = entry[e_i+1]
        
        #export info in the order of the self.info_titles
        entry_info = [info_dict[title] for title in self.info_titles]
        return entry_info

    def info(self) -> list:
        infos = []
        # {should i add in the self._empty_infodict() and info_announcers here}
        for entry in self.entries():
            infos.append(self._info_per_entry(entry))
        return infos

                    


#================================ FUNCTIONS ===================================
#================================== input =====================================

def readr(filename):
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

def sdf_entry_parser(lines):
    entries = []
    current_entry = []
    for line in lines:
        current_entry.append(line)
        if line == '$$$$':
            entries.append(current_entry)
            current_entry = []  #new entry
    return entries

#================================ converter ====================================

def suppl_to_mol_er(sdf_file):
    return [x for x in sdf_file]

def propertifier(sdf_entries): #22 sec
    """lst of properties in entry as follows:
    every entry is [[zinc_id],[smiles]]
    #could be method in class 'MySupplier'
    """
    entries_properties = []
    for entry in sdf_entries:
        i = 0
        zinc_id,smiles = '',''

        for i in range(len(entry)):
            if entry[i].startswith('>  <zinc_id>'):
                zinc_id = entry[i+1]
            elif entry[i].startswith('> <smiles>'):
                smiles = entry[i+1]
        entries_properties.append([zinc_id,smiles])
    print('{} sdf entries parsed'.format(len(entries_properties))
    return entries_properties


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

def outfile_namer(filename_root: str, addition = '') -> str:
    """gives outfile names in the right formats (with the date 'MMDD_' prefix)
    input: (str) filename_root -- main name (e.g. the cdk sdf  filename)
    requires datetime -- 'from datetime import date' 
    """
    today = date.today()
    my_date = today.strftime("%m%d")
    
    if addition:
        outfile_name = "{}_{}_{}".format(my_date, filename_root, addition)
    else: 
        outfile_name = "{}_{}".format(my_date, filename_root)

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
    #read in all smiles into two collections: np and all
        #read(np database) --> collection
        #if mol in np-database:
            #do not add in the read(all database) --> collection
    #write zinc_id to smiles file for small-data storage
    #use zinc_id selection to select the mol-files (rdkit) to keep

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
