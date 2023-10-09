"""
||||||||||||  ㅇㅅㅇ  |||||||||||||||
____________________________________

title: wiki_combine                 |
author: Lucina-May Nollen           | 
institute: WUR Bioinformatics       |
____________________________________

||||||||||||  ()()()  |||||||||||||||

[1] wikidata json [2] filename taxonomy
"""
# %%
from sys import argv
from datetime import date
import os
import json
import subprocess
import pandas as pd

from inoutput import picklr, jaropener, outfile_namer

# init.
global TAX_TREE
TAX_TREE = "taxonomy.tsv"  # change to path of choice in main


def dictr(jsonpath):
    with open(jsonpath, 'r') as jsonfile:
        dicted = json.load(jsonfile)
    return dicted


def dfr(jsonpath):
    with open(jsonpath, 'r'):
        dfed = pd.read_json(jsonpath, orient='records')
    return dfed


class MyTaxon:
    def __init__(self, name, ncbi_id, tax_tree='',
                 taxinfo_header=('domain',
                                 'kingdom',
                                 'phylum',
                                 'subphylum',
                                 'order')):
        if not tax_tree:
            global TAX_TREE
            tax_tree = TAX_TREE
        self.ncbi_id = ncbi_id
        self.name = name
        self.tax_tree = tax_tree  # make sure absolute path
        self.taxinfo_header = taxinfo_header

    def __checkpath__(self) -> bool:
        if not os.path.exists(self.tax_tree):
            raise NameError('Taxonomical tree path {} does not exist'.format(
                            self.tax_tree))
        return True

    def __ncbi_to_id__(self) -> int:
        """ using grep for speed and memory saving
        """
        # error handling
        self.__checkpath__()

        # main functionality
        ncbi_cmd = ['grep',
                    'ncbi:{}[^0-9+\-]'.format(self.ncbi_id),
                    self.tax_tree]
        caught = subprocess.run(ncbi_cmd, stdout=subprocess.PIPE)
        taxon_line = caught.stdout.decode()
        if not taxon_line:
            return -1
        else:
            return int(taxon_line.split('|')[0].strip())

    def getparent(self, child_id: int = 0) -> int:
        # input error handling
        if child_id == -1:
            return ('', '', '')

        # main functionality
        if not child_id:
            child_id = self.__ncbi_to_id__()
        par_cmd = ['grep',
                   '^{}[^0-9+\-]'.format(child_id),
                   self.tax_tree]
        caught = subprocess.run(par_cmd, stdout=subprocess.PIPE)
        taxon_line = caught.stdout.decode()

        # output error handling
        if not taxon_line:
            return ('', '', '')
        else:
            parent_id = int(taxon_line.split('|')[1].strip())
            parent_name = taxon_line.split('|')[2].strip()
            parent_level = taxon_line.split('|')[3].strip()
            return parent_id, parent_name, parent_level

    def parent_listr(self) -> tuple:

        # it would be nice to 'combine' the searches within the different
        # compound trees i.e. merging searches after finding common parent to
        # decrease computational load.

        entire_tree = []
        # first parent
        parent = self.getparent()
        if not parent[2]:
            return ()
        else:
            e_i = 0
            while parent[2] != 'domain' and e_i < 1000:
                parent = self.getparent(child_id=parent[0])
                entire_tree.append(parent)
                e_i += 1
            return tuple(entire_tree)

    def taxlvl_getr(self, tree=(), level_oi: str = 'kingdom') -> str:
        if not tree:
            tree = self.parent_listr()
        tax_name = ''  # init.
        for i in range(len(tree)):
            if tree[-i-1][2] == level_oi:
                # gets the highest
                tax_name = tree[1-i][1]
                break
        return tax_name

    def tax_info(self) -> tuple:
        ind_tree = self.parent_listr()
        all_info = []
        header = self.taxinfo_header
        for lvl in header:
            all_info.append(self.taxlvl_getr(ind_tree, get=lvl))
        return tuple(all_info)


def loopr(df):
    for df_i in range(len(df)):
        # for each compound
        one = MyTaxon(
            name=df['taxon_name'][df_i],
            ncbi_id=df['nbic_taxon_id'][df_i],
        )
        all_info = one.tax_info()
        info_header = one.taxinfo_header

        for lvl in range(len(info_header)):
            #write in value
            df.loc[df_i:df_i, (info_header[lvl])] = all_info[lvl]
    return df


def dir_loopr(wikifiles_directory):
    print("started dir loopr")
    initialwd = os.getcwd()
    os.chdir('{}'.format(wikifiles_directory))
    file_list = os.listdir('.')

    with open("loopr.log", 'w') as log:
        log.write("\t".join(file_list))  # to know order just in case

    #loop and write
    tsvname = outfile_namer(wikifiles_directory, 'all')+'.tsv'
    first = True
    for fil in file_list:
        # selecting only df's
        if fil.endswith('.json') or fil.endswith('.pickle'):
            ind_df = perfile(fil)
            if first == True:
                print("newtsv with", fil)
                ind_df.to_csv(tsvname, sep="\t", header=True)
            else:
                print("added to", tsvname, "with", fil)
                ind_df.to_csv(tsvname, sep="\t",
                              mode='a', index=False, header=False)
            first = False

    os.chdir(initialwd)  # set back to original wd
    print("finished loopr")
    return None


def perfile(wikidatafile):
    """returns individual df's
    """
    # check if not looped
    picklename = outfile_namer(wikidatafile, 'taxonomy')
    picklepath = "{}.pickle".format(picklename)
    if os.path.exists(picklepath):
        jaropener(picklepath)

    # using the wikidata jsons
    df = dfr(wikidatafile)

    # looping over all molecules in df
    end_df = loopr(df)

    # output
    picklr(end_df, picklename)
    return end_df


def main():
    # using the wikidata jsons
    jsonsdirectory = argv[1]

    # defining taxonomy file
    taxonomyfile = argv[2]
    global TAX_TREE
    TAX_TREE = taxonomyfile  # move to main later

    dir_loopr(jsonsdirectory)


if __name__ == "__main__":
    main()
