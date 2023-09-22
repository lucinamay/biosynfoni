#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

||||||||||||  ㅇㅅㅇ  |||||||||||||||
____________________________________

title: inoutput                     ||
creaetd: 2023.09.13 14:40           ||
language: python                    ||
author: lucina-may nollen           || 
institute: WUR Bioinformatics       ||
student no: 1197142                 ||
____________________________________

||||||||||||  ()()()  |||||||||||||||

description:    general functions needed in the other file

- readr(filename:list)
- entry_parser(lines:list, entry_sep:str = "$$$$")

"""
import pickle
import os
from datetime import date
from typing import Any

def readr(filename:str) -> list: 
    """
    reads in file in ls of lines (stripped)
    input: (str) filename -- file to be read
           (None / int) len_stop -- stop value for shorter reading
    """
    #dev function:
    #len_stop = None
    
    #main functionality
    all_lines = []
    with open (filename,'r') as fil:
        for line in fil:
            all_lines.append(line.strip())
    return all_lines

def entry_parser(lines:list, entry_sep:str = "$$$$") -> list[list]:
    entries = []
    current_entry = []
    for line in lines:
        current_entry.append(line)
        if line == entry_sep:
            entries.append(current_entry)
            current_entry = []  #new entry
    return entries

def picklr(cucumber:Any, title:str)-> str:
    picklename = "{}.pickle".format(title)
    if os.path.exists(picklename):
        print("Pickle exists at this date, will write to {}_1")
        picklename = "{}_1".format(picklename)
        
    with open(picklename, 'bw') as pkl:
        pickle.dump(cucumber, pkl, protocol=-1)
    return picklename

def outfile_namer(filename_root:str,
                  addition:str = '',
                  extraroot:bool=True) -> str:
    """gives outfile names with the date prefix 'MMDD_' 
    input: (str) filename_root -- main name (e.g. the cdk sdf  filename)
    output: (str) filename without extension
    * requires datetime -- 'from datetime import date'
    """
    #get date prefix --------------------------------------------------
    today = date.today()
    my_date = today.strftime("%m%d")
    
    #cleaning up root name --------------------------------------------
    ifrem = "{}_".format(my_date)       #for redundant date removal
    if extraroot:
        #in case not rooty enough:
        realroot = filename_root.split('/')[-1].split('.')[0].strip(ifrem)
    else:
        realroot = filename_root.strip(ifrem)

    #main functionality -----------------------------------------------
    if addition:
        outfile_name = "{}_{}_{}".format(my_date, realroot, addition)
    else:
        outfile_name = "{}_{}".format(my_date, realroot)

    return outfile_name

def output_direr(dirname:str='./output') -> tuple[str]:
    """moves to the right directory, if dir doesn't exist yet, will make dir
    (default uses relative location!: './output')
    * dependency: os
    """
    init_dir = os.getcwd()
    
    if not os.path.exists(dirname):
        os.mkdir(dirname)
    os.chdir(dirname)
    return dirname, init_dir

def flatten_tuplelist(list_of_tuples:list[tuple[tuple]])->list[list]:
    per_bit = []
    for bit in list_of_tuples:
        per_bit.append(list(sum(bit, ())))
    return per_bit
    
def fp_writr(fp_array:list, outfile_name: str)->None:
    """writes csv of fingerprint
    """
    with open("{}.csv".format(outfile_name),'w') as biosyn:
        for row in fp_array:
            biosyn.write(','.join(str(col) for col in row))
            biosyn.write('\n')
    return None



