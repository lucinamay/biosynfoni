# -*- coding: utf-8 -*-
"""

||||||||||||  ㅇㅅㅇ  |||||||||||||||
____________________________________

title: inoutput                     ||
creaetd: 2023.09.13 14:40           ||
author: lucina-may nollen           || 
institute: WUR Bioinformatics       ||
____________________________________

||||||||||||  ()()()  |||||||||||||||

description:    general functions needed in the other file

======= input =======
- readr(filename:str,ignore_start='',
          encoding='UTF-8',errors='ignore')->list
- entry_parser(lines:list, sep:str = "$$$$")->list[list]
- extractr(lines:list, starts:list=[], remove_starts:bool= True, 
            strip_extra:str=' - ')->dict[list[str]]
======= general =======
- per_entry(entries,fnx,*args,**kwargs)->list
- dictify(listoflists:list[list])->dict
- tuplist_flattenr(list_of_tuples:list[tuple[tuple]])->list[list]
- clean_dict(dic:dict)
======= storage =======
- picklr(cucumber:Any, title:str)->str
- jaropener(picklepath='')-> pickle
======= output =======
- outfile_namer(filename_root:str,
                  addition:str = '',
                  extraroot:bool=True)->str
- output_direr(dirname:str='./output')->tuple[str]
- fp_writr(fp_array:list, outfile_name: str)->None

"""
import pickle
import os
import json
from datetime import date
from typing import Any


# ============================== input handling  ==============================


def readr(filename: str, ignore_start="", encoding="UTF-8", errors="ignore") -> list:
    """
    reads in file in ls of lines (stripped)
    input: (str) filename -- file to be read
           (None / int) len_stop -- stop value for shorter reading
    """
    # dev function:
    # len_stop = None

    # main functionality
    all_lines = []
    with open(filename, "r", encoding=encoding, errors=errors) as fil:
        for line in fil:
            if ignore_start:
                if line.startswith(ignore_start):
                    continue
            all_lines.append(line.strip())
    return all_lines


def entry_parser(lines: list, sep: str = "$$$$") -> list[list]:
    entries = []
    current_entry = []
    for line in lines:
        current_entry.append(line)
        if line == sep:
            entries.append(current_entry)
            current_entry = []  # new entry
    return entries


def extractr(
    lines: list, starts: list = [], remove_starts: bool = True, strip_extra: str = " - "
) -> dict[list[str]]:
    """within a collection of lines, only extracts the lines starting with
    terms in starts, returning them as a list to accomodate multiple values
    per entry"""
    extraction = {}
    for start in starts:
        extraction[start] = []
    for line in lines:
        for start in starts:
            if line.startswith(start):
                extraction[start].append(
                    line.strip().replace(start, "").replace(strip_extra, "").strip()
                )
    return extraction


# ============================= general handling  =============================


def per_entry(entries, fnx, *args, **kwargs) -> list:
    """loops given function over given entries, passes any additional args to
    the function as follows: fnx(entries,*args,**kwargs)"""
    all_vals = []
    for entry in entries:
        val = fnx(entry, *args, **kwargs)
        all_vals.append(val)
    return all_vals


def dictify(listoflists: list[list]) -> dict:
    """for a list of lists, turns it into a dictionary where the first item is
    the key and the second the value. consecutive items are ignored. lists of
    one item are also ignored"""
    assert isinstance(listoflists, list), "input error for dictify, not list"
    dictionary = {}
    for item in listoflists:
        if isinstance(item, list) and len(item) >= 2:
            dictionary[item[0]] = item[1]
    return dictionary


def get_twovals(
    entry: list[str], start1: str, start2: str, start_val_sep=" - "
) -> list[str]:
    val1, val2 = "", ""
    for line in entry:
        if line.startswith(start1):
            val1 = line.split(start_val_sep)[-1]
        elif line.startswith(start2):
            val2 = line.split(start_val_sep)[-1]
    return [val1, val2]


def entryfile_dictify(
    filename: str,
    keyvals: tuple[str],
    start_val_sep: str,
    entry_sep: str = "//",
    encoding: str = "UTF-8",
) -> dict:
    """makes dictionary out of files"""
    annot_entries = entry_parser(readr(filename, encoding=encoding), sep=entry_sep)
    annot = dictify(
        per_entry(
            annot_entries,
            get_twovals,
            start1=keyvals[0],
            start2=keyvals[1],
            start_val_sep=start_val_sep,
        )
    )
    return annot


def tuplist_flattenr(list_of_tuples: list[tuple[tuple]]) -> list[list]:
    per_bit = []
    for bit in list_of_tuples:
        per_bit.append(list(sum(bit, ())))
    return per_bit


def clean_dict(dic: dict):
    """cleans dictionary entries with empty values"""
    clean = {}
    for key, val in dic.items():
        if val:
            clean[key] = val
    return clean


# =========================== temporary storage etc ===========================


def picklr(cucumber: Any, title: str) -> str:
    picklename = "{}.pickle".format(outfile_namer(title))
    if os.path.exists(picklename):
        print("Pickle exists at this date, will write to {}_1")
        picklename = "{}_1".format(picklename)

    with open(picklename, "bw") as pkl:
        pickle.dump(cucumber, pkl, protocol=-1)
    return picklename


def jaropener(picklepath=""):
    if not picklepath:
        picklepath = outfile_namer("wikifile", "pickle")
    with open(picklepath, "rb") as pkl:
        return pickle.load(pkl)


def dump_json(data, filename: str) -> None:
    with open(filename, "w") as f:
        json.dump(data, f)
    return None


def open_json(filename):
    with open(filename, "r") as f:
        data = json.load(f)
    return data


# ================================== output ===================================


def outfile_namer(
    filename_root: str, addition: str = "", extraroot: bool = True
) -> str:
    """gives outfile names with the date prefix 'MMDD_'
    input: (str) filename_root -- main name (e.g. the cdk sdf  filename)
    output: (str) filename without extension
    * requires datetime -- 'from datetime import date'
    """
    # get date prefix --------------------------------------------------
    today = date.today()
    my_date = today.strftime("%m%d")

    # cleaning up root name --------------------------------------------
    ifrem = "{}_".format(my_date)  # for redundant date removal
    if extraroot:
        # in case not rooty enough:
        realroot = filename_root.split("/")[-1].split(".")[0].replace(ifrem, "")
    else:
        realroot = filename_root.replace(ifrem, "")

    # main functionality -----------------------------------------------
    if addition:
        outfile_name = "{}_{}_{}".format(my_date, realroot, addition)
    else:
        outfile_name = "{}_{}".format(my_date, realroot)

    return outfile_name


def output_direr(dirname: str = "./output") -> tuple[str]:
    """moves to the right directory, if dir doesn't exist yet, will make dir
    (default uses relative location!: './output')
    * dependency: os
    """
    init_dir = os.getcwd()

    if not os.path.exists(dirname):
        os.mkdir(dirname)
    os.chdir(dirname)
    return dirname, init_dir


def csv_writr(lst: list[list], outfile: str, sep: str = ",") -> None:
    """writes csv of list of lists"""
    with open(outfile, "w") as of:
        for np in lst:
            if isinstance(np, list) or isinstance(np, tuple):
                of.write(sep.join([str(x) for x in np]))
            elif isinstance(np, str):
                of.write(np)
            elif isinstance(np, int):
                of.write(str(np))
            of.write("\n")
    return None


# ============================== recording ===================================


