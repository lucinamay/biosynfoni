#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 14 10:42:22 2023

@author: lucina-may
"""
import subprocess

#treemaker
#order is not known, so have to search throughout file 

def getchilds(parent_uid = 805080, tree='../ott3.5/taxonomy.tsv') -> list:
    # try: grep -E '\t805080[^0-9+\-]' taxonomy.tsv
    par_cmd = ['grep','-E','\t{}[^0-9+\-]'.format(parent_uid), tree]  
    caught = subprocess.run(par_cmd, stdout= subprocess.PIPE)
    results= caught.stdout.decode().split('\n')
    
    children=[]
    #output error handling
    if not results: #gives false so empty list
        return children
    else:
        for child in results: 
            if child:
                child_uid = int(child.split('|')[0].strip('\t').strip())
                child_name = child.split('|')[2].strip()
                child_level = child.split('|')[3].strip()
                children.append([child_uid, child_name, child_level])
        return children
   
def make_tree_dict(parent_uid = 805080)->dict:
    tree_dict={}
    tree_dict[parent_uid]={}
    children = getchilds(parent_uid)
    if not children:
        return tree_dict
    else:
        for child in children:
            child_uid = child[0]
            child_name = child[1]
            child_level = child[2]
            tree_dict[parent_uid][child_uid]={}
            tree_dict[parent_uid][child_uid]['name']=child_name
            tree_dict[parent_uid][child_uid]['level']=child_level
            print('recursion', child_name)
            tree_dict[parent_uid][child_uid]['children']=make_tree_dict(
                                                                child_uid)
        return tree_dict
    
tree = make_tree_dict(805080)