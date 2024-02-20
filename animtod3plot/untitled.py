# -*- coding: utf-8 -*-
"""
Created on Sun Feb 18 20:29:20 2024

@author: PC
"""

from lasso.dyna.d3plot import D3plot
import numpy as np


d3 = [D3plot("C:/Users/PC/Downloads/test2/d3plot"), D3plot("C:/Users/PC/Downloads/test4/CRA2AV4.d3plot")]

hmm = ['part_ids', 'part_ids_unordered', 'part_ids_cross_references', 'part_titles_ids', 'part_titles']



    
for ii in list(d3[0].arrays.keys()):
    for i in d3:
        
        try:
            print(ii)
            print(i.arrays[ii])
            print(i.arrays[ii])
            print(i.arrays[ii].dtype)
        except:
            print("no found: " + str(ii))

        
    ["release_version","runtime","soleng","source_version", "title", "tsheng", "version"]
