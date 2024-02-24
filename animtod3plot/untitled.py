# -*- coding: utf-8 -*-
"""
Created on Sun Feb 18 20:29:20 2024

@author: PC
"""

from lasso.dyna.d3plot import D3plot
import numpy as np


D3plot("C:/Users/PC/Downloads/test2/reduce.d3plot")



"""
d3 = [D3plot("C:/Users/PC/Downloads/test2/reduce.d3plot"), D3plot("C:/Users/PC/Downloads/test4/CRA2AV4.d3plot")]

hmm = ['part_ids', 'part_ids_unordered', 'part_ids_cross_references', 'part_titles_ids', 'part_titles']


"""
    
#print(d3[0].arrays.keys())


d3[0].arrays["part_ids"] = d3[1].arrays["part_ids"][:-1]


j = list(d3[0].arrays.keys())

for i in j:
    if i not in d3[1].arrays:
        print(i)
        del d3[0].arrays[i]
    
    
    
# Replacing all the part definitions - works    
d3[0].arrays["part_ids"]=d3[1].arrays["part_ids"][:-1]
d3[0].arrays["part_ids_unordered"]=d3[1].arrays["part_ids_unordered"][:-1]
d3[0].arrays["part_ids_cross_references"]=d3[1].arrays["part_ids_cross_references"][:-1]
d3[0].arrays["part_titles_ids"]=d3[1].arrays["part_titles_ids"][:-1]
d3[0].arrays["part_titles"] = d3[1].arrays["part_titles"][:-1]    
    

d3[0].arrays["part_internal_energy"]=d3[0].arrays["part_internal_energy"]
d3[0].arrays["part_kinetic_energy"]=d3[0].arrays["part_kinetic_energy"]
d3[0].arrays["part_velocity"]=d3[0].arrays["part_velocity"]
d3[0].arrays["part_mass"]=d3[0].arrays["part_mass"]
d3[0].arrays["part_hourglass_energy"]=d3[0].arrays["part_hourglass_energy"]

d3[0].arrays["node_displacement"]=d3[0].arrays["node_displacement"]
d3[0].arrays["node_velocity"]=d3[0].arrays["node_velocity"]
d3[0].arrays["node_acceleration"]=d3[0].arrays["node_acceleration"]




for i in d3[1].arrays:
    
    
    
    
    if "shell" in i:

        d3[0].arrays[i]= d3[1].arrays[i]
        
        
for i in d3[1].arrays:
    
    
    if "solid" in i:

        del d3[0].arrays[i]      
        
        
    #if "global" in i:

     #   d3[0].arrays[i]= d3[0].arrays[i][:3]  




d3[0].write_d3plot("C:/Users/PC/Downloads/test4/Fake_CRA2AV4.d3plot")"""


d3_new = D3plot()

d3_new.arrays = d3[0].arrays


d3_new.header.title = d3[0].header.title
d3_new.header.runtime = d3[0].header.runtime
d3_new.header.source_version = d3[0].header.source_version
d3_new.header.release_version = d3[0].header.release_version
d3_new.header.version = d3[0].header.version


d3_new.write_d3plot("C:/Users/PC/Downloads/test2/reduce_created.d3plot")

d3_new_2 = D3plot("C:/Users/PC/Downloads/test2/reduce_created.d3plot")

#print(d3[0].header.raw_header)

#print(d3_new_2.header.raw_header)



for i in d3[0].header.raw_header:
    
    
    if d3[0].header.raw_header[i]  != d3_new_2.header.raw_header[i]:
        
        print(i)
        print(d3[0].header.raw_header[i])
        print(d3_new_2.header.raw_header[i])
        
        

for i in d3[0].arrays:

    
    if not np.array_equal(d3[0].arrays[i].dtype, d3_new_2.arrays[i].dtype):

        print(i)
        print(d3[0].arrays[i])
        print(d3_new_2.arrays[i])



"""




    