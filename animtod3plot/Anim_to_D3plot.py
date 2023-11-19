# -*- coding: utf-8 -*-
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

from Vortex_Radioss.lasso.dyna.d3plot import D3plot
from Vortex_Radioss.lasso.dyna.array_type import ArrayType
from Vortex_Radioss.animtod3plot.RadiossTHReader import RadiossTHReader
from Vortex_Radioss.animtod3plot.RadiossReader import RadiossReader

import numpy as np
import os
import time
import math
from tqdm import tqdm
import glob


class readAndConvert:
    
    def __init__(
        self,
        filepath: str = None,
    ):
        """Constructor for a readAndConvert

        Parameters
        ----------
        filepath: str
            path to a animation files

        Examples
        --------
            >>> from animtod3plot import Anim_to_D3plot
            >>> # read animation files and build d3plot file
            >>> a2d = Anim_to_D3plot.readAndConvert("path_to_animation_files")"""

        self._start=time.time()
        self._d3plot = D3plot()
        self.A_2_D(filepath)



    def sequential(input_array):
        "Radioss puts a load of zeros at the end, increase by 1 to keep Primer happy for sequential node numbering"    
  
        zero_num = max(input_array).astype('int64')
        output = []
        for i in input_array:
            if i != 0:
                output.append(i)
            else:
                zero_num +=1
                output.append(zero_num)
        output = np.array(output)
        return output 

    def generate_sorter(input_array):
       __index_tracker = [i for i in range(0,len(input_array))]
       _index_tracker = [y for x, y in sorted(zip(list(input_array), list(__index_tracker)))]      
       index_tracker = [0] * len(_index_tracker)
       
       for i_a, a in enumerate(_index_tracker):
           index_tracker[a] = i_a    
           
       return np.array(_index_tracker)
   
    def apply_sorter(input_array, tracker_array):
        return input_array[tracker_array]   
    
    def generate_random_name(length, ifile):
        import random
        import string
        # choose from all lowercase letter
        letters = string.ascii_lowercase
        result_str = ''.join(random.choice(letters) for i in range(length))
        return "_" + result_str + str(ifile)    

    def A_2_D(self,file_stem):
    
        if os.path.isfile(file_stem + "d3plot"):
            return True

        state_index = 1
        file_list = glob.glob(file_stem + "A*[0-9]")
        
        shell_output = {}
        shell_output["element_shell_thickness"] = "element_shell_thickness"
        shell_output["element_shell_specific_energy"] = "element_shell_internal_energy"

        original_node_coordinates       = None
        
               
        print("Re-sorting the indexes")
        rr = RadiossReader(file_list[0]) 
        if rr.raw_header["nbFacets"] > 0:
            shell_ids_tracker = readAndConvert.generate_sorter(readAndConvert.sequential(rr.arrays["element_shell_ids"]))

        if rr.raw_header["nbElts1D"] > 0:
            beam_ids_tracker = readAndConvert.generate_sorter(readAndConvert.sequential(rr.arrays["element_beam_ids"]))

        if rr.raw_header["nbElts3D"] > 0:
            solid_ids_tracker = readAndConvert.generate_sorter(readAndConvert.sequential(rr.arrays["element_solid_ids"]))
            
        original_node_coordinates = rr.arrays["node_coordinates"]       
        
        node_ids_out = readAndConvert.sequential(rr.arrays["node_ids"])
        element_shell_ids_out = readAndConvert.apply_sorter(readAndConvert.sequential(rr.arrays["element_shell_ids"]), shell_ids_tracker)
        element_beam_ids_out = readAndConvert.apply_sorter(readAndConvert.sequential(rr.arrays["element_beam_ids"]), beam_ids_tracker)
        element_solid_ids_out = readAndConvert.apply_sorter(readAndConvert.sequential(rr.arrays["element_solid_ids"]), solid_ids_tracker)
        
        del rr
        
        for ifile, file in enumerate(tqdm(file_list)):
            
            
            rr = RadiossReader(file) 
            
            node_acceleration               = []
            node_velocity                   = []
            node_displacement               = []
            
            
            element_beam_is_alive           = []
            element_shell_is_alive          = []
            element_solid_is_alive          = []
            
            element_shell_thickness         = []
            element_shell_specific_energy   = []
            element_shell_plastic_strain    = []
            element_shell_stress            = []
            
            
            part_names                      = []
            
            timesteps                       = []            
        
            "Node Updates"

            if "node_acceleration" in rr.arrays:  
                node_acceleration.append(rr.arrays["node_acceleration"])

            if "node_velocity" in rr.arrays:  
                node_velocity.append(rr.arrays["node_velocity"])
               
            node_displacement.append(rr.arrays["node_coordinates"])
            

            "Beam Updates"    
            
            if rr.raw_header["nbElts1D"] > 0 and "element_beam_is_alive" in rr.arrays:
                element_beam_is_alive.append(readAndConvert.apply_sorter(rr.arrays["element_beam_is_alive"], beam_ids_tracker).astype("<f"))
            
            
            "shell updates"

            if rr.raw_header["nbFacets"] > 0 and "element_shell_is_alive" in rr.arrays:
                element_shell_is_alive.append(readAndConvert.apply_sorter(rr.arrays["element_shell_is_alive"], shell_ids_tracker).astype("<f"))

            element_shell_output = [] 
            if "nbEFunc" in rr.raw_header:        
                for iefun in range(0, rr.raw_header["nbEFunc"]):
                    element_shell_output.append( "element_shell_" + str(rr.raw_arrays["fTextA"][iefun + rr.raw_header["nbFunc"]]).lower().replace(" ", "_").strip())

            #print(element_shell_output)
            
            if "element_shell_thickness" in element_shell_output:
                element_shell_thickness.append(readAndConvert.apply_sorter(rr.arrays["element_shell_thickness"], shell_ids_tracker))

            if "element_shell_specific_energy" in element_shell_output:
                element_shell_specific_energy.append(readAndConvert.apply_sorter(rr.arrays["element_shell_specific_energy"], shell_ids_tracker))

            if "element_shell_plastic_strain" in element_shell_output:
                element_shell_plastic_strain.append(readAndConvert.apply_sorter(rr.arrays["element_shell_plastic_strain"], shell_ids_tracker))

            element_shell_tens_output = [] 
            if "nbTens" in rr.raw_header:        
                for iefun in range(0, rr.raw_header["nbTens"]):
                    element_shell_tens_output.append( "element_shell_" + str(rr.raw_arrays["tTextA"][iefun]).lower().replace(" ", "_").strip())
                    
            #print(element_shell_tens_output)

            if "Stress (upper)" in rr.arrays or "Stress (lower)" in rr.arrays :   
                _element_shell_stress = np.zeros(shape=(rr.raw_header["nbFacets"], 2, 6))
                    
                if "Stress (upper)" in rr.arrays:
                    _element_shell_stress[:, 0, 0] = rr.arrays["Stress (upper)"][:, 0]
                    _element_shell_stress[:, 0, 1] = rr.arrays["Stress (upper)"][:, 1]
                    _element_shell_stress[:, 0, 2] = np.zeros(shape=(rr.raw_header["nbFacets"]))
                    _element_shell_stress[:, 0, 3] = rr.arrays["Stress (upper)"][:, 2]
                    _element_shell_stress[:, 0, 5] = np.zeros(shape=(rr.raw_header["nbFacets"]))    
            
                if "Stress (lower)" in rr.arrays:
                    _element_shell_stress[:, 1, 0] = rr.arrays["Stress (lower)"][:, 0]
                    _element_shell_stress[:, 1, 1] = rr.arrays["Stress (lower)"][:, 1]
                    _element_shell_stress[:, 1, 2] = np.zeros(shape=(rr.raw_header["nbFacets"]))
                    _element_shell_stress[:, 1, 3] = rr.arrays["Stress (lower)"][:, 2]
                    _element_shell_stress[:, 1, 5] = np.zeros(shape=(rr.raw_header["nbFacets"]))  
                    
                if "Stress (upper)" in rr.arrays or "Stress (lower)" in rr.arrays :   
                    element_shell_stress.append(readAndConvert.apply_sorter(_element_shell_stress, shell_ids_tracker))
        

            "Solid Updates"
            
            if rr.raw_header["nbElts3D"] > 0 and "element_solid_is_alive" in rr.arrays:
                element_solid_is_alive.append(readAndConvert.apply_sorter(rr.arrays["element_solid_is_alive"], solid_ids_tracker).astype("<f"))

            

            "Part Updates"
            
            #"Extract and map part mass from TH file"

            #if part_names == []:
            #    part_mass = []
            #    for part_name in rr.raw_arrays["pTextA"]:
            #        part_name = part_name[10:].strip()
            #        part_names.append(part_name)
                    #if part_name in rr_th.array["part"]:
                    #    part_mass.append(rr_th.array["part"][part_name]["MASS"][0][0])
                    #else:
                    #    part_mass.append(0)
            #        part_mass.append(0)
                        
                #if "element_shell_thickness" in rr.arrays:        
                #    part_shell_volume_tmp=[[] for item in range(0,len(rr.raw_arrays["pTextA"]))]
                #    
                #    for i in range (0, rr.raw_header["nbFacets"]):
                #        part_shell_volume_tmp[rr.arrays["element_shell_part_indexes"][i]].append(element_shell_area[i] * rr.arrays["element_shell_thickness"][i])
                #    part_volume=[sum(item) for item in part_shell_volume_tmp]                
                #    
                #    part_density = [a/b if b!= 0 else 0 for a,b in zip(part_mass, part_volume)]
                
            #if "element_shell_thickness" in rr.arrays:
            #    element_shell_specific_energy.append(readAndConvert.apply_sorter([rr.arrays["element_shell_specific_energy"][i_shell]/part_density[rr.arrays["element_shell_part_indexes"][i_shell]] \
            #                                          if part_density[rr.arrays["element_shell_part_indexes"][i_shell]] != 0 else 0 \
            #                                              for i_shell in range(0, rr.raw_header["nbFacets"])], shell_ids_tracker))
                        

            "Timestep Updates"
            timesteps.append(rr.arrays["timesteps"])
        
            "Assign the arrays to the D3PLOT class for writing"
            
            "Nodes"
            self._d3plot.arrays[ArrayType.node_displacement]              = np.array(node_displacement)
            self._d3plot.arrays[ArrayType.node_coordinates]               = original_node_coordinates                        
            
            
            self._d3plot.arrays[ArrayType.node_ids]                       = node_ids_out
    
            if "node_acceleration" in rr.arrays: 
                self._d3plot.arrays[ArrayType.node_acceleration]              = np.array(node_acceleration)
            if "node_velocity" in rr.arrays: 
                self._d3plot.arrays[ArrayType.node_velocity]                  = np.array(node_velocity)
            
            "Shells"    
            if rr.raw_header["nbFacets"] > 0:
                self._d3plot.arrays[ArrayType.element_shell_ids]              = element_shell_ids_out
                self._d3plot.arrays[ArrayType.element_shell_node_indexes]     = readAndConvert.apply_sorter(rr.arrays["element_shell_node_indexes"], shell_ids_tracker)
                self._d3plot.arrays[ArrayType.element_shell_part_indexes]     = readAndConvert.apply_sorter(rr.arrays["element_shell_part_indexes"], shell_ids_tracker)
                self._d3plot.arrays[ArrayType.element_shell_is_alive]         = np.array(element_shell_is_alive)
                
            v=np.array(rr.raw_arrays["pTextA"])
            
            shell_part_ids                                          =  np.array(rr.raw_arrays["pTextA"]).astype("U9").astype(int)
            shell_part_titles                                       =  np.array(rr.raw_arrays["pTextA"]) 
            
            shell_num = rr.raw_header["nbFacets"]
            shell_part_num = len(shell_part_ids)
            
            if element_shell_thickness != []:
                self._d3plot.arrays[ArrayType.element_shell_thickness]        = np.array(element_shell_thickness)
    
            if element_shell_specific_energy != []:
                self._d3plot.arrays[ArrayType.element_shell_internal_energy]  = np.array(element_shell_specific_energy)
    
            #if element_shell_plastic_strain != []:
            #    print("element_shell_plastic_strain != []")
    
                #_element_shell_plastic_strain = np.zeros(shape=(rr.raw_header["nbFacets"], 2))
                #_element_shell_plastic_strain[0 , :] = np.zeros(shape=rr.raw_header["nbFacets"])
                #_element_shell_plastic_strain[:,1] = rr.arrays["element_shell_plastic_strain"]
                #element_shell_plastic_strain.append(readAndConvert.apply_sorter(_element_shell_plastic_strain, shell_ids_tracker))
    
                #self._d3plot.arrays[ArrayType.element_shell_effective_plastic_strain]           = np.array(element_shell_plastic_strain)
    
                
                
                #self._d3plot.arrays[ArrayType.element_shell_effective_plastic_strain,0]  = np.zeros(shape=(rr.raw_header["nbFacets"]))
                #self._d3plot.arrays[ArrayType.element_shell_effective_plastic_strain,1]  = np.array(element_shell_plastic_strain)
    
                #print(element_shell_plastic_strain)
            
                #self._d3plot.arrays[ArrayType.element_shell_effective_plastic_strain] = np.zeros(shape=(states, rr.raw_header["nbFacets"], 2))
            
            
            if element_shell_stress != []:   
                self._d3plot.arrays[ArrayType.element_shell_stress]           = np.array(element_shell_stress)
            
            "Beams"
            
            if rr.raw_header["nbElts1D"] > 0: 
                additional_beam_number                              = rr.raw_header["nbElts1D"] - len(rr.arrays["element_beam_part_indexes"])
                additional_beams                                    = np.zeros(shape=(additional_beam_number))
                additional_beams.fill(max(rr.arrays["element_beam_part_indexes"])+1)
                element_beam_part_indexes                           = np.concatenate([rr.arrays["element_beam_part_indexes"], additional_beams])
            
                #self._d3plot.arrays[ArrayType.element_beam_ids] = rr.arrays["element_beam_ids"]
                self._d3plot.arrays[ArrayType.element_beam_ids]              = element_beam_ids_out
            
                element_beam_node_indexes=np.zeros(shape=(rr.raw_header["nbElts1D"],5))
            
                element_beam_node_indexes[:,0]                      = rr.arrays["element_beam_node_indexes"][:,0]
                element_beam_node_indexes[:,1]                      = rr.arrays["element_beam_node_indexes"][:,1]
                self._d3plot.arrays[ArrayType.element_beam_node_indexes]  = element_beam_node_indexes.astype(int)
                self._d3plot.arrays[ArrayType.element_beam_part_indexes]  = element_beam_part_indexes.astype(int) + shell_part_num
            
                self._d3plot.arrays[ArrayType.element_beam_is_alive]      = np.array(element_beam_is_alive)
            
            
                #self._d3plot.arrays[ArrayType.element_beam_axial_force]  = d3.arrays["element_beam_axial_force"]
                #self._d3plot.arrays[ArrayType.element_beam_bending_moment] = d3.arrays["element_beam_bending_moment"]
                #self._d3plot.arrays[ArrayType.element_beam_shear_force] = d3.arrays["element_beam_shear_force"]
                #self._d3plot.arrays[ArrayType.element_beam_torsion_moment] = d3.arrays["element_beam_torsion_moment"]
            
                "Partless beams exists so we have to generate a dummy one 1000000000 - this could be improved"
            
                beam_part_ids                                       =   np.concatenate([np.array(rr.raw_arrays["pText1DA"])\
                                                                                        .astype("U9").astype(int),\
                                                                                        np.array([1000000000])])
                beam_part_titles                                    =   np.concatenate([np.array(rr.raw_arrays["pText1DA"]),
                                                                                        np.array([1000000000])])
                beam_part_num                                       =   len(beam_part_ids)
            
            else:
                beam_part_ids                                       =   []
                beam_part_titles                                    =   []
                beam_part_num                                       =   0
          
            "Solids"
            
            if rr.raw_header["nbElts3D"] > 0: 
                self._d3plot.arrays[ArrayType.element_solid_node_indexes] = rr.arrays["element_solid_node_indexes"]
                self._d3plot.arrays[ArrayType.element_solid_part_indexes] = rr.arrays["element_solid_part_indexes"] + shell_part_num + beam_part_num
            
                self._d3plot.arrays[ArrayType.element_solid_is_alive]      = np.array(element_solid_is_alive)
                self._d3plot.arrays[ArrayType.element_solid_ids]            = element_solid_ids_out
            
                #self._d3plot.arrays[ArrayType.element_solid_effective_plastic_strain] = d3.arrays["element_solid_effective_plastic_strain"]
                #self._d3plot.arrays[ArrayType.element_solid_ids] = d3.arrays["element_solid_ids"]
                #self._d3plot.arrays[ArrayType.element_solid_node_indexes] = rr.arrays["element_solid_node_indexes"]
                #self._d3plot.arrays[ArrayType.element_solid_part_indexes] = rr.arrays["element_solid_part_indexes"]
                #self._d3plot.arrays[ArrayType.element_solid_stress] = np.zeros(shape=(states, rr.raw_header["nbElts3D"], 9))
            
                solid_part_ids                                      =  np.array(rr.raw_arrays["pText3DA"]).astype("U9").astype(int)
                solid_part_titles                                    =  np.array(rr.raw_arrays["pText3DA"])
            
            else:
                solid_part_ids                                      =   []
                solid_part_titles                                   =   []
    
    
            "Global"
            #self._d3plot.arrays[ArrayType.global_internal_energy] = np.array([0] * states)
            #self._d3plot.arrays[ArrayType.global_kinetic_energy] = np.array([0] * states)
            #self._d3plot.arrays[ArrayType.global_total_energy] = np.array([0] * states)
            #self._d3plot.arrays[ArrayType.global_velocity] = np.array([[0, 0, 0]] * states).astype(float)
            
            "Parts"
    
            self._d3plot.arrays[ArrayType.part_ids]                   = np.concatenate([shell_part_ids, beam_part_ids, solid_part_ids])
            #self._d3plot.arrays[ArrayType.part_mass]                  = np.array([np.concatenate([np.array(part_mass), np.zeros(shape=(rr.raw_header["nbParts1D"] + 1)), np.zeros(shape=(rr.raw_header["nbParts3D"]))])] * states)
            
            #self._d3plot.arrays[ArrayType.part_hourglass_energy]      = np.array([np.concatenate([np.array(part_mass), np.zeros(shape=(rr.raw_header["nbParts1D"] + 1)), np.zeros(shape=(rr.raw_header["nbParts3D"]))])] * states)
            #self._d3plot.arrays[ArrayType.part_ids_cross_references] = d3.arrays["part_ids_cross_references"]
            #self._d3plot.arrays[ArrayType.part_ids_unordered] = np.concatenate([shell_part_ids, beam_part_ids, solid_part_ids])
            #self._d3plot.arrays[ArrayType.part_internal_energy] = np.array([np.concatenate([np.array(part_mass), np.zeros(shape=(rr.raw_header["nbParts1D"] + 1)), np.zeros(shape=(rr.raw_header["nbParts3D"]))])] * states)
            
            #self._d3plot.arrays[ArrayType.part_kinetic_energy] = np.array([np.concatenate([np.array(part_mass), np.zeros(shape=(rr.raw_header["nbParts1D"] + 1)), np.zeros(shape=(rr.raw_header["nbParts3D"]))])] * states)
           
            self._d3plot.arrays[ArrayType.part_titles] = np.array(np.concatenate([shell_part_titles, beam_part_titles, solid_part_titles])).astype('S')
            self._d3plot.arrays[ArrayType.part_titles_ids] = np.concatenate([shell_part_ids, beam_part_ids, solid_part_ids])
            #self._d3plot.arrays[ArrayType.part_velocity] = np.array([np.concatenate([np.array(part_mass), np.zeros(shape=(rr.raw_header["nbParts1D"] + 1)), np.zeros(shape=(rr.raw_header["nbParts3D"]))])] * states)
            #self._d3plot.arrays[ArrayType.rigid_wall_force] = d3.arrays["rigid_wall_force"]
    
            self._d3plot.arrays[ArrayType.global_timesteps]           = np.array(timesteps)
            
            
            
            random_name = readAndConvert.generate_random_name(10, ifile)
            temp_d3plot_name = os.path.dirname(file_stem) + random_name
            
            self._d3plot.write_d3plot(temp_d3plot_name, single_file = False)
            
            "Rename or remove header"
            if ifile == 0:
                _ = file_stem + ".d3plot"
                if os.path.isfile(_):
                    os.remove(_)
                os.rename(temp_d3plot_name, _)
            else:
                os.remove(temp_d3plot_name)
                
            "Rename plot state"
            d3plot_index = ifile+1
            if d3plot_index <100: 
                padding = 2
            else:
                padding = 3
                
            _ = file_stem + ".d3plot" + str(d3plot_index).zfill(padding)
            if os.path.isfile(_):
                os.remove(_)
            os.rename(temp_d3plot_name + "01", _)
                         
         
        return True
        

if __name__ == '__main__':        
    
    file_stem = "O:/model"

    a2d = readAndConvert(file_stem)
