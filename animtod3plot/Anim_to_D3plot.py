# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

from vortex_radioss.lasso.dyna.d3plot import D3plot
from vortex_radioss.lasso.dyna.array_type import ArrayType
import RadiossReader

import numpy as np
import os
import time
from tqdm import tqdm
import glob

class convert:
    @staticmethod
    def global_internal_energy(*data):
        # conversion logic for global_internal_energy
        return data

    @staticmethod
    def part_internal_energy(*data):
        # conversion logic for part_internal_energy
        return data

    @staticmethod
    def global_kinetic_energy(*data):
        # conversion logic for global_kinetic_energy
        return data

    @staticmethod
    def part_kinetic_energy(*data):
        # conversion logic for part_kinetic_energy
        return data

    @staticmethod
    def rigid_body_velocity(*data):
        # conversion logic for rigid_body_velocity
        return data
    
    @staticmethod
    def sum2Scalars(*data):
        # sum of 2 scalar results
        half_size = len(data) // 2
        data1 = data[:half_size]
        data2 = data[half_size:]
        print(data1)
        print(data2)
        return sum(data1) + sum(data2)
    
    @staticmethod
    def element_shell_internal_energy(*data):
            
        # data[0] is shell energy per unit mass (Radioss units)
        # data[1] is shell density
        # out is shell energy per unit volume (LS-Dyna units)

        out = np.multiply(data[0], data[1])
                     
        return out   
    
    @staticmethod
    def element_shell_stress(*data):
        
        # Mid-surface stresses if present are not converted
        # Only Upper and Lower stresses are converted 
        # Out of plane stresses not converted
        # Co-ordinate systems not corrected for - should be OK for Von-Mises, Tresca etc and \
        # stress principles
        # [σx, σy, σz, σxy, σyz, σxz]
        
        shell_num      = len(data[0])
        nip            = data[2][0] 
        
        out            = np.zeros(shape=(shell_num, nip, 6))              
        
        # Top integration point
        out[:, -1, 0]   = data[0][:, 0]
        out[:, -1, 1]   = data[0][:, 1]
        out[:, -1, 3]   = data[0][:, 2]
            
        # Bottom integration point
        out[:, 0, 0]   = data[1][:, 0]
        out[:, 0, 1]   = data[1][:, 1]
        out[:, 0, 3]   = data[1][:, 2]
        
        return out
    
    @staticmethod
    def element_shell_effective_plastic_strain(*data):
        
        # Mid-surface stresses if present are not converted
        # Only Upper and Lower stresses are converted         
        # [σx, σy, σz, σxy, σyz, σxz]
        
        shell_num      = len(data[0])
        nip            = data[2][0] 
        
        out            = np.zeros(shape=(shell_num, nip))
            
        # Top integration point
        out[:, -1]   = data[0]

        # Bottom integration point
        out[:, 0]    = data[1]
        
        return out

    @staticmethod

    def element_shell_history_vars(*data):

        shell_num = len(data[0])
        nip = data[2][0]

        out = np.zeros(shape=(shell_num, nip, 1))

        # Top integration point
        out[:, -1, 0] = data[0]

        # Bottom integration point
        out[:, 0, 0]   = data[1]

        return out

    def element_solid_history_vars(*data):

        solid_num = len(data[0])
        #nip = data[2][0]
        out = np.zeros(shape=(solid_num, 1, 1))
        out[:, -1, 0] = data[0]

        return out

class readAndConvert:
    
    def __init__(
        self,
        filepath: str = None,
        use_shell_mask=True,
        use_solid_mask=False,
        use_beam_mask=False,
        silent=False
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
        self.A_2_D(filepath, use_shell_mask, use_solid_mask, use_beam_mask, silent)

    def sequential(input_array):
        
        "IDs defined as zero need unique, non-zero renumbering"    
  
        zero_num_start = max(input_array).astype('int64')
        mask = (input_array == np.zeros(len(input_array)))
        new_id = np.cumsum(mask) + zero_num_start
        output = np.where(mask, new_id, input_array)
        
        return output 
    
    def generate_mask_map(input_array):
        
        unique, indices = np.unique(input_array, return_index=True)
        
        if unique[0]:
            mask_map = indices
        else:
            mask_map = indices[1:]
        
        return np.array(mask_map).astype(int)
    
    def generate_display_entity(indexes, names, search_str_arr, show):
        
        if show:
            _ = np.zeros(shape=len(names))
            __ = 1
        else:
            _ = np.ones(shape=len(names))
            __ = 0

        for search_str in search_str_arr:       
            ___ = np.flatnonzero(np.core.defchararray.find(names,search_str)!=-1)

            _[___] = __

        display_array = np.array(_[indexes]).astype(bool)
        
        return display_array           

    def generate_sorter(input_array):
           
       return np.argsort(input_array)
   
    def invert_sorter(input_array):
        index_tracker = np.zeros((np.max(input_array) + 1,))
        index_tracker[input_array] = np.arange(0, input_array.size)    
           
        return np.array(index_tracker).astype(int)  
   
    def apply_sorter(input_array, tracker_array):
        
        return input_array[tracker_array]   
    
    def generate_random_name(length, ifile):
        import random
        import string
        # choose from all lowercase letter
        letters = string.ascii_lowercase
        result_str = ''.join(random.choice(letters) for i in range(length))
        return "_" + result_str + str(ifile)    
    
    def LOGGER(self, string, silent):
        if not silent:
            print(string)
    

    def A_2_D(self,file_stem, use_shell_mask, use_solid_mask, use_beam_mask, silent):
    
        if os.path.isfile(file_stem + "d3plot"):
            return True

        file_list = glob.glob(file_stem + "A*[0-9]")
        file_list.sort()
        
        original_node_coordinates       = None
                       
        self.LOGGER("Mapping database", silent)
        
        if not file_list:
            print("No files found..\nPlease check file stem:\ne.g\nC:/Folder/ModelA00*\nWould be:\nC:/Folder/Model*")
            return
        else:
            rr = RadiossReader.RadiossReader(file_list[0])  
        
        #rr.raw_header["nbElts1D"] = 0
        
        "We can optionally generate mask arrays that filter out elements that have an ID of zero - These are non-structural and are generated by Radioss for visualisations"
        "We sort the ids into ascending order and generate a mapping function for the scalar and tensor data"

        n_nodes = rr.raw_header["nbNodes"]
        
        n_shell = 0
        nip_shell = 2
        n_solid = 0
        n_solid_layers = 1
        if rr.raw_header["nbFacets"] > 0:
             
            shell_ids_tracker = readAndConvert.generate_sorter(readAndConvert.sequential(rr.arrays["element_shell_ids"]))
            
            if use_shell_mask:
                
                # Masking can remove useful entities such as rigid walls
                # This will search a text array to match a string and return the boolean array to indicate should be shown
                # Changing True to False will generate the inverse array
                part_strings = ["_rigid_wall_",  ": RIGIDWALL_"]
                display_shell = readAndConvert.generate_display_entity(readAndConvert.apply_sorter(rr.arrays["element_shell_part_indexes"], shell_ids_tracker), rr.raw_arrays["pTextA"], part_strings, True) 

                # The shell mask function filters out any entity that has a id of zero
                # In order to keep entities from the bool array we employ a trick below
                # we create dummy unique ids that are non-zero for all elements that are active in the display shell array
                # This dummy array is then given to the masking function

                _ = max(rr.arrays["element_shell_ids"])
                dummy_id_offset = (np.cumsum(display_shell) + _)*display_shell.astype(bool)
                shell_mask = readAndConvert.generate_mask_map(readAndConvert.apply_sorter(rr.arrays["element_shell_ids"], shell_ids_tracker) + dummy_id_offset)
    
            else:
                shell_mask = readAndConvert.apply_sorter(np.arange(0, len(rr.arrays["element_shell_ids"])), shell_ids_tracker)              
            
            
            shell_ids_tracker       = shell_ids_tracker[shell_mask]
            
            element_shell_ids_out   =   readAndConvert.apply_sorter(readAndConvert.sequential(rr.arrays["element_shell_ids"]), shell_ids_tracker)
            
            n_shell = len(element_shell_ids_out)
            


        if rr.raw_header["nbElts1D"] > 0:                                         
            
            beam_ids_tracker = readAndConvert.generate_sorter(readAndConvert.sequential(rr.arrays["element_beam_ids"]))
            element_beam_ids_out    =   readAndConvert.apply_sorter(readAndConvert.sequential(rr.arrays["element_beam_ids"]), beam_ids_tracker)

        if rr.raw_header["nbElts3D"] > 0:
            
            solid_ids_tracker = readAndConvert.generate_sorter(readAndConvert.sequential(rr.arrays["element_solid_ids"]))
            element_solid_ids_out   =   readAndConvert.apply_sorter(readAndConvert.sequential(rr.arrays["element_solid_ids"]), solid_ids_tracker)
            n_solid = len(element_solid_ids_out)
            n_solid_layers = 1
               
        _                                                           = readAndConvert.sequential(rr.arrays["node_ids"])
        node_ids_tracker                                            = readAndConvert.generate_sorter(_)                
        inverted_node_ids_tracker                                   = readAndConvert.invert_sorter(node_ids_tracker)
        original_node_coordinates                                   = readAndConvert.apply_sorter(rr.arrays["node_coordinates"], node_ids_tracker)
        
        node_ids_out = readAndConvert.apply_sorter(_, node_ids_tracker)
        
        self._d3plot.arrays[ArrayType.node_ids]                       = node_ids_out        
        
        "Parts and geometry"

        # We need to order the PIDs in numerical order - we generate tracking functions to map indexes correctly 
        if rr.raw_header["nbFacets"] > 0:
            
            shell_part_ids                                              =  np.array(rr.raw_arrays["pTextA"]).astype("U9").astype(int)
            shell_part_names                                            =  np.char.strip(np.array(list(np.char.split(rr.raw_arrays["pTextA"], sep=":")))[:,1])
            shell_part_num                                              = len(shell_part_ids)
        
        else:
            shell_part_ids                                              = []
            shell_part_names                                            = []
            shell_part_num                                              = 0
            
        if rr.raw_header["nbElts3D"] > 0:            
            element_solid_node_indexes = rr.arrays["element_solid_node_indexes"]            
            for i_solid, solid in enumerate(rr.arrays["element_solid_node_indexes"]):
                # Correct for tetra node connectivity
                _ = len(np.unique(solid))
                if _ == 4:
                    __ = np.array([solid[0], solid[2], solid[4], solid[5], solid[5], solid[5], solid[5], solid[5]])
                    if len(np.unique(__)) == 4:
                        element_solid_node_indexes[i_solid] = __
            solid_part_ids                                       =  np.array(rr.raw_arrays["pText3DA"]).astype("U9").astype(int)
            solid_part_names                                     =  np.char.strip(np.array(list(np.char.split(rr.raw_arrays["pText3DA"], sep=":")))[:,1])
        
        else:
            solid_part_ids                                      =   []             
            solid_part_names                                    =   []   
                        
        if rr.raw_header["nbElts1D"] > 0: 
           
            additional_beam_number                                      = rr.raw_header["nbElts1D"] - len(rr.arrays["element_beam_part_indexes"])
            additional_beams                                            = np.zeros(shape=(additional_beam_number))
            additional_beams.fill(max(rr.arrays["element_beam_part_indexes"])+1)
            element_beam_part_indexes                                   = np.concatenate([rr.arrays["element_beam_part_indexes"], additional_beams])
            self._d3plot.arrays[ArrayType.element_beam_ids]             = element_beam_ids_out
            element_beam_node_indexes=np.zeros(shape=(rr.raw_header["nbElts1D"],5))
            element_beam_node_indexes[:,0]                              = rr.arrays["element_beam_node_indexes"][:,0]
            element_beam_node_indexes[:,1]                              = rr.arrays["element_beam_node_indexes"][:,1]
            
            
            "Partless beams exists so we have to generate a dummy one"
            
            dummy_beam_pid = max(np.concatenate([shell_part_ids, np.array(rr.raw_arrays["pText1DA"]).astype("U9").astype(int), solid_part_ids])) +1
            beam_part_ids                                       =   np.concatenate([np.array(rr.raw_arrays["pText1DA"]).astype("U9").astype(int),np.array([dummy_beam_pid])])
            beam_part_num                                       =   len(beam_part_ids) 
           
            beam_part_names                                     =   np.array(list(np.char.split(rr.raw_arrays["pText1DA"], sep=":", maxsplit=1)))
            beam_part_names                                     =   np.char.strip(np.concatenate([beam_part_names[:,1], np.array(["ADDITIONAL VISUALISATION BEAMS"])]))
        else:
            beam_part_ids                                       =   []
            beam_part_names                                     =   []
            beam_part_num                                       =   0              
                      
# Generate Part trackers
        # All part ids
        _                                                           = np.concatenate([shell_part_ids, beam_part_ids, solid_part_ids])
        # All part names 
        __                                                          = np.concatenate([shell_part_names, beam_part_names, solid_part_names])        
        # Set all parts with invalid names' PIDS to zero
        invalid_name_strings = ["_rigid_wall_",  ": RIGIDWALL_"]
        ____ = np.array([])
        for search_str in invalid_name_strings:       
            ___ = np.flatnonzero(np.core.defchararray.find(__,search_str)!=-1)
            _[___] = 0
            ____ = np.concatenate([____,___]).astype(int)

        # Generate new unique PIDS with smallest value
        currently_valid_pids = (set(_))
        possible_pids = set(np.arange(len(_)+1))
        available_pids = list(possible_pids.difference(currently_valid_pids))
        available_pids.sort()
        # Assign and update the new PIDS
        #print(available_pids)
        #print(____)
        _[____] = available_pids[:len(____)]
        

        part_ids_tracker                                            = readAndConvert.generate_sorter(_)
        inverted_part_ids_tracker                                   = readAndConvert.invert_sorter(part_ids_tracker)
        
        # Part arrays these all need defining
        self._d3plot.arrays[ArrayType.part_ids]                     = readAndConvert.apply_sorter(_, part_ids_tracker).astype(int)
        #self._d3plot.arrays[ArrayType.part_ids_unordered]           = readAndConvert.apply_sorter(_, part_ids_tracker).astype(int) 
        
        #self._d3plot.arrays[ArrayType.part_ids]                     = np.arange(1, len(_)+1)
        # These are not IDS, this is a FORTRAN INDEX System
        self._d3plot.arrays[ArrayType.part_ids_unordered]           = readAndConvert.apply_sorter(_, part_ids_tracker).astype(int) 
        
        # Part titles are encoded into ascii binary
        __                                                          = np.array(readAndConvert.apply_sorter(__, part_ids_tracker))
        # Replace some strings with new strings
        __[np.where(__=="rby")[0]] = "Rigid elements"        
        
        #__ = ["a" for i in range(0,len(__))]
        
        self._d3plot.arrays[ArrayType.part_titles]                  = np.char.encode(np.array(__), encoding = "UTF-8")
        # Part masses are floats
        self._d3plot.arrays[ArrayType.part_mass]                    = np.array([np.zeros(len(__))]).astype("<f")
        self._d3plot.arrays[ArrayType.part_hourglass_energy]                    = np.array([np.zeros(len(__))]).astype("<f")
        

        self._d3plot.arrays[ArrayType.part_titles_ids]              = readAndConvert.apply_sorter(_, part_ids_tracker).astype(int)
        # Not an ID this is an index in Fortran starting at 1
        self._d3plot.arrays[ArrayType.part_ids_cross_references]    = np.arange(1,len(_) +1).astype("int64")
        
        if rr.raw_header["nbFacets"] > 0:
            # Assign the shell part indexes
            _   = readAndConvert.apply_sorter(rr.arrays["element_shell_part_indexes"], shell_ids_tracker)
            self._d3plot.arrays[ArrayType.element_shell_part_indexes]     = inverted_part_ids_tracker[_]
            self._d3plot.arrays[ArrayType.element_shell_ids]              = element_shell_ids_out.astype(int) 
            self._d3plot.arrays[ArrayType.element_shell_node_indexes]     = inverted_node_ids_tracker[readAndConvert.apply_sorter(rr.arrays["element_shell_node_indexes"], shell_ids_tracker).astype(int)  ]
        
        if rr.raw_header["nbElts1D"] > 0:
            # Assign the beam indexes
            _   = readAndConvert.apply_sorter(element_beam_part_indexes, beam_ids_tracker).astype(int) + shell_part_num            
            self._d3plot.arrays[ArrayType.element_beam_part_indexes]    = inverted_part_ids_tracker[_]
            self._d3plot.arrays[ArrayType.element_beam_node_indexes]    = inverted_node_ids_tracker[readAndConvert.apply_sorter(element_beam_node_indexes, beam_ids_tracker).astype(int)]
        
        if rr.raw_header["nbElts3D"] > 0:            
            # Assign the solid part indexes
            _   = readAndConvert.apply_sorter(rr.arrays["element_solid_part_indexes"], solid_ids_tracker) + shell_part_num + beam_part_num              
            self._d3plot.arrays[ArrayType.element_solid_part_indexes]   = inverted_part_ids_tracker[_]   
            self._d3plot.arrays[ArrayType.element_solid_ids]            = element_solid_ids_out.astype(int)  
            self._d3plot.arrays[ArrayType.element_solid_node_indexes]   = inverted_node_ids_tracker[readAndConvert.apply_sorter(element_solid_node_indexes, solid_ids_tracker)]
            
        self.LOGGER("Processing states", silent)
        
        for ifile, file in enumerate(tqdm(file_list, disable = silent)):            
            
            if ifile:
                rr = RadiossReader.RadiossReader(file) 
                #rr.raw_header["nbElts1D"] = 0
                
            
            self._d3plot.arrays[ArrayType.node_coordinates]               = original_node_coordinates          
            
            "Node Updates"
                                                            
            self._d3plot.arrays[ArrayType.node_displacement]              = np.array([readAndConvert.apply_sorter(rr.arrays["node_coordinates"], node_ids_tracker)]).astype("<f")
               
            "Timestep Updates"
            timesteps                       = []              
            timesteps.append(rr.arrays["timesteps"])            
            
            "Shells"    
            if rr.raw_header["nbFacets"] > 0:
                element_shell_is_alive                                        = []
                element_shell_is_alive.append(readAndConvert.apply_sorter(rr.arrays["element_shell_is_alive"], shell_ids_tracker).astype("<f"))
                self._d3plot.arrays[ArrayType.element_shell_is_alive]         = np.array(element_shell_is_alive).astype("<f") 
                            
            "Beams"
            
            if rr.raw_header["nbElts1D"] > 0: 
                element_beam_is_alive           = []
                element_beam_is_alive.append(readAndConvert.apply_sorter(rr.arrays["element_beam_is_alive"], beam_ids_tracker).astype("<f"))
                self._d3plot.arrays[ArrayType.element_beam_is_alive]        = np.array(element_beam_is_alive).astype("<f") 

            "Solids"
            
            if rr.raw_header["nbElts3D"] > 0: 
                element_solid_is_alive                                      = []
                element_solid_is_alive.append(readAndConvert.apply_sorter(rr.arrays["element_solid_is_alive"], solid_ids_tracker))
                self._d3plot.arrays[ArrayType.element_solid_is_alive]       = np.array(element_solid_is_alive).astype("<f")                         
            
            "Global"
            self._d3plot.arrays[ArrayType.global_timesteps]             = np.array(timesteps)
                                                           
            if True:
                
                database_extent_binary ={}
                array_requirements = {}

                if rr.raw_header["nbNodes"] > 0:   
                    
                    flag = "NODES"
                    
                    database_extent_binary[flag] = {}      
                    _ = database_extent_binary[flag]
                    
                    database_extent_binary[flag][0] = []
                    database_extent_binary[flag][1] = _[0] + [ArrayType.node_velocity]  
                    database_extent_binary[flag][2] = _[1] + [ArrayType.node_acceleration]  
                    
                    # Dyna output
                    array_requirements[ArrayType.node_velocity] = {}
                    _ = array_requirements[ArrayType.node_velocity]
                    # Radioss outputs needed to compute Dyna output
                    _["dependents"]     = ["node_velocity"]
                    _["shape"]          = (1,n_nodes,3)
                    _["convert"]        = None
                    _["tracker"]        = node_ids_tracker
                    _["additional"]     = []
                    
                    # Dyna output
                    array_requirements[ArrayType.node_acceleration] = {}
                    _ = array_requirements[ArrayType.node_acceleration]
                    # Radioss outputs needed to compute Dyna output
                    _["dependents"]     = ["node_acceleration"]
                    _["shape"]          = (1,n_nodes,3)
                    _["convert"]        = None
                    _["tracker"]        = node_ids_tracker
                    _["additional"]     = []
                
                if rr.raw_header["nbEFunc"] > 0:   
                    
                    flag = "SHELLS"
                    
                    database_extent_binary[flag] = {}      
                    _ = database_extent_binary[flag]
                    
                    # [0] are essential arrays and need creating even if no data available
                    
                    database_extent_binary[flag][0] = []
                    database_extent_binary[flag][0] = _[0] + [ArrayType.element_shell_is_alive]  
                    database_extent_binary[flag][0] = _[0] + [ArrayType.element_shell_stress]    
                    database_extent_binary[flag][0] = _[0] + [ArrayType.element_shell_effective_plastic_strain]
                    database_extent_binary[flag][0] = _[0] + [ArrayType.element_shell_history_vars]

                    database_extent_binary[flag][1] = _[0] + [ArrayType.element_shell_thickness]        
                    database_extent_binary[flag][1] = _[1] + [ArrayType.element_shell_internal_energy]                        
                    
                    # Dyna output
                    array_requirements[ArrayType.element_shell_thickness] = {}
                    _ = array_requirements[ArrayType.element_shell_thickness]
                    # Radioss outputs needed to compute Dyna output
                    _["dependents"]     = ["element_shell_thickness"]
                    _["shape"]          = (1,n_shell)
                    _["convert"]        = None
                    _["tracker"]        = shell_ids_tracker
                    _["additional"]     = []
                    
                    # Dyna output
                    array_requirements[ArrayType.element_shell_is_alive] = {}
                    _ = array_requirements[ArrayType.element_shell_is_alive]
                    # Radioss outputs needed to compute Dyna output
                    _["dependents"]     = ["element_shell_is_alive"]
                    _["shape"]          = (1,n_shell)
                    _["convert"]        = None
                    _["tracker"]        = shell_ids_tracker
                    _["additional"]     = []
                    
                    # Dyna output
                    array_requirements[ArrayType.element_shell_internal_energy] = {}
                    _ = array_requirements[ArrayType.element_shell_internal_energy]
                    # Radioss outputs needed to compute Dyna output
                    _["dependents"]     = ["element_shell_specific_energy", "element_shell_density"]
                    _["shape"]          = (1,n_shell)
                    _["convert"]        = convert.element_shell_internal_energy
                    _["tracker"]        = shell_ids_tracker
                    _["additional"]     = []
                    
                   # Dyna output
                    array_requirements[ArrayType.element_shell_stress] = {}
                    _ = array_requirements[ArrayType.element_shell_stress]
                    # Radioss outputs needed to comptute Dyna output
                    _["dependents"]     = ["Stress (upper)","Stress (lower)"]
                    _["shape"]          = (1,n_shell, nip_shell, 6)
                    _["convert"]        = convert.element_shell_stress
                    _["tracker"]        = shell_ids_tracker
                    _["additional"]     = [nip_shell]
                    
                   # Dyna output
                    array_requirements[ArrayType.element_shell_effective_plastic_strain] = {}
                    _ = array_requirements[ArrayType.element_shell_effective_plastic_strain]
                    # Radioss outputs needed to comptute Dyna output
                    _["dependents"]     = ["element_shell_plastic_strain_upper", 'element_shell_plastic_strain_lower']
                    _["shape"]          = (1, n_shell, nip_shell)
                    _["convert"]        = convert.element_shell_effective_plastic_strain
                    _["tracker"]        = shell_ids_tracker  
                    _["additional"]     = [nip_shell]

                   # Dyna output
                    array_requirements[ArrayType.element_shell_history_vars] = {}
                    _ = array_requirements[ArrayType.element_shell_history_vars]
                    # Radioss outputs needed to comptute Dyna output
                    _["dependents"]     = ["element_shell_damage,(layer___1)__", "element_shell_damage,(layer___5)__"]
                    _["shape"]          = (1, n_shell, nip_shell, 1)
                    _["convert"]        = convert.element_shell_history_vars
                    _["tracker"]        = shell_ids_tracker
                    _["additional"]     = [nip_shell]
                try:
                    if rr.raw_header["nbEFunc3D"] > 0:
                        flag = "SOLIDS"

                        database_extent_binary[flag] = {}
                        _ = database_extent_binary[flag]

                        # [0] are essential arrays and need creating even if no data available

                        database_extent_binary[flag][0] = []
                        # database_extent_binary[flag][0] = _[0] + [ArrayType.element_solid_is_alive]
                        database_extent_binary[flag][0] = _[0] + [ArrayType.element_solid_history_variables]

                        array_requirements[ArrayType.element_solid_history_variables] = {}
                        _ = array_requirements[ArrayType.element_solid_history_variables]
                        # Radioss outputs needed to compute Dyna output
                        _["dependents"]     = ["element_solid_max_damage_element"]
                        _["shape"]          = (1,n_solid,n_solid_layers,1)
                        _["convert"]        = convert.element_solid_history_vars
                        _["tracker"]        = solid_ids_tracker
                        _["additional"]     = []
                except:
                    pass
            "Assign the arrays to the D3PLOT class for writing"
                        
            "Generate the availability check"
            
            dependency_check = {}
            flag_max        =  {}

            for flag in database_extent_binary:
                flag_max[flag] = -float('inf')
                _ = database_extent_binary[flag]
                for array_group in _:
                    __ = _[array_group]
                    for array_output in __:
                        # Check all dependencies exist
                        all_dependents_exist = True
                        ___ = array_requirements[array_output]["dependents"]
                        for check_dependent in ___:
                            if check_dependent not in rr.arrays:
                                all_dependents_exist = False
                                break
                    
                        dependency_check[array_output] = all_dependents_exist
                        if all_dependents_exist:
                            flag_max[flag] = max(flag_max[flag],array_group) 
            
            "Generate the output arrays"
  
            for flag in database_extent_binary:
                _ = database_extent_binary[flag]
                for array_group in _:
                    __ = _[array_group]
                    
                    if array_group <= flag_max[flag]:
                        for array_output in __:
                            if dependency_check[array_output]:
                                conversion_function = array_requirements[array_output]["convert"]
                                tracker             = array_requirements[array_output]["tracker"]
                                if conversion_function:
                                    conversion_inputs   = [rr.arrays[i] for i in array_requirements[array_output]["dependents"]]
                                    additional_inputs   = [i for i in array_requirements[array_output]["additional"]]
                                    ___                   = conversion_function(*conversion_inputs, additional_inputs)
                                    if tracker is not None:
                                        self._d3plot.arrays[array_output] = np.array(readAndConvert.apply_sorter(___, tracker)[np.newaxis, :]).astype("<f")
                                    else:
                                        self._d3plot.arrays[array_output] = np.array(___[np.newaxis, :]).astype("<f")
                                else:
                                    ___ = rr.arrays[array_requirements[array_output]["dependents"][0]]
                                    if tracker is not None:
                                        self._d3plot.arrays[array_output] = np.array(readAndConvert.apply_sorter(___, tracker)[np.newaxis, :]).astype("<f")
                                    else:
                                        self._d3plot.arrays[array_output] = np.array(___[np.newaxis, :]).astype("<f")
                            else:
                                self._d3plot.arrays[array_output] = np.zeros(shape=(array_requirements[array_output]["shape"])).astype("<f")
                                                                            
            "Write out the d3plot state"
                                        
            random_name = readAndConvert.generate_random_name(10, ifile)
            temp_d3plot_name = os.path.dirname(file_stem) + "/" + random_name
            
            d3plot_index = ifile+1
            if d3plot_index <100: 
                padding = 2
            else:
                padding = 3            
            
            "Remove old files"
            if ifile == 0:
                _ = file_stem + ".d3plot"
                if os.path.isfile(_):
                    os.remove(_)
            _ = file_stem + ".d3plot" + str(d3plot_index).zfill(padding)
            if os.path.isfile(_):
                os.remove(_)           
            
            self._d3plot.header.itype = np.int64
            self._d3plot.header.ftype = np.float64
            self._d3plot.header.wordsize = 8

            
            self._d3plot.write_d3plot(temp_d3plot_name, single_file = False)
            
            "Rename new files and remove temp files"
            if ifile == 0:
                _ = file_stem + ".d3plot"
                os.rename(temp_d3plot_name, _)
            else:
                os.remove(temp_d3plot_name)
                
            _ = file_stem + ".d3plot" + str(d3plot_index).zfill(padding)
            os.rename(temp_d3plot_name + "01", _)        
                                                 
        return True
        
if __name__ == '__main__':   
             
    #file_stem = "C:/Users/PC/Downloads/test/DynaOpt"    
    #file_stem = "P:/Optimisation/A001a/Baseline/FFB_0/DynaOpt"
    file_stem = "./yaris"
    #file_stem = "C:/Users/PC/Downloads/test4/CRA2AV4"
    
    a2d = readAndConvert(file_stem, use_shell_mask=False)
    
      
