# -*- coding: utf-8 -*-
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

from lasso.dyna.d3plot import D3plot
from lasso.dyna.array_type import ArrayType
#from .RadiossTHReader import RadiossTHReader
import RadiossReader

import numpy as np
import os
import time
from tqdm import tqdm
import glob

#########################################
# New part
#########################################
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
    def element_shell_stress(*data):
            
        #print(data)

        #data1[:, 0, 0] = data[0][:, 0]
        #data1[:, 0, 1] = data[0][:, 1]
        #data1[:, 0, 2] = np.zeros(shape=(rr.raw_header["nbFacets"]))
        #data1[:, 0, 3] = data[0][:, 2]
        #data1[:, 0, 5] = np.zeros(shape=(rr.raw_header["nbFacets"]))    
            
            
        #data1[:, 1, 0] = data[1][:, 0]
        #data1[:, 1, 1] = data[1][:, 1]
        #data1[:, 1, 2] = np.zeros(shape=(rr.raw_header["nbFacets"]))
        #data1[:, 1, 3] = data[1][:, 2]
        #data1[:, 1, 5] = np.zeros(shape=(rr.raw_header["nbFacets"]))  
             
        return data
                
#########################################


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

        file_list = glob.glob(file_stem + "A*[0-9]")
        file_list.sort()
        
        original_node_coordinates       = None
        
               
        print("Re-sorting the indexes")
        rr = RadiossReader.RadiossReader(file_list[0]) 

        
        n_parts = rr.raw_header["nbParts"]
        n_nodes = rr.raw_header["nbNodes"]
        n_beams = rr.raw_header["nbElts1D"]
        n_shell = rr.raw_header["nbFacets"]
        n_solid = rr.raw_header["nbElts3D"]
        #n_contact = rr.raw_header["nbContact"]
        #n_rigid_body = rr.raw_header["nbRigidBodies"]
        #nip_beams = rr.raw_header["nbIPBeam"]
        #nip_shell = rr.raw_header["nbIPShell"]
        #nip_solid = rr.raw_header["nbIPSolid"]

        if rr.raw_header["nbFacets"] > 0:
            shell_ids_tracker = readAndConvert.generate_sorter(readAndConvert.sequential(rr.arrays["element_shell_ids"]))
            element_shell_ids_out   =   readAndConvert.apply_sorter(readAndConvert.sequential(rr.arrays["element_shell_ids"]), shell_ids_tracker)

        if rr.raw_header["nbElts1D"] > 0:
            beam_ids_tracker = readAndConvert.generate_sorter(readAndConvert.sequential(rr.arrays["element_beam_ids"]))
            element_beam_ids_out    =   readAndConvert.apply_sorter(readAndConvert.sequential(rr.arrays["element_beam_ids"]), beam_ids_tracker)

        if rr.raw_header["nbElts3D"] > 0:
            solid_ids_tracker = readAndConvert.generate_sorter(readAndConvert.sequential(rr.arrays["element_solid_ids"]))
            element_solid_ids_out   =   readAndConvert.apply_sorter(readAndConvert.sequential(rr.arrays["element_solid_ids"]), solid_ids_tracker)
            

        original_node_coordinates = rr.arrays["node_coordinates"]       
        node_ids_out = readAndConvert.sequential(rr.arrays["node_ids"])
        
        del rr
        
        for ifile, file in enumerate(tqdm(file_list)):
            
            rr = RadiossReader.RadiossReader(file) 
            
            self._d3plot.arrays[ArrayType.node_coordinates]               = original_node_coordinates          
            
            "Node Updates"
 
            self._d3plot.arrays[ArrayType.node_displacement]              = np.array([rr.arrays["node_coordinates"]])            
            self._d3plot.arrays[ArrayType.node_ids]                       = node_ids_out
            
            if False:
                
                
                database_extent_binary ={}
                array_requirements = {}
                
                IGLB = False
                
                if IGLB:                
                    """
                    Output flag for global and part history variables such as internal
                    energy, kinetic energy, rigid body velocity, etc., that is, the sort of
                    data that can also be output to glstat and matsum.
                    
                    """
                    
                    # There are more but unsure of names...needs adding from an example d3plot
                    
                    flag = "IGLB"
                    
                    database_extent_binary[flag] = {} 
                    _ = database_extent_binary[flag]
                    _[0] = set([])
                    _[1] = _[0].add([   ArrayType.global_internal_energy,
                                        ArrayType.part_internal_energy,
                                        ArrayType.global_kinetic_energy,
                                        ArrayType.part_kinetic_energy,
                                        ArrayType.rigid_body_velocity])
                    
                    array_requirements[ArrayType.global_internal_energy] = {}
                    _ = array_requirements[ArrayType.global_internal_energy]
                    _["dependents"]     = ["global_internal_energy"]
                    _["shape"]          = (1,)
                    _["convert"]        = None
                    _["tracker"]        = None
                    
                    array_requirements[ArrayType.part_internal_energy] = {}
                    _ = array_requirements[ArrayType.part_internal_energy]
                    _["dependents"]     = ["part_internal_energy"]
                    _["shape"]          = (1,n_parts)   
                    _["convert"]        = convert.part_internal_energy
                    _["tracker"]        = None
                    
                    
                    array_requirements[ArrayType.global_kinetic_energy] = {}
                    _ = array_requirements[ArrayType.global_kinetic_energy]
                    _["dependents"]     = ["global_kinetic_energy"]
                    _["shape"]          = (1,)
                    _["convert"]        = convert.global_kinetic_energy
                    _["tracker"]        = None
                                                            
                    array_requirements[ArrayType.part_kinetic_energy] = {}
                    _ = array_requirements[ArrayType.part_kinetic_energy]
                    _["dependents"]     = ["part_kinetic_energy"]
                    _["shape"]          = (1,n_parts)    
                    _["convert"]        = convert.part_kinetic_energy
                    _["tracker"]        = None                    
                    
                    
                    array_requirements[ArrayType.rigid_body_velocity] = {}
                    _ = array_requirements[ArrayType.rigid_body_velocity]
                    _["dependents"]     = ["rigid_body_velocity"]
                    _["shape"]          = (1,n_rigid_body)     
                    _["convert"]        = convert.rigid_body_velocity
                    _["tracker"]        = None                       
               
                IXYZ = True
                
                if IXYZ:
                    """
                    Output flag for nodal coordinates
                    """
                    
                    database_extent_binary["IXYZ"]      = {}            
                    _                                   = database_extent_binary["IXYZ"]
                    database_extent_binary["IXYZ"][0]   = []
                    database_extent_binary["IXYZ"][1]   = _[0] + [ArrayType.node_displacement
                                                                       ]            
                
                    array_requirements[ArrayType.node_displacement] = {}
                    _ = array_requirements[ArrayType.node_displacement]
                    
                    _["dependents"]     = ["node_coordinates"]
                    _["shape"]          = (rr.raw_header["nbNodes"],3)
                    _["convert"]        = None
                    _["tracker"]        = None
                    
                
                IVEL = True
                
                if IVEL:
                    """
                    Output flag for nodal velocity
                    """
                    
                    flag = "IVEL"
                    
                    database_extent_binary[flag] = {}      
                    _ = database_extent_binary[flag]
                    
                    database_extent_binary[flag][0] = []
                    database_extent_binary[flag][1] = _[0] + [ArrayType.node_velocity
                                                                   ]           
                
                    array_requirements[ArrayType.node_velocity] = {}
                    _ = array_requirements[ArrayType.node_velocity]
                    _["dependents"]     = ["node_velocity"]
                    _["shape"]          = (rr.raw_header["nbNodes"],3)
                    _["convert"]        = None
                    _["tracker"]        = None         
                    
               
                IACC = True 
                if IACC:
                     """
                     Output flag for nodal acceleration
                     """
                     
                     flag = "IACC"
                     
                     database_extent_binary[flag] = {}      
                     _ = database_extent_binary[flag]
                     
                     database_extent_binary[flag][0] = []
                     database_extent_binary[flag][1] = _[0] + [ArrayType.node_acceleration
                                                                    ]          
                 
                     array_requirements[ArrayType.node_acceleration] = {}
                     _ = array_requirements[ArrayType.node_acceleration]
                     _["dependents"]     = ["node_acceleration"]
                     _["shape"]          = (rr.raw_header["nbNodes"],3)
                     _["convert"]        = None
                     _["tracker"]        = None        
                     
                ISTRS = False
                
                if ISTRS:
                     """
                     Output flag for stress tensor and “plastic strain”
                     """
                     
                     flag = "ISTRS"
                     
                     database_extent_binary[flag] = {}      
                     _ = database_extent_binary[flag]
                     
                     database_extent_binary[flag][0] = []
                     database_extent_binary[flag][1] = _[0] +      [    ArrayType.element_beam_effective_plastic_strain,
                                                                        ArrayType.element_shell_effective_plastic_strain,
                                                                        ArrayType.element_solid_effective_plastic_strain,
                                                                        ArrayType.element_beam_stress,
                                                                        ArrayType.element_shell_stress,
                                                                        ArrayType.element_solid_stress]          
                               
                 
                     array_requirements[ArrayType.element_beam_effective_plastic_strain] = {}
                     _ = array_requirements[ArrayType.element_beam_effective_plastic_strain]
                     _["dependents"]     = ["beam_plastic_strain_tensor"]
                     _["shape"]          = (1,n_beams, nip_beam)
                     _["convert"]        = convert.beam_plastic_strain_tensor
                     _["tracker"]        = beam_ids_tracker      
                     
                     array_requirements[ArrayType.element_shell_effective_plastic_strain] = {}
                     _ = array_requirements[ArrayType.element_shell_effective_plastic_strain]
                     _["dependents"]     = ["shell_plastic_strain_tensor"]
                     _["shape"]          = (1,n_shell, nip_shell)
                     _["convert"]        = convert.shell_plastic_strain_tensor
                     _["tracker"]        = shell_ids_tracker                   
                
                     array_requirements[ArrayType.element_solid_effective_plastic_strain] = {}
                     _ = array_requirements[ArrayType.element_solid_effective_plastic_strain]
                     _["dependents"]     = ["solid_plastic_strain_tensor"]
                     _["shape"]          = (1,n_solid, nip_solid)
                     _["convert"]        = convert.solid_plastic_strain_tensor
                     _["tracker"]        = solid_ids_tracker              
                
                     array_requirements[ArrayType.element_beam_stress] = {}
                     _ = array_requirements[ArrayType.element_beam_stress]
                     _["dependents"]     = ["beam_stress_tensor"]
                     _["shape"]          = (1,n_beams, nip_beam, 3)
                     _["convert"]        = convert.beam_stress_tensor
                     _["tracker"]        = beam_ids_tracker      
                     
                     array_requirements[ArrayType.element_shell_stress] = {}
                     _ = array_requirements[ArrayType.element_shell_stress]
                     _["dependents"]     = ["shell_stress_tensor"]
                     _["shape"]          = (1,n_shell, nip_shell, 4)
                     _["convert"]        = convert.shell_stress_tensor
                     _["tracker"]        = shell_ids_tracker                   
                
                     array_requirements[ArrayType.element_solid_stress] = {}
                     _ = array_requirements[ArrayType.element_solid_stress]
                     _["dependents"]     = ["solid_stress_tensor"]
                     _["shape"]          = (1,n_solid, nip_solid, 6)
                     _["convert"]        = convert.solid_stress_tensor
                     _["tracker"]        = solid_ids_tracker            
                     
                
                ISTRA = False
                if ISTRA:
                     """
                     Output flag for strain tensor:
                     """
                     
                     flag = "ISTRA"
                     
                     database_extent_binary[flag] = {}      
                     _ = database_extent_binary[flag]
                     
                     database_extent_binary[flag][0] = []
                     database_extent_binary[flag][1] = _[0] +      [    ArrayType.element_beam_strain,
                                                                        ArrayType.element_shell_strain,
                                                                        ArrayType.element_solid_strain
                                                                        ]
                
                     array_requirements[ArrayType.element_beam_strain] = {}
                     _ = array_requirements[ArrayType.element_beam_strain]
                     _["dependents"]     = ["beam_strain_tensor"]
                     _["shape"]          = (1,n_beams, nip_beam, 3)
                     _["convert"]        = convert.beam_strain_tensor
                     _["tracker"]        = beam_ids_tracker      
                     
                     array_requirements[ArrayType.element_shell_strain] = {}
                     _ = array_requirements[ArrayType.element_shell_strain]
                     _["dependents"]     = ["shell_strain_tensor"]
                     _["shape"]          = (1,n_shell, nip_shell, 4)
                     _["convert"]        = convert.shell_strain_tensor
                     _["tracker"]        = shell_ids_tracker                   
                
                     array_requirements[ArrayType.element_solid_strain] = {}
                     _ = array_requirements[ArrayType.element_solid_strain]
                     _["dependents"]     = ["solid_strain_tensor"]
                     _["shape"]          = (1,n_solid, nip_solid, 6)
                     _["convert"]        = convert.solid_strain_tensor
                     _["tracker"]        = solid_ids_tracker      
                     
                ISED = False 
                if ISED:
                     """
                     Output flag for strain_energy_density:
                     """
                     
                     flag = "ISED"
                     
                     database_extent_binary[flag] = {}      
                     _ = database_extent_binary[flag]
                     
                     database_extent_binary[flag][0] = []
                     database_extent_binary[flag][1] = _[0]        [    ArrayType.element_beam_strain_energy_density,
                                                                        ArrayType.element_shell_strain_energy_density,
                                                                        ArrayType.element_solid_strain_energy_density
                                                                        ]                                                
                
                     array_requirements[ArrayType.element_beam_strain_energy_density] = {}
                     _ = array_requirements[ArrayType.element_beam_strain_energy_density]
                     _["dependents"]     = ["beam_strain_energy_density"]
                     _["shape"]          = (1,n_beams)
                     _["convert"]        = convert.beam_strain_energy_density
                     _["tracker"]        = beam_ids_tracker      
                     
                     array_requirements[ArrayType.element_shell_strain_energy_density] = {}
                     _ = array_requirements[ArrayType.element_shell_strain_energy_density]
                     _["dependents"]     = ["shell_strain_energy_density"]
                     _["shape"]          = (1,n_shell)
                     _["convert"]        = convert.shell_strain_energy_density
                     _["tracker"]        = shell_ids_tracker                   
                
                     array_requirements[ArrayType.element_solid_strain_energy_density] = {}
                     _ = array_requirements[ArrayType.element_solid_strain_energy_density]
                     _["dependents"]     = ["solid_strain_energy_density"]
                     _["shape"]          = (1,n_solid)
                     _["convert"]        = convert.solid_strain_energy_density
                     _["tracker"]        = solid_ids_tracker       
                     
                NEIPH = False 
                if NEIPH:                 
                     
                    """                 
                    NEIPH Number of additional integration point history variables written to
                    the binary databases (d3plot, d3part, d3drlf) for solid elements and
                    SPH particles. The integration point data is written in the same
                    order that it is stored in memory; each material model has its own
                    history variables that are stored. For user defined materials the
                    history data that is needed for plotting should be stored before the
                    data which is not of interest. See also *DEFINE_MATERIAL_HISTORIES.
                    For output of additional integration point history
                    variables for solid elements to the elout database, see the field OPTION1
                    for *DATABASE_ELOUT.:
                    """                    
                    
                    # Assume to be unused for now
                    # Additional Radioss ouputs that are not easily converted could be added here
                    
                    None
                    
                NEIPS = False 
                if NEIPS:                
                        
                    """    
                    NEIPS Number of additional integration point history variables written to
                    the binary databases (d3plot, d3part, d3drlf) for both shell and
                    thick shell elements for each integration point; see NEIPH above
                    and *DEFINE_MATERIAL_HISTORIES. For output of additional
                    integration point history variables for shell and thick shell
                    elements to the elout database, see the fields OPTION2    
                    """
                    
                    # Assume to be unused for now
                    # Additional Radioss ouputs that are not easily converted could be added here
                    
                    None                
                    
                MAXINT = False 
                if MAXINT:                
                        
                    """                                      
                    MAXINT Number of shell and thick shell through-thickness integration
                    points for which output is written to d3plot. This does not apply
                    to the strain tensor output flagged by STRFLG.
                    MAXINT
                    (def = 3)
                    Number of
                    Integration
                    Points Description
                    3 > 3
                    (even & odd)
                    results are output for the outermost
                    (top) and innermost (bottom)
                    integration points together with
                    results for the neutral axis.
                    3 1 All three results are identical.
                    > 3 ≤ MAXINT
                    Results for the first MAXINT
                    integration points in the element will
                    be output.
                    ≠ 3 Even
                    See above. This will exclude midsurface
                    results, whereas when MAXINT
                    = 3 mid-surface results are
                    calculated and reported.
                    < 0 Any
                    MAXINT integration points are
                    output for each in plane integration
                    point location and no averaging is
                    used. This can greatly increase the
                    size of the binary databases d3plot,
                    d3thdt, and d3part.
                    See Remark 1 for more information.                 
                                     
                    1. MAXINT Field. If MAXINT is set to 3, then mid-surface, inner-surface and
                    outer-surface stresses are output at the center of the element. For an even
                    number of integration points, the points closest to the center are averaged to
                    obtain the midsurface values. If multiple integration points are used in the
                    shell plane, the stresses at the center of the element are found by computing
                    the average of these points. For MAXINT equal to 3, LS-DYNA assumes that
                    the data for the user defined integration rules are ordered from bottom to top
                    even if this is not the case. If MAXINT is not equal to 3, then the stresses at the
                    center of the element are output in the order that they are stored for the selected
                    integration rule. If multiple points are used in plane, the stresses are first
                    averaged.
                    """         
                    
                    # This will be whatever if output in Radioss
                    
                    None                 
                
                STRFLG = False 
                if STRFLG:            
                    """
                    Flag for output of strain tensors. STRFLG is interpreted digit-wise
                    STRFLG = [𝑁𝑀𝐿],
                    STRFLG = 𝐿 +𝑀 × 10 + 𝑁 × 100
                    L.EQ.1: Write strain tensor data to d3plot, elout, and dynain.
                    For shell and thick shell elements two tensors are written,
                    one at the innermost and one at the outermost integration
                    point. For solid elements a single strain
                    tensor is written.
                    M.EQ.1: Write plastic strain data to d3plot.
                    N.EQ.1: Write thermal strain data to d3plot.  
                    Examples. For STRFLG = 11 (011) LS-DYNA will write both strain
                    and plastic strain tensors, but no thermal strain tensors. Whereas
                    for STRFLG = 110, LS-DYNA will write plastic and thermal strain
                    tensors but no strain tensors. For more information and supported
                    elements and materials, see Remark 10.
                    """
                    # Skip this asume all will be created for now
                    # See ISTRA
                    None
                    
                SIGFLG = False 
                if SIGFLG:   
                    """
                    Flag for including the stress tensor for shells and solids.
                    EQ.1: include (default),
                    EQ.2: exclude for shells, include for solids.
                    EQ.3: exclude for shells and solids.             
                    """
                    
                    # Skip this asume all will be created for now
                    # See ISTRS
                    None
                    
                EPSFLG = False 
                if EPSFLG:   
                    """
                    Flag for including the effective plastic strains for shells and solids:
                    EQ.1: include (default),
                    EQ.2: exclude for shells, include for solids.
                    EQ.3: exclude for shells and solids.          
                    """
                    
                    # Skip this asume all will be created for now
                    # See ISTRS
                    None                
                    
                ENGFLG = False 
                if ENGFLG:   
                    """
                    Flag for including shell and tshell internal energy density and shell
                    thickness:
                    EQ.1: include (default),
                    EQ.2: exclude.       
                    """
                    
                    # Skip this asume all will be created for now
                    # See HYDRO
                    None                   
                    
                CMPFLG = False 
                if CMPFLG:  
                    """
                    Flag to indicate the coordinate system for output of stress and
                    strain of solids, shells and thick shells comprised of orthotropic or
                    anisotropic materials. CMPFLG affects d3plot, d3part, eloutdet,
                    and elout, with exceptions as noted below.
                    EQ.-1: Same as 1, but for *MAT_FABRIC (FORM = 14 or -14)
                    and *MAT_FABRIC_MAP the stress and strain is in engineering
                    quantities instead of Green-Lagrange strain
                    and 2nd Piola-Kirchhoff stress.
                    
                    EQ.0: global coordinate system with exception of elout for
                    shells (see EOCS in *CONTROL_OUTPUT).
                    
                    EQ.1: local material coordinate system (as defined by AOPT
                    and associated parameters in the *MAT input, and if applicable,
                    by angles B1, B2, etc. in *SECTION_SHELL,
                    *SECTION_TSHELL, or *PART_COMPOSITE, and by
                    optional input in the *ELEMENT data). CMPFLG = 1
                    affects both d3plot and elout databases with exception
                    that EOCS = 1 or 2 overrides CMPFLG for shell output in
                    elout.                
                    """
                    
                    # This might resolve the discepency between Radioss and LS-Dyna
                    # Check - Does Radioss output in local or global coordinate system?
                    # If local set header flag to indicate EQ.1 (2 Dyna runs to find the flag & value)
                    None          
                    
                IEVERP = False 
                if IEVERP:  
                    """
                    Every output state for the d3plot database is written to a separate
                    file.
                    EQ.0: more than one state can be on each plot file,
                    EQ.1: one state only on each plot file.               
                    """
                    
                    # Assume EQ.1  
                    None                 
                
                BEAMIP = False 
                if BEAMIP: 
                            
                    """
                    Number of beam integration points for output. This option does
                    not apply to beams that use a resultant formulation. See Remark
                    2.            
                    
                    
                    Remark 2. BEAMIP Field. Beam stresses are output if and only if BEAMIP is greater than
                    zero. In this latter case the data that is output is written in the same order that
                    the integration points are defined. The data at each integration point consists
                    of the following five values for elastic-plastic Hughes-Liu beams: the normal
                    stress, 𝜎𝑟𝑟; the transverse shear stresses, 𝜎𝑟𝑠 and 𝜎𝑡𝑟; the effective plastic strain;
                    and the axial strain which is logarithmic. For beams that are not elastic-plastic,
                    the first history variable, if any, is output instead of the plastic strain. For the
                    beam elements of Belytschko and his co-workers, the transverse shear stress
                    components are not used in the formulation. No data is output for the Belytschko-
                    Schwer resultant beam.
                    """
                    
                    # This will be whatever Radioss ouputs
                    # It might be necessary to set the header flag manually
                    
                    None
                    
                DCOMP = False 
                if DCOMP:                 
                    """    
                    DCOMP Data compression to eliminate rigid body data:
                    EQ.1: off (default), no rigid body data compression,
                    EQ.2: on, rigid body data compression active,
                    EQ.3: off, no rigid body data compression, but all nodal
                    velocities and accelerations are eliminated from the database.
                    EQ.4: on, rigid body data compression active and all nodal
                    velocities and accelerations are eliminated from the database.
                    EQ.5: on, rigid body data compression active and rigid nodal
                    data are eliminated from the database. Only 6 DOF rigid
                    body motion is written.
                    EQ.6: on, rigid body data compression active, rigid nodal data, 
                    and all nodal velocities and accelerations are eliminated
                    from the database. Only 6 DOF rigid body motion is written.               
                    """
                    
                    # Assume no compression
                    
                    None
                
                SHGE = False
                if SHGE:
                    """
                    Flag for including shell hourglass energy density:
                    EQ.1: off (default), no hourglass energy written,
                    EQ.2: on.
                    """
                     
                    flag = "SHGE"
                    
                    database_extent_binary[flag] = {}      
                    _ = database_extent_binary[flag]
                    
                    database_extent_binary[flag][1] = []
                    database_extent_binary[flag][2] = _[1] +      [ArrayType.element_shell_hourglass_energy_density
                                                                   ]
                                              
                    # Check the units here I *think* Dyna is energy per unit volume whereas Radioss is energy per unit mass
                    array_requirements[ArrayType.element_shell_hourglass_energy_density] = {}
                    _ = array_requirements[ArrayType.element_shell_hourglass_energy_density]
                    _["dependents"]     = ["element_shell_hourglass_energy_density", "material density"]
                    _["shape"]          = (1,n_shell)
                    _["convert"]        = convert.element_shell_hourglass_energy_density
                    _["tracker"]        = shell_ids_tracker  
                    
                STSSZ = False 
                if STSSZ:
                    """
                    Flag for including shell element time step, mass, or added mass:
                    EQ.1: off (default),
                    EQ.2: output time step size,
                    EQ.3: output mass, added mass, or time step size.
                    Mass Scaling. If mass scaling is active, the output of the time step size reveals
                    little information about the calculation. If global mass scaling is used for a
                    constant time step, the total element mass is output; however, if the mass is
                    increased so that a minimum time step size is maintained (DT2MS is negative),
                    the added mass is output. Also, see the control card *CONTROL_TIMESTEP.
                    """
                     
                    flag = "STSSZ"
                    
                    database_extent_binary[flag] = {}      
                    _ = database_extent_binary[flag]
                    
                    database_extent_binary[flag][1] = []
                    database_extent_binary[flag][2] = _[1] + [ ArrayType.element_shell_timestep_size
                                                                    ]            
                    database_extent_binary[flag][3] = _[2] + [ ArrayType.element_shell_mass, 
                                                                    ArrayType.element_shell_added_mass]
                                                                       
                                              
                    array_requirements[ArrayType.element_shell_timestep_size] = {}
                    _ = array_requirements[ArrayType.element_shell_timestep_size]
                    _["dependents"]     = ["element_shell_hourglass_timestep_size"]
                    _["shape"]          = (1,n_shell)
                    _["convert"]        = convert.element_shell_hourglass_timestep_size
                    _["tracker"]        = shell_ids_tracker    
                    
                    array_requirements[ArrayType.element_shell_mass] = {}
                    _ = array_requirements[ArrayType.element_shell_mass]
                    _["dependents"]     = ["element_shell_mass"]
                    _["shape"]          = (1,n_shell)
                    _["convert"]        = convert.element_shell_mass
                    _["tracker"]        = shell_ids_tracker       
                    
                    array_requirements[ArrayType.element_shell_added_mass] = {}
                    _ = array_requirements[ArrayType.element_shell_added_mass]
                    _["dependents"]     = ["element_shell_added_mass"]
                    _["shape"]          = (1,n_shell)
                    _["convert"]        = convert.element_shell_added_mass
                    _["tracker"]        = shell_ids_tracker                  
               
                N3THDT = False 
                if N3THDT:                 
                    """    
                    Flag for including material energy in d3thdt database:
                    EQ.1: off, energy is not written to d3thdt database,
                    EQ.2: on (default), energy is written to d3thdt database.              
                    """
                    
                    # Assume no compression
                    
                    None     
                    
                    
                IALEMAT = False 
                if IALEMAT:
                    """
                    Output solid part ID list containing ALE materials.
                    EQ.1: on (default)
                    """
                     
                    # This will need to be defined if ALE is present
                    
                    flag = "IALEMAT"
                    
                    database_extent_binary[flag] = {}      
                    _ = database_extent_binary[flag]
                            
                    database_extent_binary[flag][1] = [ArrayType.element_solid_ale_parts]
      
                    array_requirements[ArrayType.element_solid_ale_parts] = {}
                    _ = array_requirements[ArrayType.element_solid_ale_parts]
                    _["dependents"]     = ["element_solid_ale_parts"]
                    _["shape"]          = (1,n_ale_solid_parts)
                    _["convert"]        = convert.element_solid_ale_parts
                    _["tracker"]        = ale_solid_ids_tracker
                    
    
                NINTSLD = False 
                if NINTSLD:                 
                    """    
                    Number of solid element integration points written to the LS-DYNA
                    database. When NINTSLD is set to 1 (default) or to any value
                    other than 8, integration point values are averaged and only those
                    averages are written output. To obtain values for individual
                    integration points, set NINTSLD to 8, even if the multi-integration
                    point solid has fewer than 8 integration points.             
                    """
                    
                    # This will be whatever is output in Radioss
                    
                    None           
                     
                PKP_SEN = False 
                if PKP_SEN:
                    """
                    Flag to output the peak pressure and surface energy (energy per
                    unit area) computed by each contact interface into the interface
                    force database. To obtain the surface energy, FRCENG, must be
                    sent to 1 on the control contact card. When PKP_SEN = 1, it is
                    possible to identify the energies generated on the upper and lower
                    shell surfaces, which is important in metal forming applications.
                    This data is mapped after each h-adaptive remeshing.
                    EQ.0: No data is written
                    EQ.1: Output the peak pressures and surface energy by contact
                    interface
                    """
                    
                    # Double check units, should just be pressure and energy 
                     
                    flag = "PKP_SEN"
                    
                    database_extent_binary[flag] = {}      
                    _ = database_extent_binary[flag]
                    
                    database_extent_binary[flag][0] = []
                    database_extent_binary[flag][1] = _[0] + [ ArrayType.contact_peak_pressures,
                                                                    ArrayType.contact_surface_energy
                                                                    ]            
                                          
                    array_requirements[ArrayType.contact_peak_pressures] = {}
                    _ = array_requirements[ArrayType.contact_peak_pressures]
                    _["dependents"]     = ["contact_peak_pressures"]
                    _["shape"]          = (1,n_contact)
                    _["convert"]        = convert.element_shell_contact_peak_pressures
                    _["tracker"]        = contact_ids_tracker
                    
                    array_requirements[ArrayType.contact_surface_energy] = {}
                    _ = array_requirements[ArrayType.contact_surface_energy]
                    _["dependents"]     = ["contact_surface_energy"]
                    _["shape"]          = (1,n_contact)
                    _["convert"]        = convert.element_shell_contact_surface_energy
                    _["tracker"]        = contact_ids_tracker 
                    
                SCLP = False 
                if SCLP:                 
                    """    
                    A scaling parameter used in the computation of the peak pressure.
                    This parameter is generally set to unity (the default), but it must be
                    greater than 0.            
                    """
                    
                    # This will be whatever is output in Radioss
                    
                    None         
                    
                HYDRO = False 
                if HYDRO:
                    """
                    Either 3, 5 or 7 additional history variables useful to shock physics
                    are output as the last history variables to d3plot (does not apply to
                    elout). For HYDRO = 1, the internal energy per reference volume,
                    the reference volume, and the pressure from bulk viscosity are
                    added to the database; for HYDRO = 2, the relative volume and
                    current density are also added; and for HYDRO = 4, volumetric
                    strain (defined as relative volume – 1.0) and hourglass energy per
                    unit initial volume are additionally added. These history variables
                    are not valid for ALE elements.
                    interface
                    """
                    
                    # Double check units, should just be pressure and energy 
                     
                    flag = "HYDRO"
                    
                    database_extent_binary[flag] = {}      
                    _ = database_extent_binary[flag]
                    
                    database_extent_binary[flag][0] = []
                    database_extent_binary[flag][1] = _[0] +      [ArrayType.element_shell_internal_energy_density,
                                                                   ArrayType.element_shell_volume,
                                                                   ArrayType.element_shell_bulk_viscocity_pressure,
                                                                   ArrayType.element_solid_internal_energy_density,
                                                                   ArrayType.element_solid_volume,
                                                                   ArrayType.element_solid_bulk_viscocity_pressure                                                                                                                              
                                                                   ]
                    database_extent_binary[flag][2] = _[1] +      [ArrayType.element_shell_relative_volume,
                                                                   ArrayType.element_shell_density,
                                                                   ArrayType.element_solid_relative_volume,
                                                                   ArrayType.element_solid_density                                                                                                                            
                                                                   ]  
    
                    database_extent_binary[flag][4] = _[2] +      [ArrayType.element_shell_volumetric_strain,
                                                                   ArrayType.element_shell_hourglass_energy_per_initial_volume,
                                                                   ArrayType.element_solid_volumetric_strain,
                                                                   ArrayType.element_solid_hourglass_energy_per_initial_volume,                                                                                                                           
                                                                   ]                
                    
                                          
                    array_requirements[ArrayType.element_shell_internal_energy_density] = {}
                    _ = array_requirements[ArrayType.element_shell_internal_energy_density]
                    _["dependents"]     = ["shell_internal_energy_density", "shell_density"]
                    _["shape"]          = (1,n_shell)
                    _["convert"]        = convert.element_shell_internal_energy_density
                    _["tracker"]        = shell_ids_tracker    
                    
                    array_requirements[ArrayType.element_shell_volume] = {}
                    _ = array_requirements[ArrayType.element_shell_volume]
                    _["dependents"]     = ["shell_volume"]
                    _["shape"]          = (1,n_shell)
                    _["convert"]        = convert.element_shell_volume
                    _["tracker"]        = shell_ids_tracker                  
                
                    array_requirements[ArrayType.element_shell_bulk_viscocity_pressure] = {}
                    _ = array_requirements[ArrayType.element_shell_bulk_viscocity_pressure]
                    _["dependents"]     = ["shell_viscocity_pressure"]
                    _["shape"]          = (1,n_shell)
                    _["convert"]        = convert.element_viscocity_pressure
                    _["tracker"]        = shell_ids_tracker       
                    
                    array_requirements[ArrayType.element_shell_relative_volume] = {}
                    _ = array_requirements[ArrayType.element_shell_bulk_relative_volume]
                    _["dependents"]     = ["shell_relative_volume"]
                    _["shape"]          = (1,n_shell)
                    _["convert"]        = convert.element_shell_bulk_relative_volume
                    _["tracker"]        = shell_ids_tracker       
                    
                    array_requirements[ArrayType.element_shell_density] = {}
                    _ = array_requirements[ArrayType.element_shell_density]
                    _["dependents"]     = ["shell_density"]
                    _["shape"]          = (1,n_shell)
                    _["convert"]        = convert.element_shell_volume
                    _["tracker"]        = shell_ids_tracker                     
                                
                    array_requirements[ArrayType.element_shell_volumetric_strain] = {}
                    _ = array_requirements[ArrayType.element_shell_volumetric_strain]
                    _["dependents"]     = ["shell_volumetric_strain"]
                    _["shape"]          = (1,n_shell)
                    _["convert"]        = convert.element_shell_volumetric_strain
                    _["tracker"]        = shell_ids_tracker    
                    
                    array_requirements[ArrayType.element_shell_hourglass_energy_per_initial_volume] = {}
                    _ = array_requirements[ArrayType.element_shell_hourglass_energy_per_initial_volume]
                    _["dependents"]     = ["shell_hourglass_energy_per_initial_volume"]
                    _["shape"]          = (1,n_shell)
                    _["convert"]        = convert.element_shell_hourglass_energy_per_initial_volume
                    _["tracker"]        = shell_ids_tracker    
                    
                    array_requirements[ArrayType.element_solid_internal_energy_density] = {}
                    _ = array_requirements[ArrayType.element_solid_internal_energy_density]
                    _["dependents"]     = ["solid_internal_energy_density", "solid_density"]
                    _["shape"]          = (1,n_solid)
                    _["convert"]        = convert.element_solid_internal_energy_density
                    _["tracker"]        = solid_ids_tracker    
                    
                    array_requirements[ArrayType.element_solid_volume] = {}
                    _ = array_requirements[ArrayType.element_solid_volumey]
                    _["dependents"]     = ["solid_volume"]
                    _["shape"]          = (1,n_solid)
                    _["convert"]        = convert.element_solid_volume
                    _["tracker"]        = solid_ids_tracker                  
                
                    array_requirements[ArrayType.element_solid_bulk_viscocity_pressure] = {}
                    _ = array_requirements[ArrayType.element_solid_bulk_viscocity_pressure]
                    _["dependents"]     = ["solid_viscocity_pressure"]
                    _["shape"]          = (1,n_solid)
                    _["convert"]        = convert.element_viscocity_pressure
                    _["tracker"]        = solid_ids_tracker       
                    
                    array_requirements[ArrayType.element_solid_relative_volume] = {}
                    _ = array_requirements[ArrayType.element_solid_bulk_relative_volume]
                    _["dependents"]     = ["solid_relative_volume"]
                    _["shape"]          = (1,n_solid)
                    _["convert"]        = convert.element_solid_bulk_relative_volume
                    _["tracker"]        = solid_ids_tracker       
                    
                    array_requirements[ArrayType.element_solid_density] = {}
                    _ = array_requirements[ArrayType.element_solid_density]
                    _["dependents"]     = ["solid_density"]
                    _["shape"]          = (1,n_solid)
                    _["convert"]        = convert.element_solid_volume
                    _["tracker"]        = solid_ids_tracker                     
                                
                    array_requirements[ArrayType.element_solid_volumetric_strain] = {}
                    _ = array_requirements[ArrayType.element_solid_volumetric_strain]
                    _["dependents"]     = ["solid_volumetric_strain"]
                    _["shape"]          = (1,n_solid)
                    _["convert"]        = convert.element_solid_volumetric_strain
                    _["tracker"]        = solid_ids_tracker    
                    
                    array_requirements[ArrayType.element_solid_hourglass_energy_per_initial_volume] = {}
                    _ = array_requirements[ArrayType.element_solid_hourglass_energy_per_initial_volume]
                    _["dependents"]     = ["solid_hourglass_energy_per_initial_volume"]
                    _["shape"]          = (1,n_solid)
                    _["convert"]        = convert.element_solid_hourglass_energy_per_initial_volume
                    _["tracker"]        = solid_ids_tracker                 
                
                MSSCL = False 
                if MSSCL:
                    """
                    Output nodal information related to mass scaling into the d3plot
                    database. This option can be activated if and only if DT2MS < 0.0.
                    See control keyword *CONTROL_TIMESTEP.
                    EQ.0: No data is written
                    EQ.1: Output incremental nodal mass
                    EQ.2: Output percentage increase in nodal mass
                    """
                     
                    flag = "MSSCL"
                    
                    database_extent_binary[flag] = {}      
                    _ = database_extent_binary[flag]
                    
                    database_extent_binary[flag][0] = []
                    database_extent_binary[flag][1] = _[0] +      [ArrayType.nodal_mass
                                                                   ]     
                    database_extent_binary[flag][2] = _[1] +      [ArrayType.percentage_increase_nodal_mass
                                                                   ] 
                              
                
                    array_requirements[ArrayType.nodal_mass] = {}
                    _ = array_requirements[ArrayType.nodal_mass]
                    _["dependents"]     = ["nodal_mass"]
                    _["shape"]          = (1,n_nodes)
                    _["convert"]        = convert.nodal_mass
                    _["tracker"]        = nodes_ids_tracker           
                    
                    array_requirements[ArrayType.percentage_increase_nodal_mass] = {}
                    _ = array_requirements[ArrayType.percentage_increase_nodal_mass]
                    _["dependents"]     = ["percentage_increase_nodal_mass"]
                    _["shape"]          = (1,n_nodes)
                    _["convert"]        = convert.percentage_increase_nodal_mass
                    _["tracker"]        = nodes_ids_tracker      
                
                THERM = False 
                if THERM:
                    """
                    Output of thermal data to d3plot. The use of this option
                    (THERM > 0) may make the database incompatible with other 3rd
                    party software.
                    EQ.0: (default) output temperature
                    EQ.1: output temperature
                    EQ.2: output temperature and flux
                    EQ.3: output temperature, flux, and shell bottom and top
                    surface temperature
                    """
                    
                    # Double check units, should just be pressure and energy 
                     
                    flag = "THERM"
                    
                    database_extent_binary[flag] = {}      
                    _ = database_extent_binary[flag]
                    
                    database_extent_binary[flag][0] = []
                    database_extent_binary[flag][1] = _[0] +      [ ArrayType.element_shell_temperature,
                                                                    ArrayType.element_solid_temperature
                                                                    ]          
                    database_extent_binary[flag][2] = _[1] +       [ArrayType.element_shell_flux,
                                                                    ArrayType.element_solid_flux
                                                                    ]
                    database_extent_binary[flag][3] = _[2] +      [ ArrayType.element_shell_top_bottom_temperature,
                                                                    ArrayType.element_solid_top_bottom_temperature
                                                                    ]
                                                                  
                    array_requirements[ArrayType.element_shell_temperature] = {}
                    _ = array_requirements[ArrayType.element_shell_temperature]
                    _["dependents"]     = ["element_shell_temperature"]
                    _["shape"]          = (1,n_shell)
                    _["convert"]        = convert.element_shell_temperature
                    _["tracker"]        = shell_ids_tracker    
                    
                    array_requirements[ArrayType.element_shell_flux] = {}
                    _ = array_requirements[ArrayType.element_shell_flux]
                    _["dependents"]     = ["element_shell_flux"]
                    _["shape"]          = (1,n_shell)
                    _["convert"]        = convert.element_shell_flux
                    _["tracker"]        = shell_ids_tracker  
                    
                    array_requirements[ArrayType.element_shell_top_bottom_temperature] = {}
                    _ = array_requirements[ArrayType.element_shell_top_bottom_temperature]
                    _["dependents"]     = ["element_shell_top_bottom_temperature"]
                    _["shape"]          = (1,2, n_shell)
                    _["convert"]        = convert.shell_top_bottom_temperature
                    _["tracker"]        = shell_ids_tracker          
                    
                    array_requirements[ArrayType.element_solid_temperature] = {}
                    _ = array_requirements[ArrayType.element_solid_temperature]
                    _["dependents"]     = ["element_solid_temperature"]
                    _["shape"]          = (1,n_solid)
                    _["convert"]        = convert.element_solid_temperature
                    _["tracker"]        = solid_ids_tracker    
                    
                    array_requirements[ArrayType.element_solid_flux] = {}
                    _ = array_requirements[ArrayType.element_solid_flux]
                    _["dependents"]     = ["element_solid_flux"]
                    _["shape"]          = (1,n_solid)
                    _["convert"]        = convert.element_solid_flux
                    _["tracker"]        = solid_ids_tracker  
                    
                    array_requirements[ArrayType.element_solid_top_bottom_temperature] = {}
                    _ = array_requirements[ArrayType.element_solid_top_bottom_temperature]
                    _["dependents"]     = ["element_solid_top_bottom_temperature"]
                    _["shape"]          = (1,2, n_solid)
                    _["convert"]        = convert.solid_top_bottom_temperature
                    _["tracker"]        = solid_ids_tracker                          
                
                
                INTOUT = False 
                if INTOUT:
                    """
                    Output stress/strain at all integration points for detailed element
                    output in the ASCII file eloutdet. DT and BINARY of *DATABASE_
                    ELOUT apply to eloutdet. See Remark 4.
                    EQ.STRESS: when stress output is required
                    EQ.STRAIN: when strain output is required
                    EQ.ALL: when both stress and strain output are required
                    """            
                    
                    # For Binout - Not relevent for D3plot 
                    
                NODOUT = False
                if NODOUT:
                    """
                    Output extrapolated stress/strain at connectivity nodes for
                    detailed element output in the ASCII file eloutdet. DT and BINARY
                    of *DATABASE_ELOUT apply to eloutdet. See Remark 4.
                    EQ.STRESS: when stress output is required
                    EQ.STRAIN: when strain output is required
                    EQ.ALL: when both stress and strain output are
                    required
                    EQ.STRESS_GL: when nodal averaged stress output along the
                    global coordinate system is required
                    EQ.STRAIN_GL: when nodal averaged strain output along the
                    global coordinate system is required
                    EQ.ALL_GL: for global nodal averaged stress and strain
                    output
                    """            
                    
                    # For Binout - Not relevent for D3plot     
                    
                RESPLT = False
                if RESPLT:
                    """
                    Output of translational and rotational residual forces to d3plot and
                    d3iter.
                    EQ.0: No output
                    EQ.1: Output residual
                    """            
                    
                    # For Inplicit as a convergence check - unclear what this involves               
                
                NEIPB = False
                if NEIPB:
                    """
                    Number of additional element or integration point history
                    variables written to the binary databases (d3plot, d3part, d3drlf)
                    for beam elements; see NEIPH above, BEAMIP and *DEFINE_MATERIAL_
                    HISTORIES. See also Remark 12. For output of
                    additional integration point history variables for beam elements to
                    the elout database, see the variable OPTION4 in *DATABASE_-
                    ELOUT.
                    
                    Remark 12. History Variables for Beams (NEIPB). In general, NEIPB follows the same
                    conventions as NEIPH and NEIPS do for solid and shell elements and is supported
                    in LS-PrePost v4.3 or later. Average, min and max values for each element
                    are output, including data for resultant elements. If BEAMIP is nonzero,
                    then element data is complemented with BEAMIP integration point values that
                    can be examined individually. Beam history data is post-processed similarly to
                    that of solid and shell element history data.                
                    """            
                    
                    # This will be whatever Radioss outputs       
                    
                QUADSLD = False
                if QUADSLD:
                    """
                    (Under development) Output flag for quadratic solid types 24, 25,
                    and 26:
                    EQ.0: average stress and strain and rendered with 8 nodes
                    EQ.1: average stress and strain and rendered with all edge and
                    face nodes
                    EQ.2: all integration points written and rendered with all edge
                    and face nodes               
                    """            
                    
                    # Leave this for now   
                    
                CUBSLD = False 
                if CUBSLD:
                    """
                    (Under development) Output flag for cubic solids types 27, 28, and
                    29:
                    EQ.0: average stress and strain and rendered with 8 nodes
                    EQ.1: average stress and strain and rendered with all edge and
                    face nodes
                    EQ.2: all integration              
                    """            
                    
                    # Leave this for now                 
            
#########################################
# New part
#########################################
            if True:
                
                database_extent_binary ={}
                array_requirements = {}

                flag = "NODES"
                
                if rr.raw_header["nbNodes"] > 0:   
                    
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
                    _["shape"]          = (1,n_nodes)
                    _["convert"]        = None
                    _["tracker"]        = None
                    
                    # Dyna output
                    array_requirements[ArrayType.node_acceleration] = {}
                    _ = array_requirements[ArrayType.node_acceleration]
                    # Radioss outputs needed to compute Dyna output
                    _["dependents"]     = ["node_acceleration"]
                    _["shape"]          = (1,n_nodes)
                    _["convert"]        = None
                    _["tracker"]        = None




                flag = "SHELLS"
                
                if rr.raw_header["nbEFunc"] > 0:   
                    
                    database_extent_binary[flag] = {}      
                    _ = database_extent_binary[flag]
                    
                    database_extent_binary[flag][0] = []
                    database_extent_binary[flag][1] = _[0] + [ArrayType.element_shell_is_alive]  
                    database_extent_binary[flag][2] = _[1] + [ArrayType.element_shell_thickness]        
                    database_extent_binary[flag][3] = _[2] + [ArrayType.element_shell_internal_energy]    
                    #database_extent_binary[flag][4] = _[3] + [ArrayType.element_shell_stress]    
                    
                    # Dyna output
                    array_requirements[ArrayType.element_shell_thickness] = {}
                    _ = array_requirements[ArrayType.element_shell_thickness]
                    # Radioss outputs needed to compute Dyna output
                    _["dependents"]     = ["element_shell_thickness"]
                    _["shape"]          = (1,n_shell)
                    _["convert"]        = None
                    _["tracker"]        = shell_ids_tracker
                    
                    # Dyna output
                    array_requirements[ArrayType.element_shell_is_alive] = {}
                    _ = array_requirements[ArrayType.element_shell_is_alive]
                    # Radioss outputs needed to compute Dyna output
                    _["dependents"]     = ["element_shell_is_alive"]
                    _["shape"]          = (1,n_shell)
                    _["convert"]        = None
                    _["tracker"]        = shell_ids_tracker
                    
                    # Dyna output
                    array_requirements[ArrayType.element_shell_internal_energy] = {}
                    _ = array_requirements[ArrayType.element_shell_internal_energy]
                    # Radioss outputs needed to compute Dyna output
                    _["dependents"]     = ["element_shell_specific_energy"]
                    _["shape"]          = (1,n_shell)
                    _["convert"]        = None
                    _["tracker"]        = shell_ids_tracker
                    
                   # Dyna output
                    #array_requirements[ArrayType.element_shell_stress] = {}
                    #_ = array_requirements[ArrayType.element_shell_stress]
                    ## Radioss outputs needed to comptute Dyna output
                    ##_["dependents"]     = ["element_shell_stress_(upper)","element_shell_stress_(lower)"]
                    #_["dependents"]     = ["Stress (upper)","Stress (lower)"]
                    #_["shape"]          = (1,n_shell, 1, 4)
                    #_["convert"]        = convert.element_shell_stress
                    #_["tracker"]        = shell_ids_tracker
#########################################

            
            element_beam_is_alive           = []
            element_shell_is_alive          = []
            element_solid_is_alive          = []
            
            timesteps                       = []                        

            "Beam Updates"    
            
            if rr.raw_header["nbElts1D"] > 0 and "element_beam_is_alive" in rr.arrays:
                element_beam_is_alive.append(readAndConvert.apply_sorter(rr.arrays["element_beam_is_alive"], beam_ids_tracker).astype("<f"))
            
            
            "shell updates"

            if rr.raw_header["nbFacets"] > 0 and "element_shell_is_alive" in rr.arrays:
                element_shell_is_alive.append(readAndConvert.apply_sorter(rr.arrays["element_shell_is_alive"], shell_ids_tracker).astype("<f"))

            "Solid Updates"
            
            if rr.raw_header["nbElts3D"] > 0 and "element_solid_is_alive" in rr.arrays:
                element_solid_is_alive.append(readAndConvert.apply_sorter(rr.arrays["element_solid_is_alive"], solid_ids_tracker).astype("<f"))

            "Timestep Updates"
            timesteps.append(rr.arrays["timesteps"])
        
            "Assign the arrays to the D3PLOT class for writing"
                        
            "Generate the availability check"
            
            dependency_check = {}
            flag_max        =  {}

            #print(rr.arrays)

            for flag in database_extent_binary:
                flag_max[flag] = -float('inf')
                _ = database_extent_binary[flag]
                for array_group in _:
                    __ = _[array_group]
                    for array_output in __:
                        # Check all dependencies exist
                        all_dependents_exist = True
                        ___ = array_requirements[array_output]["dependents"]
                        #print(___)
                        for check_dependent in ___:
                            if check_dependent not in rr.arrays:
                                all_dependents_exist = False
                                break
                    
                        dependency_check[array_output] = all_dependents_exist
                        if all_dependents_exist:
                            flag_max[flag] = max(flag_max[flag],array_group) 
            
            
            
            "Shells"    
            if rr.raw_header["nbFacets"] > 0:
                self._d3plot.arrays[ArrayType.element_shell_ids]              = element_shell_ids_out
                self._d3plot.arrays[ArrayType.element_shell_node_indexes]     = readAndConvert.apply_sorter(rr.arrays["element_shell_node_indexes"], shell_ids_tracker)
                self._d3plot.arrays[ArrayType.element_shell_part_indexes]     = readAndConvert.apply_sorter(rr.arrays["element_shell_part_indexes"], shell_ids_tracker)
                self._d3plot.arrays[ArrayType.element_shell_is_alive]         = np.array(element_shell_is_alive)
                
            
            shell_part_ids                                          =  np.array(rr.raw_arrays["pTextA"]).astype("U9").astype(int)
            shell_part_num = len(shell_part_ids)
            
            "Beams"
            
            if rr.raw_header["nbElts1D"] > 0: 
                additional_beam_number                                      = rr.raw_header["nbElts1D"] - len(rr.arrays["element_beam_part_indexes"])
                additional_beams                                            = np.zeros(shape=(additional_beam_number))
                additional_beams.fill(max(rr.arrays["element_beam_part_indexes"])+1)
                element_beam_part_indexes                                   = np.concatenate([rr.arrays["element_beam_part_indexes"], additional_beams])
                self._d3plot.arrays[ArrayType.element_beam_ids]             = element_beam_ids_out
                element_beam_node_indexes=np.zeros(shape=(rr.raw_header["nbElts1D"],5))
                element_beam_node_indexes[:,0]                              = rr.arrays["element_beam_node_indexes"][:,0]
                element_beam_node_indexes[:,1]                              = rr.arrays["element_beam_node_indexes"][:,1]
                self._d3plot.arrays[ArrayType.element_beam_node_indexes]    = readAndConvert.apply_sorter(element_beam_node_indexes, beam_ids_tracker).astype(int)
                self._d3plot.arrays[ArrayType.element_beam_part_indexes]    = readAndConvert.apply_sorter(element_beam_part_indexes, beam_ids_tracker).astype(int) + shell_part_num
                self._d3plot.arrays[ArrayType.element_beam_is_alive]        = np.array(element_beam_is_alive)
            
                "Partless beams exists so we have to generate a dummy one 1000000000 - this could be improved"
            
                beam_part_ids                                       =   np.concatenate([np.array(rr.raw_arrays["pText1DA"])\
                                                                                        .astype("U9").astype(int),\
                                                                                        np.array([1000000000])])
                beam_part_num                                       =   len(beam_part_ids)            
            
            else:
                beam_part_ids                                       =   []
                beam_part_num                                       =   0            
            
            "Solids"
            
            if rr.raw_header["nbElts3D"] > 0: 
                                
                element_solid_node_indexes = rr.arrays["element_solid_node_indexes"]
                
                for i_solid, solid in enumerate(rr.arrays["element_solid_node_indexes"]):
                    # Correct for tetra node connectivity
                    _ = len(set(solid))
                    if _ == 4:
                        element_solid_node_indexes[i_solid] = np.array([solid[0], solid[2], solid[4], solid[5], solid[5], solid[5], solid[5], solid[5]])
                
                self._d3plot.arrays[ArrayType.element_solid_node_indexes] = readAndConvert.apply_sorter(element_solid_node_indexes, solid_ids_tracker)
                self._d3plot.arrays[ArrayType.element_solid_part_indexes] = readAndConvert.apply_sorter(rr.arrays["element_solid_part_indexes"], solid_ids_tracker) + shell_part_num + beam_part_num
            
                self._d3plot.arrays[ArrayType.element_solid_is_alive]       = np.array(element_solid_is_alive)
                self._d3plot.arrays[ArrayType.element_solid_ids]            = element_solid_ids_out            
            
                solid_part_ids                                       =  np.array(rr.raw_arrays["pText3DA"]).astype("U9").astype(int)
            
            else:
                solid_part_ids                                      =   [] 
                
            "Parts"
    
            self._d3plot.arrays[ArrayType.part_ids]                     = np.concatenate([shell_part_ids, beam_part_ids, solid_part_ids])
            self._d3plot.arrays[ArrayType.part_mass]                    = np.array([np.zeros(shape=(np.concatenate([shell_part_ids, beam_part_ids, solid_part_ids]).shape))])
            self._d3plot.arrays[ArrayType.global_timesteps]             = np.array(timesteps)
                                                           
            "Zero Arrays - For compatibility"  
            
            self._d3plot.arrays[ArrayType.element_shell_stress]                         = np.zeros(shape=(1,rr.raw_header["nbFacets"], 2, 6))
            self._d3plot.arrays[ArrayType.element_shell_effective_plastic_strain]       = np.zeros(shape=(1, rr.raw_header["nbFacets"], 2))
          
            
            "Generate the output arrays"
            
            dependency_check["IACC"] = False            
            for flag in database_extent_binary:
                _ = database_extent_binary[flag]
                for array_group in _:
                    __ = _[array_group]
                    
                    if array_group <= flag_max[flag]:
                        for array_output in __:
                            #print(array_output)
                            #print(dependency_check[array_output])
                            if dependency_check[array_output]:
                                conversion_function = array_requirements[array_output]["convert"]
                                tracker             = array_requirements[array_output]["tracker"]
                                if conversion_function:
                                    conversion_inputs = [rr.arrays[i] for i in array_requirements[array_output]["dependents"]]
                                    if tracker is not None:
                                        self._d3plot.arrays[array_output] = readAndConvert.apply_sorter(conversion_function(*conversion_inputs), tracker)[np.newaxis, :]
                                    else:
                                        self._d3plot.arrays[array_output] = conversion_function(*conversion_inputs)[np.newaxis, :]
                                else:
                                    if tracker is not None:
                                        self._d3plot.arrays[array_output] = readAndConvert.apply_sorter(rr.arrays[array_requirements[array_output]["dependents"][0]], tracker)[np.newaxis, :]
                                    else:
                                        self._d3plot.arrays[array_output] = rr.arrays[array_requirements[array_output]["dependents"][0]][np.newaxis, :]
                                        
            "Write out the d3plot state"
                                        
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
             
    file_stem = "path/model"
    a2d = readAndConvert(file_stem)
