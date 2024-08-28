#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

from typing import Union

import numpy as np
import rich

from lasso.io.binary_buffer import BinaryBuffer
from lasso.logging import get_logger

from .RadiossReader import RadiossReader

LOGGER = get_logger(__file__)

class RadiossTHReader:
    ''' Class for reading a Radioss TH File
    '''

    # meta
    filepath: str = ""  #: filepath of file

    # file info
    itype: np.dtype = np.dtype("int32").newbyteorder('>') #: integer type of Radioss Animation File
    itype_s: np.dtype = np.dtype("int16").newbyteorder('>') #: short integer type of Radioss Animation File
    ftype: np.dtype = np.dtype("float32").newbyteorder('>') #: floating point type of Radioss Animation File
    
    wordsize: int = 4  #: size of words in bytes (4 = single precision, 8 = double precision)
    wordsize_s: int = 2  

    def __init__(self, filepath: Union[str, BinaryBuffer, None] = None):
        """ Create a RadiossReader instance

        Parameters
        ----------
        filepath: Union[str, BinaryBuffer, None]
            path to a Radioss Animation File file or a buffer holding Radioss Animation File memory

        Returns
        -------
        RadiossReader: RadiossReader
            radioss state arrays

        Examples
        --------

            print(rr.raw_header.keys())
            
            print(rr.raw_arrays.keys())
            
            print(rr.arrays.keys())

        Notes
        -----
            Reads one Radioss Anim state into memory and produces:
                Header      - Raw header data
                Raw_Array   - Packed array data
                Array       - Unpacked data aims to be similar to LS-Dyna
        """

        if filepath is not None:
            self.load_file(filepath)

    def print(self) -> None:
        """ Print the header
        """
        rich.print(self.__dict__)

    def _read_file_buffer(self, filepath: str) -> BinaryBuffer:
        """ Reads a Radioss Animation Files header

        Parameters
        ----------
        filepath: str
            path to Radioss Animation File

        Returns
        -------
        bb: BinaryBuffer
            buffer holding the exact header data in binary form
        """

        LOGGER.debug("_read_file_buffer start")
        LOGGER.debug(f"filepath: {filepath}")

        bb = BinaryBuffer(filepath)

        LOGGER.debug("_read_file_buffer end")

        return bb
    
    def add_array(self, length: int, width: int, position: int, bb: BinaryBuffer, dtype: (np.dtype, str)):
        """ Load a numpy array into dictionary & offset the position index
        
            Returns
            -------
            array_from_buffer: np.array
                The new numpy array
                
            position: int
                The updated word position
            
            """
        if dtype != str:
            array_from_buffer=[]*width
            if length:
                # short and normal floats and ints
                if isinstance(dtype, np.dtype):
                    wordlength = length*width*self.wordsize
                elif dtype == np.dtype("int16").newbyteorder('>'): 
                    wordlength = length*width*self.wordsize_s
                else:
                    raise RuntimeError("Unknown Data Type: %d" % dtype)
                array_from_buffer = \
                                bb.read_ndarray(position,
                                 wordlength,
                                 1,dtype)                      
                if width != 1:
                    array_from_buffer=array_from_buffer.reshape((length, width))         
                position += wordlength   
        else:
            array_from_buffer=[]
            for l in range(0, length):
                try:
                    string_from_buffer=str(bb.read_text(position,width))\
                    .replace("\x01", "").replace("\x00", "")
                    array_from_buffer.append(string_from_buffer)
                except:
                    array_from_buffer.append("")
                position +=width
        return array_from_buffer, position    
    
    def slicer(self, position: int, thickness: int, dtype: np.dtype):
        _all_data = self.all_data[:, position: position + thickness].flatten()
        bytes = _all_data.tobytes()
        position += thickness
        return np.frombuffer(bytes, dtype, -1, 0).reshape(int(len(_all_data)/thickness), thickness), position  
           
    def load_file(self, file: Union[str, BinaryBuffer]) -> 'RadiossReader':
        """ Load Radioss Animation File header from a Radioss Animation File file

        Parameters
        ----------
        file: Union[str, BinaryBuffer]
            path to Radioss Animation File or `BinaryBuffer` holding memory of Radioss Animation File

        Returns
        -------
        self: ReadRadioss
            returning self on success

        """

        LOGGER.debug("_load_file start")
        LOGGER.debug(f"file: {file}")

        if not isinstance(file, (str, BinaryBuffer)):
            err_msg                                     = "Argument 'file' must have type 'str' or 'lasso.io.BinaryBuffer'."
            raise ValueError(err_msg)

        # get the memory
        # TODO rewrite this
        if isinstance(file, str):
            bb                                          = self._read_file_buffer(file)
            self.n_header_bytes                         = len(bb)
        else:
            bb = file
            self.wordsize, self.itype, self.ftype       = self._determine_file_settings(bb)
            self.n_header_bytes                         = self._determine_n_bytes(bb, self.wordsize)

        LOGGER.debug(f"n_header_bytes: {self.n_header_bytes}")
        
        # Create a static array of case names
        case=                                           {    1 : "IE"   , 2 :     "KE", 3 : "XMOM", 4 : "YMOM", 5 : "ZMOM", 6 :  "MASS",
                                                             7 : "HE"   , 8 : "TURBKE", 9 :  "XCG",10 :  "YCG",11 :  "ZCG",12 : "XXMOM",
                                                             13:"YYMOM" , 14:  "ZZMOM",15 :  "IXX",16 :  "IYY",17 :  "IZZ",18 :   "IXY",
                                                             19:"IYZ"   , 20:    "IZX",21 :  "RIE",22 : "KERB",23 :"RKERB",24 :   "RKE",
                                                             25:"HEAT"  , 26:  "ZZMOM",27 :"ZZMOM",28 :"HEAT"   }    
        # Dict to list
        case_array                                      = [None] * (max(case.keys()) + 1)
        for case_index in case:
            case_array[case_index]                          = case[case_index]

        #/*-----------------------------
        #     Header Parsing  
        #-----------------------------*/                
        # Read in 2D Data
        self.array={}
        self.raw_header = {}
        position                                    = 0
        #/*-------TITRE------------*/
        __, position                                = self.read_single_word(bb, self.itype, position)
        # Time text
        thicode, position                           = self.read_single_word(bb, self.itype, position)        
        title_text, position                        = self.read_single_word(bb,[str, 80], position)
        __, position                                = self.read_single_word(bb, self.itype, position)
        # Determine width of strings
        if thicode                                  >= 4021:
            titleLength                                 = 100
        elif thicode                                >= 3041:
            titlelength                                 = 80
        else:
            titlelength                                 = 40 
        # Version and Date
        length, position                            = self.read_single_word(bb, self.itype, position)  
        ccode, position                             = self.read_single_word(bb,[str, 80], position)
        self.raw_header["date and version"]         = ccode
        length, position                            = self.read_single_word(bb, self.itype, position)         
        # Additional records
        if thicode                                  >= 3050:
            # Additional records
            __, position                                = self.read_single_word(bb, self.itype, position)  
            __, position                                = self.read_single_word(bb, self.itype, position)  
            __, position                                = self.read_single_word(bb, self.itype, position)  
            # 1st record
            __, position                                = self.read_single_word(bb, self.itype, position)  
            __, position                                = self.read_single_word(bb, self.itype, position)  
            __, position                                = self.read_single_word(bb, self.itype, position)  
            # 2nd record
            length, position                            = self.read_single_word(bb, self.itype, position) 
            rcode, position                             = self.read_single_word(bb, self.ftype, position)
            rcode, position                             = self.read_single_word(bb, self.ftype, position)
            rcode, position                             = self.read_single_word(bb, self.ftype, position)
            length, position                            = self.read_single_word(bb, self.itype, position)            
        # Heirarchy    
        __, position                                    = self.read_single_word(bb, self.itype, position) 
        npart_nthpart, position                         = self.read_single_word(bb, self.itype, position)
        nummat, position                                = self.read_single_word(bb, self.itype, position)
        numgeo, position                                = self.read_single_word(bb, self.itype, position)
        nsubs, position                                 = self.read_single_word(bb, self.itype, position)
        nthgrp2, position                               = self.read_single_word(bb, self.itype, position)
        nglob, position                                 = self.read_single_word(bb, self.itype, position)
        __, position                                    = self.read_single_word(bb, self.itype, position)  
        
        if nglob>0:            
            __, position                                    = self.read_single_word(bb, self.itype, position)         
            __, position                                    = self.add_array(nglob, 1, position, bb, self.itype)            
            __, position                                    = self.read_single_word(bb, self.itype, position)            
        #-------PART DESCRIPTION------------*/
        nvar_part                                       = []
        nvar_part_tot                                   = 0
        self.raw_header["case_names"]                   = {}
        self.raw_header["part_names"]                   = {}
        self.raw_header["case_names"]["part"]           = {}
        part_names=[]
        for i in range(0,npart_nthpart):
            __, position                                    = self.read_single_word(bb, self.itype, position)         
            __, position                                    = self.read_single_word(bb, self.itype, position)       
            part_name, position                             = self.read_single_word(bb,[str, titlelength], position)
            part_name                                       = part_name.strip()
            part_names.append(part_name)
            __, position                                    = self.read_single_word(bb, self.itype, position)      
            __, position                                    = self.read_single_word(bb, self.itype, position)      
            __, position                                    = self.read_single_word(bb, self.itype, position)      
            nvar, position                                  = self.read_single_word(bb, self.itype, position)  
            nvar_part_tot                                   += nvar
            __, position                                    = self.read_single_word(bb, self.itype, position) 
            nvar_part.append(nvar)
            nvar_part_tot                                   += nvar
            case_names                                      = [] 
            for j in range(0, nvar):
                if                                              j==0: 
                    __, position                                    = self.read_single_word(bb, self.itype, position)                    
                case_flag, position                                 = self.read_single_word(bb, self.itype, position)                
                if                                              j == nvar - 1:
                    __, position                                    = self.read_single_word(bb, self.itype, position)
                try:   
                    case_names.append(case_array[case_flag])
                except:
                    case_names.append("empty")
            self.raw_header["case_names"]\
                ["part"][i]                         = case_names                        
        self.raw_header["part_names"]               = part_names 
        
        #/*-------MATER DESCRIPTION------------*/
        if nummat > 0:   
            for i in range(0, nummat):
                __, position                                = self.read_single_word(bb, self.itype, position)
                __, position                                = self.read_single_word(bb, self.itype, position)
                __, position                                = self.read_single_word(bb,[str, titlelength], position)
                __, position                                = self.read_single_word(bb, self.itype, position)                  
        # Geo Description            
        if nummat > 0:   
            for i in range(0, numgeo):
                __, position                                = self.read_single_word(bb, self.itype, position)
                __, position                                = self.read_single_word(bb, self.itype, position)
                __, position                                = self.read_single_word(bb,[str, titlelength], position)
                __, position                                = self.read_single_word(bb, self.itype, position)                      
        # Heirarchy 
        if nsubs                                    > 0:   
            self.raw_header["case_names"]["sub"]        = {}          
            nvar_subs                                   = 0    
            sub_names                                   = []
            for i in range(0, nsubs):    
                __, position                                = self.read_single_word(bb, self.itype, position)
                idrequest, position                         = self.read_single_word(bb, self.itype, position)
                __, position                                = self.read_single_word(bb, self.itype, position)        
                nbsubsf, position                           = self.read_single_word(bb, self.itype, position)
                nbpartf, position                           = self.read_single_word(bb, self.itype, position)
                nvar, position                              = self.read_single_word(bb, self.itype, position)
                sub_name, position                          = self.read_single_word(bb,[str, titlelength], position)
                sub_name                                    = sub_name.strip()
                sub_names.append(sub_name)
                __, position                                = self.read_single_word(bb, self.itype, position)
                case_names                                  = []
                if nbsubsf                                  > 0:
                    for j in range(0, nbsubsf):
                        if                                      j == 0:
                            __, position                            = self.read_single_word(bb, self.itype, position)
                        __, position                            = self.read_single_word(bb, self.itype, position)
                        if                                      j == nbsubsf -1:
                            __, position                            = self.read_single_word(bb, self.itype, position)
                if nbpartf                                      > 0:
                    for j in range(0, nbpartf):
                        if j                                            == 0:
                            __, position                                    = self.read_single_word(bb, self.itype, position)
                        __, position                                    = self.read_single_word(bb, self.itype, position)
                        if j                                            == nbpartf -1:
                            __, position                                    = self.read_single_word(bb, self.itype, position)                                      
                if nvar                                         > 0:
                    for j in range(0, nvar):
                        if j                                            == 0:
                            __, position                                    = self.read_single_word(bb, self.itype, position)                            
                        case_flag, position                             = self.read_single_word(bb, self.itype, position)                        
                        if j                                            == nvar -1:
                            __, position                                    = self.read_single_word(bb, self.itype, position)                                
                        nvar_subs                                       += 1
                        case_names.append(case.get(case_flag, "empty"))
                    self.raw_header["case_names"]["sub"][i]         = case_names  
                    self.raw_header["sub_names"]                    = sub_names                        
        # TH
        nbelem_thgrp                                = []
        nvar_thgrp                                  = []
        self.raw_header["case_names"]["nthgrp2"]    = {}
        nthgrp2_names                               = []
        if nthgrp2                                  > 0:
            # Ensure no duplicate names
            name_set = set([])
            for i in range(0, nthgrp2):
                __, position                                = self.read_single_word(bb, self.itype, position)
                __, position                                = self.read_single_word(bb, self.itype, position)
                __, position                                = self.read_single_word(bb, self.itype, position)
                __, position                                = self.read_single_word(bb, self.itype, position)
                nbelem, position                            = self.read_single_word(bb, self.itype, position)
                nbelem_thgrp.append(nbelem)
                nvar, position                              = self.read_single_word(bb, self.itype, position)
                nvar_thgrp.append(nvar)
                _nthgrp2, position                          = self.read_single_word(bb,[str, titlelength], position)
                nthgrp2_name                                = _nthgrp2.strip()                   
                self.raw_header["case_names"]["nthgrp2"][i] = {}
                nthgrp2_names.append(nthgrp2_name)  
                __, position                                = self.read_single_word(bb, self.itype, position)
                
                _names=[]
                
                for j in range(0, nbelem):
                    
                    __, position                                = self.read_single_word(bb, self.itype, position)
                    _id, position                               = self.read_single_word(bb, self.itype, position)
                    _name, position                             = self.read_single_word(bb,[str, titlelength], position)
                    _names.append([_id, _name, nvar])
                    
                    __, position        = self.read_single_word(bb, self.itype, position)
                    
                self.raw_header["case_names"]["nthgrp2"][i]  = _names
                for j in range(0, nvar):    
                        if j == 0:
                            __, position        = self.read_single_word(bb, self.itype, position)
                        __, position         = self.read_single_word(bb, self.itype, position)
                        if j == nvar -1:
                            __, position        = self.read_single_word(bb, self.itype, position)     

        # Calculate the number of timesteps by calculating the width of a timestep and the remaining filesize
        
        # Get the width
        width = 0
        #TIME   
        width += 3
        #GLOBAL VARIABLES 
        width += 2+nglob                    
        #PART VARIABLES           
        if npart_nthpart > 0:
            if nvar_part_tot >0:
                width += 1
                for i in range(0, npart_nthpart):
                    width += nvar_part[i]
                width += 1
        #SUBSET VARIABLES 
        if nvar_subs > 0:
            width += 2 + nvar_subs
        #TH GROUP
        for i in range(0, nthgrp2):
            width += 1
            for j in range(0, nbelem_thgrp[i]):
                width += nvar_thgrp[i]
            width += 1
        n_timesteps=int((len(bb) - position)/(4*width))

        # Read everything in and create 2D float array        
        self.all_data = bb.read_ndarray(position, width * n_timesteps * 4, 1, self.ftype).reshape((n_timesteps, width))
        
        # Slice the array and build main dictionary

        position = 0
        
        # Build the global dict
        # Time
        __, position                                                    = self.slicer(position, 1, self.ftype)
        self.array["time"], position                                    = self.slicer(position, 1, self.ftype)
        __, position                                                    = self.slicer(position, 1, self.ftype)
        # Global
        self.array["global"]                                            = {}
        __, position                                                    = self.slicer(position, 1, self.ftype)
        self.array["global"]["internal_energy"], position               = self.slicer(position, 1, self.ftype)
        self.array["global"]["kinetic_energy"], position                = self.slicer(position, 1, self.ftype)
        self.array["global"]["x_momentum"], position                    = self.slicer(position, 1, self.ftype)
        self.array["global"]["y_momentum"], position                    = self.slicer(position, 1, self.ftype)
        self.array["global"]["z_momentum"], position                    = self.slicer(position, 1, self.ftype)
        self.array["global"]["mass"], position                          = self.slicer(position, 1, self.ftype)
        self.array["global"]["time_step"], position                     = self.slicer(position, 1, self.ftype)
        self.array["global"]["rotational_energy"], position             = self.slicer(position, 1, self.ftype)
        self.array["global"]["external_work"], position                 = self.slicer(position, 1, self.ftype)
        self.array["global"]["spring_energy"], position                 = self.slicer(position, 1, self.ftype)
        self.array["global"]["contact_energy"], position                = self.slicer(position, 1, self.ftype)
        self.array["global"]["hourglass_energy"], position              = self.slicer(position, 1, self.ftype)
        if nglob                                                         == 15:
            self.array["global"]["elastic contact energy"], position        = self.slicer(position, 1, self.ftype)
            self.array["global"]["frictional contact energy"], position     = self.slicer(position, 1, self.ftype)
            self.array["global"]["damping contact energy"], position        = self.slicer(position, 1, self.ftype)
        __, position                                                     = self.slicer(position, 1, self.ftype)
        # Part variables
        if npart_nthpart                                                > 0:
            if nvar_part_tot                                                > 0:
                __, position                                                    = self.slicer(position, 1, self.ftype)
            self.array["part"]                                              = {}
            for i, part_name in enumerate(self.raw_header["part_names"]):
                if part_name not in self.array["part"]:
                    self.array["part"][part_name]                                              = {}
                case_names = self.raw_header["case_names"]["part"][i]
                for case_name in case_names:
                    self.array["part"][part_name][case_name], position              = self.slicer(position, 1, self.ftype)
            __, position                                                    = self.slicer(position, 1, self.ftype)
        # Build the subset dict      
        if nvar_subs                                                    > 0:
            __, position                                                    = self.slicer(position, 1, self.ftype)
            self.array["subs"]                                              = {}
            for i, sub_name in enumerate(sub_names):
                if sub_name not in self.array["subs"]:
                    self.array["subs"][sub_name]                                    = {}
                case_names = self.raw_header["case_names"]["sub"][i]
                for case_name in case_names:
                    self.array["subs"][sub_name][case_name], position               = self.slicer(position, 1, self.ftype)
            __, position                                                    = self.slicer(position, 1, self.ftype)
        
        # Build the TH dict

        self.array["group"]                                                     = {}
        for i, nthgrp2_name in enumerate(nthgrp2_names):
            __, position                                                    = self.slicer(position, 1, self.ftype)
            name                                                                    = nthgrp2_name.strip().lower()
            if name not in self.array["group"]:
                self.array["group"][name]                                               = {}
            for TH in self.raw_header["case_names"]["nthgrp2"][i]:
                self.array["group"][name][TH[0]]                                        = {}
                self.array["group"][name][TH[0]]\
                    ["name"]                                                            = TH[1].strip().lower()
                self.array["group"][name][TH[0]]\
                    ["time_history"], position                                          = self.slicer(position, TH[2], self.ftype)
            __, position                                                    = self.slicer(position, 1, self.ftype)
        
        return self
    
    def read_words(self,
                   bb: BinaryBuffer,
                   words_to_read: dict,
                   storage_dict: dict = None):
        ''' Read several words described by a dict

        Parameters
        ----------
        words_to_read: dict
            this dict describes the words to be read. One entry
            must be a tuple of len two (byte position and dtype)
        storage_dict: dict
            in this dict the read words will be saved

        Returns
        -------
        storage_dict: dict
            the storage dict given as arg or a new dict if none was given
        '''

        if storage_dict is None:
            storage_dict = {}

        for name, data in words_to_read.items():

            # check buffer length
            if data[0] >= len(bb):
                continue

            # read data
            if data[1] == self.itype:
                storage_dict[name] = int(bb.read_number(
                    data[0], data[1]))
            elif data[1] == self.ftype:
                storage_dict[name] = float(bb.read_number(
                    data[0], data[1]))   
            elif data[1] == np.dtype('V4'):
                storage_dict[name] = str(bb.read_number(
                    data[0], data[1]))             
            elif data[1] == str:
                try:
                    storage_dict[name] = str(bb.read_text(
                        data[0], data[2])).replace("\x00", "")
                except UnicodeDecodeError:
                    storage_dict[name] = ""

            else:
                raise RuntimeError(
                    "Encountered unknown dtype {} during reading.".format(str(data[1])))

        return storage_dict
    
    def read_single_word(self, bb: BinaryBuffer, dtype: (np.dtype, list), position: int):
        
        ''' Readsingle words described by a dtype

        Parameters
        ----------
        dtype: np.dtype or list for string
        position: int

        Returns
        -------
        val, position
            val is str, int or float
            position is int
        '''
        # if dtype is list, string type list[1] is word width
        if isinstance(dtype, list):
            #try:
            val = bb.read_text(position,dtype[1])
            #except:
                #val=""
            wordlength = dtype[1]    
        elif dtype== self.itype:
            val = int(bb.read_number(
                position, dtype))
            wordlength = self.wordsize
        elif dtype == self.ftype:
            val = float(bb.read_number(
                position, dtype))   
            wordlength = self.wordsize
        elif dtype == np.dtype('V4'):
            val = str(bb.read_number(
                position, dtype))    
            wordlength = self.wordsize
        else:
            raise RuntimeError(
                "Encountered unknown dtype {} during reading.".format(str(dtype)))

        position += wordlength
          
            
        return val, position
