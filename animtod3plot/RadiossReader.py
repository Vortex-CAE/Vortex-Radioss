#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

from typing import Union

import numpy as np
import rich

from lasso.io.binary_buffer import BinaryBuffer
from lasso.logging import get_logger

LOGGER = get_logger(__file__)



class RadiossReader:
    ''' Class for reading a Radioss Animation File
    '''

    # meta
    filepath: str = ""  #: filepath of file

    # file info
    itype: np.dtype     = np.dtype("int32").newbyteorder('>') #: integer type of Radioss Animation File
    itype_s: np.dtype   = np.dtype("int16").newbyteorder('>') #: short integer type of Radioss Animation File
    btype: np.dtype     = bool #: boolean type of Radioss Animation File
    ftype: np.dtype     = np.dtype("float32").newbyteorder('>') #: floating point type of Radioss Animation File
    wordsize: int       = 4   #: Size for ftype, itype
    wordsize_s: int     = 2 #: Size for type itype_s
    wordsize_b: int     = 1 #: Size for type bool
    
    """
    - 3D GEOMETRY
    nbElts3D                        - number of 8nodes elements
    nbParts3D                       - number of parts
    nbEFunc3D                       - number of element scalar values
    nbTens3D                        - number of tensors
    connect3DA                      - connectivity array of 8  nbElts3D
    delElt3DA                       - are the elements deleted or not, array of nbElts3D
    defPart3DA                      - part definition: array of nbParts3D
    pText3DA                        - part name: array of nbParts3D
    fText3DA                        - array of scalar name: nbEFunc3D
    eFunc3DA                        - scalar value per element array of nbEFunc3DnbElts3D
    tText3DA                        - tensor name array of nbTens3D
    tensVal3DA                      - tens value array of 6nbTens3DnbElt3D
    eMass3DA                        - nbElt3D, elt mass
    elNum3DA                        - nbElt3D, intern elt numbering

    - 2D GEOMETRY
    nbFacets                        - number of 4nodes elements
    nbNodes                         - total number of nodes
    nbParts                         - number of parts
    nbFunc                          - number of nodal scalar values
    nbEFunc                         - number of element scalar values
    nbVect                          - number of vectors
    nbTens                          - number of tensors
    nbSkew                          - number of skews
    skewValA                        - array of the skew values for each elt
    coorA                           - coordinates array of 3nbNodes
    connectA                        - connectivity array of 4  nbFacets
    delEltA                         - are the elements deleted or not, array of nbFacets
    defPartA                        - part definition: array of nbParts
    pTextA                          - part name: array of nbParts
    uint16_t normShortA             - facet normal in uint16_t : array of 3nbNodes
    normFloatA                      - facet normal in : array of 3nbNodes
    fTextA                          - array of scalar name: nbFunc+nbEFunc
    funcA                           - scalar value per node array of nbFuncnbNodes
    eFuncA                          - scalar value per element array of nbEFuncnbFacets
    vTextA                          - vect name array of nbVect
    vectValA                        - vect val array of 3nbVectnbNodes
    tTextA                          - tensor name array of nbTens
    tensValA                        - tens value array of 3nbTensnbElt
    nMassA = nullptr, eMassA        - nbNodes, nbElt, nodal, elt mass
    nodNumA = nullptr, elNumA       - nbNodes, nbElt, intern node/elt nb

    - 1D geometry
    nbElts1D                        - number of 2nodes elements
    nbParts1D                       - number of parts
    nbEFunc1D                       - number of lin. elt scalar values
    nbTors1D                        - number of torseur values
    isSkew1D                        - number of skews
    connect1DA                      - element connectivity array
    delElt1DA                       - delEltA indicates which elts are deleted or not
    defPart1DA                      - parts definition: array
    pText1DA                        - part names array
    fText1DA                        - array of scalar function names
    eFunc1DA                        - array of element scalar values
    tText1DA                        - array of tensor names
    torsVal1DA                      - read the array of 3 forces, 3 moments torsor values for each elt
    elt2Skew1DA                     - array of the skew number for each elt
    eMass1DA                        - mass :elementar mass
    elNum1DA                        - element numbering

    - hierarchy
    nbSubsets                       - number of subsets
    subsetText                      - subset name
    numParent                       - parent number
    nbSubsetSon                     - nb subset sons
    subsetSonA                      - subset son list
    nbSubPart2D                     - nb 2D subparts
    subPart2DA                      - 2D subpart list
    nbSubPart3D                     - nb 3D subparts
    subPart3DA                      - 3D subpart list
    nbSubPart1D                     - nb 1D subparts
    subPart1DA                      - 1D subpart list
    nbMaterials                     - material number
    nbProperties                    - number of Properties
    materialTextA                   - material names
    propertiesTextA                 - property names
    materialTypeA                   - material types
    propertiesTypeA                 - property types
    part2subset2DA                  - array of subset for each part
    partMaterial2DA                 - array of material for each part
    partProperties2DA               - array of properties for each part
    part2subset3DA                  - array of subset for each part
    partMaterial3DA                 - array of material for each part
    partProperties3DA               - array of properties for each part
    part2subset1DA                  - array of subset for each part
    partMaterial1DA                 - array of material for each part
    partProperties1DA               - array of properties for each part

    - NODES/ELTS FOR Time History
    nbNodesTH                       - number of Time History nodes
    nbElts2DTH                      - number of Time History 2D elements
    nbElts3DTH                      - number of Time History 3D elements
    nbElts1DTH                      - number of Time History 1D elements
    nodes2THA                       - node list
    n2thTextA                       - node names
    elt2DTHA                        - elt 2D list
    elt2DthTextA                    - elt 2D name
    elt3DTHA                        - elt 3D list
    elt3DthTextA                    - elt 3D name
    elt1DTHA                        - elt 1D list
    elt1DthTextA                    - elt 1D name

    - SPH ELTS
    nbEltsSPH                       - number of elements
    nbPartsSPH                      - number of parts
    nbEFuncSPH                      - number of scalar values
    nbTensSPH                       - number of tensors
    connecSPH                       - connectivity array
    delEltSPH                       - are the elements deleted or not
    defPartSPH                      - part definition
    eFuncSPH                        - scalar values
    scalTextSPH                     - scalar names
    tensTextSPH                     - tensor names
    tensValSPH                      - tensor values
    pTextSP                         - part names
    eMassSPH                        - element mass
    nodNumSPH                       - element numbering
    numParentSPH                    - parent number
    matPartSPH                      - material number
    propPartSPH                     - property number
    """

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
            
            
            """import chardet
            #with open(filepath, 'rb') as file:
            with open(filepath, 'rb') as x:
                line = x.readline()
                curChar = chardet.detect(line)
                print(curChar)
                while line:
                    if curChar != chardet.detect(line):
                        curChar = chardet.detect(line)
                        print(curChar)
                    line = x.readline()"""
            
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
                if dtype == np.dtype("int32").newbyteorder('>') or dtype == np.dtype("float32").newbyteorder('>'):
                    wordlength = length*width*self.wordsize
                elif dtype == np.dtype("int16").newbyteorder('>'): 
                    wordlength = length*width*self.wordsize_s
                elif dtype == np.dtype(bool):
                    wordlength = length*width*self.wordsize_b
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
                string_from_buffer  = str(bb.read_text(position,width))
                string_from_buffer  = string_from_buffer.replace("\x01", "")
                string_from_buffer  = string_from_buffer.replace("\x00", "")
                array_from_buffer.append(string_from_buffer)
                position += width
        return array_from_buffer, position    
           
    def load_file(self, file: Union[str, BinaryBuffer]) -> 'RadiossReader':
        
        import time
        time_start=time.time()
        
        """ Load Radioss Animation File header from a Radioss Animation File file

        Parameters
        ----------
        file: Union[str, BinaryBuffer]
            path to Radioss Animation File or `BinaryBuffer` holding memory of Radioss Animation File

        Returns
        -------
        self: ReadRadioss
            returning self on success

        Notes
        -----
            This routine only loads the minimal amount of data
            that is neccessary. Thus it is safe to use on huge files.

        Examples
        --------
            >>> header = ReadRadioss().load_file("path/to/Radioss Animation File")
            >>> header.n_shells
            19684
        """

        LOGGER.debug("_load_file start")
        LOGGER.debug(f"file: {file}")

        if not isinstance(file, (str, BinaryBuffer)):
            err_msg = "Argument 'file' must have type 'str' or 'lasso.io.BinaryBuffer'."
            raise ValueError(err_msg)

        # get the memory
        # TODO rewrite this
        if isinstance(file, str):
            bb                  = self._read_file_buffer(file)
            self.n_header_bytes = len(bb)
        else:
            bb = file
            self.wordsize, self.itype, self.ftype   = self._determine_file_settings(bb)
            self.n_header_bytes                     = self._determine_n_bytes(bb, self.wordsize)

        LOGGER.debug(f"n_header_bytes: {self.n_header_bytes}")
        
        # Read in 2D Data
        self.raw_arrays     = {}
        self.raw_header     = {}
        position            = 0
        
        # Process the binary
        
        # File structure version
        self.raw_header["magic"], position              = self.read_single_word(bb, self.itype, position) 
        
        # Time
        self.raw_header["time"], position               = self.read_single_word(bb, self.ftype, position) 
        
        # Time text
        self.raw_header["time_text"], position          = self.read_single_word(bb,[str, 81], position)
        
        # Mod anim text
        self.raw_header["mod_anim_text"], position      = self.read_single_word(bb,[str, 81], position)
        
        # Radioss Run text
        self.raw_header["radioss_run_text"], position   = self.read_single_word(bb,[str, 81], position)
        
        # Defines if theflagA mass is saved or not
        self.raw_header["flag_a_0"], position           = self.read_single_word(bb, self.itype, position)
        
        # Defines if the node-element numbering arrays are saved or not
        self.raw_header["flag_a_1"], position           = self.read_single_word(bb, self.itype, position)
        
        # Defines format if there is 3D geometry
        self.raw_header["flag_a_2"], position           = self.read_single_word(bb, self.itype, position)
        
        # Defines format if there is 1D geometry
        self.raw_header["flag_a_3"], position           = self.read_single_word(bb, self.itype, position)
        
        # Defines Hierarchy
        self.raw_header["flag_a_4"], position           = self.read_single_word(bb, self.itype, position)
        
        # Defines node/elt list for T/H
        self.raw_header["flag_a_5"], position           = self.read_single_word(bb, self.itype, position)
        
        # Defines if there is a new skew for tensor 2D
        self.raw_header["flag_a_6"], position           = self.read_single_word(bb, self.itype, position)
        
        # Define if there is SPH format
        self.raw_header["flag_a_7"], position           = self.read_single_word(bb, self.itype, position)
        
        # Unused
        self.raw_header["flag_a_8"], position           = self.read_single_word(bb, self.itype, position)
        
        # Unused
        self.raw_header["flag_a_9"], position           = self.read_single_word(bb, self.itype, position)
               
        #********************
        # 2D GEOMETRY
        #********************    
        
        # Node number
        self.raw_header["nbNodes"], position            = self.read_single_word(bb, self.itype, position)        
        
        # Number of shells
        self.raw_header["nbFacets"], position           = self.read_single_word(bb, self.itype, position)        
        
        # Number of parts
        self.raw_header["nbParts"], position            = self.read_single_word(bb, self.itype, position)        
        
        # Number of nodal scalar values
        self.raw_header["nbFunc"], position             = self.read_single_word(bb, self.itype, position)        
        
        # Number of shell scalar values
        self.raw_header["nbEFunc"], position            = self.read_single_word(bb, self.itype, position)        
        
        # Number of node vector values
        self.raw_header["nbVect"], position             = self.read_single_word(bb, self.itype, position)    
        
        # Number of shell tensor values
        self.raw_header["nbTens"], position             = self.read_single_word(bb, self.itype, position)        
        
        # Number of shell skew array
        self.raw_header["nbSkew"], position             = self.read_single_word(bb, self.itype, position)                

        # Read in nbSkew Array - whatever that is
        self.raw_arrays["skewValA"] , position          = self.add_array(self.raw_header["nbSkew"], 6, position, bb, self.itype_s)
        
        # Read in Nodal Coordinates
        self.raw_arrays["coorA"] , position             = self.add_array(self.raw_header["nbNodes"], 3, position, bb, self.ftype)
        
        # Read in Element Shell Connectivity        
        self.raw_arrays["connectA"] , position          = self.add_array(self.raw_header["nbFacets"], 4, position, bb, self.itype)
        
        # Read in Deleted Elements
        self.raw_arrays["delEltA"] , position           = self.add_array(self.raw_header["nbFacets"], 1, position, bb, self.btype)

        # Read in Part Index
        self.raw_arrays["defPartA"] , position          = self.add_array(self.raw_header["nbParts"], 1, position, bb, self.itype)
        
        # Read in Part Names     
        self.raw_arrays["pTextA"] , position            = self.add_array(self.raw_header["nbParts"], 50, position, bb, str)
        
        # Read in Norm Value for each Node
        self.raw_arrays["normFloatA"] , position        = self.add_array(self.raw_header["nbNodes"], 3, position, bb, self.itype_s)     
        
        # Array of the total scalar function names
        self.raw_arrays["fTextA"] , position            = self.add_array(self.raw_header["nbFunc"] + self.raw_header["nbEFunc"], 81, position, bb, str)
        
        # Read in FuncA Scalar Values
        self.raw_arrays["funcA"] , position             = self.add_array(self.raw_header["nbNodes"] * self.raw_header["nbFunc"], 1, position, bb, self.ftype)       
        
        # Read in FuncA Scalar Values
        self.raw_arrays["eFuncA"] , position            = self.add_array(self.raw_header["nbFacets"] * self.raw_header["nbEFunc"], 1, position, bb, self.ftype)     
        
        # Read in Vector names
        self.raw_arrays["vTextA"] , position            = self.add_array(self.raw_header["nbVect"], 81, position, bb, str)
        
        # Read in Vector Values
        self.raw_arrays["vectValA"] , position          = self.add_array(self.raw_header["nbNodes"] * self.raw_header["nbVect"], 3, position, bb, self.ftype)    
        
        # Tensor Names
        self.raw_arrays["tTextA"] , position            = self.add_array(self.raw_header["nbTens"], 81, position, bb, str)
        
        # Tensor Values
        self.raw_arrays["tensValA"] , position          = self.add_array(self.raw_header["nbFacets"] * self.raw_header["nbTens"], 3, position, bb, self.ftype) 
        
        # mass : elemental and nodal masses
        if self.raw_header["flag_a_0"] == 1:
            self.raw_arrays["eMassA"] , position            = self.add_array(self.raw_header["nbFacets"], 1, position, bb, self.ftype) 
            self.raw_arrays["nMassA"] , position            = self.add_array(self.raw_header["nbNodes"], 1, position, bb, self.ftype)
            
        # Node and Element Numbering
        if self.raw_header["flag_a_1"]:
            self.raw_arrays["nodNumA"] , position           = self.add_array(self.raw_header["nbNodes"], 1, position, bb, self.itype) 
            self.raw_arrays["elNumA"] , position            = self.add_array(self.raw_header["nbFacets"], 1, position, bb, self.itype) 
            
        # Heirarchy 
        if self.raw_header["flag_a_4"]:
            self.raw_arrays["part2subset2DA"], position     = self.add_array(self.raw_header["nbParts"], 1, position, bb, self.itype)
            self.raw_arrays["partMaterial2DA"] , position   = self.add_array(self.raw_header["nbParts"], 1, position, bb, self.itype)
            self.raw_arrays["partProperties2DA"] , position = self.add_array(self.raw_header["nbParts"], 1, position, bb, self.itype)
            
        # 3D Geometry
        # If there are solid elements
        if self.raw_header["flag_a_2"]:
            # Number of solids
            self.raw_header["nbElts3D"] , position          = self.read_single_word(bb,self.itype, position)
            
            # Number of solid parts
            self.raw_header["nbParts3D"] , position         = self.read_single_word(bb,self.itype, position)
            
            # Number of solid scalar values
            self.raw_header["nbEFunc3D"] , position         = self.read_single_word(bb,self.itype, position)
            
            # Number of solid tensor values
            self.raw_header["nbTens3D"] , position          = self.read_single_word(bb,self.itype, position)
            
            # Element solid connecivity
            self.raw_arrays["connect3DA"], position         = self.add_array(self.raw_header["nbElts3D"], 8, position, bb, self.itype)
            
            # Element solid eroded array
            self.raw_arrays["delElt3DA"] , position         = self.add_array(self.raw_header["nbElts3D"], 1, position, bb, self.btype)     
            
            # Element solid part indexes - needs unpacking
            self.raw_arrays["defPart3DA"] , position        = self.add_array(self.raw_header["nbParts3D"], 1, position, bb, self.itype)
            
            # Element solid part names
            self.raw_arrays["pText3DA"] , position          = self.add_array(self.raw_header["nbParts3D"], 50, position, bb, str)
                
            # Element solid scalar names
            self.raw_arrays["fText3DA"] , position          = self.add_array(self.raw_header["nbEFunc3D"], 81, position, bb, str)
            
            # Element solid scalar values
            self.raw_arrays["eFunc3DA"] , position          = self.add_array(self.raw_header["nbEFunc3D"] * self.raw_header["nbElts3D"], 1, position, bb, self.ftype)
           
            # Element solid tensor names
            self.raw_arrays["tText3DA"] , position          = self.add_array(self.raw_header["nbTens3D"], 81, position, bb, str)
            
            # Element solid tensor values
            self.raw_arrays["tensVal3DA"] , position        = self.add_array(self.raw_header["nbElts3D"] * self.raw_header["nbTens3D"], 6, position, bb, self.ftype)
            
            # If mass data was stored
            if self.raw_header["flag_a_0"] == 1:
                self.raw_arrays["eMass3DA"] , position      = self.add_array(self.raw_header["nbElts3D"], 1, position, bb, self.ftype)

            # If element numbering was stored
            if self.raw_header["flag_a_1"] == 1:
                self.raw_arrays["elNum3DA"] , position      = self.add_array(self.raw_header["nbElts3D"], 1, position, bb, self.itype)

            # If Hierchy data is stored
            if self.raw_header["flag_a_4"]:
                # Array of subset for each part
                self.raw_arrays["part2subset3DA"] , position    = self.add_array(self.raw_header["nbParts3D"], 1, position, bb, self.itype)
                # Array of material for each part
                self.raw_arrays["partMaterial3DA"] , position   = self.add_array(self.raw_header["nbParts3D"], 1, position, bb, self.itype)
                # Array of properties for each part
                self.raw_arrays["partProperties3DA"] , position = self.add_array(self.raw_header["nbParts3D"], 1, position, bb, self.itype)

        else:

            # Number of 3D elements
            self.raw_header["nbElts3D"]     = 0
            
            # Number of 3D parts
            self.raw_header["nbParts3D"]    = 0
            
            # Number of solid scalar values
            self.raw_header["nbEFunc3D"]    = 0   
            
            # Number of solid tensor values
            self.raw_header["nbTens3D"]     = 0         
        
        # ********************
        # 1D GEOMETRY
        # ********************
        
        # If there is 1D Geometry
        if self.raw_header["flag_a_3"]:
            
            # Number of 1D elements
            self.raw_header["nbElts1D"] , position          = self.read_single_word(bb,self.itype, position)
            
            # Number of 1D parts
            self.raw_header["nbParts1D"] , position         = self.read_single_word(bb,self.itype, position)
            
            # Number of 1D scalar values
            self.raw_header["nbEFunc1D"] , position         = self.read_single_word(bb,self.itype, position)
            
            # Number of 1D tensor values
            self.raw_header["nbTors1D"] , position          = self.read_single_word(bb,self.itype, position)
            
            # # Number for 1D skew arrays
            self.raw_header["isSkew1D"] , position          = self.read_single_word(bb,self.itype, position)

            # Element 1D connectivity array
            self.raw_arrays["connect1DA"] , position        = self.add_array(self.raw_header["nbElts1D"], 2, position, bb, self.itype)
            
            # Element 1D erosion array
            self.raw_arrays["delElt1DA"] , position         = self.add_array(self.raw_header["nbElts1D"], 1, position, bb, self.btype)
            
            # Element 1D part index array - needs unpacking            
            self.raw_arrays["defPart1DA"] , position        = self.add_array(self.raw_header["nbParts1D"], 1, position, bb, self.itype)
            
            # 1D part names
            self.raw_arrays["pText1DA"] , position          = self.add_array(self.raw_header["nbParts1D"], 50, position, bb, str)
            
            # 1D scalar names
            self.raw_arrays["fText1DA"] , position          = self.add_array(self.raw_header["nbEFunc1D"], 81, position, bb, str)
            
            # 1D scalar values
            self.raw_arrays["eFunc1DA"] , position          = self.add_array(self.raw_header["nbEFunc1D"] * self.raw_header["nbElts1D"], 1, position, bb, self.ftype)
            
            # 1D tensor names
            self.raw_arrays["tText1DA"] , position          = self.add_array(self.raw_header["nbTors1D"], 81, position, bb, str)
            
            # 1D tensor values
            self.raw_arrays["torsVal1DA"] , position        = self.add_array(self.raw_header["nbElts1D"] * self.raw_header["nbTors1D"], 9, position, bb, self.ftype)
            
            # 1D skew number
            if self.raw_header["isSkew1D"]:
                self.raw_arrays["elt2Skew1DA"] , position   = self.add_array(self.raw_header["nbElts1D"], 1, position, bb, self.itype)
                
            # 1D mass data array
            if self.raw_header["flag_a_0"]==1:
                self.raw_arrays["eMass1DA"] , position      = self.add_array(self.raw_header["nbElts1D"], 1, position, bb, self.ftype)
            
            # 1D element numbering array
            if self.raw_header["flag_a_1"]==1:
                self.raw_arrays["elNum1DA"] , position      = self.add_array(self.raw_header["nbElts1D"], 1, position, bb, self.itype)
                
            # If Heirarchy data is present
            if self.raw_header["flag_a_4"]:
                
                # 1D subset array
                self.raw_arrays["part2subset1DA"] , position    = self.add_array(self.raw_header["nbParts1D"], 1, position, bb, self.itype)
                
                # 1D part array
                self.raw_arrays["partMaterial1DA"] , position   = self.add_array(self.raw_header["nbParts1D"], 1, position, bb, self.itype)
                
                # 1D properties array
                self.raw_arrays["partProperties1DA"] , position = self.add_array(self.raw_header["nbParts1D"], 1, position, bb, self.itype)
        
        else:
            # Number of 1D elements
            self.raw_header["nbElts1D"]     = 0
            
            # Number of 1D parts
            self.raw_header["nbParts1D"]    = 0
            
            # Number of 1D scalar values
            self.raw_header["nbEFunc1D"]    = 0  
            
            # Number of 1D tensor values
            self.raw_header["nbTors1D"]     = 0
            
            # # Number for 1D skew arrays
            self.raw_header["isSkew1D"]     = 0   
            
            # # Number for 1D skew arrays
            self.raw_header["isSkew1D"]     = 0         
            
        # ********************
        # Parse the Heirarchy 
        # ********************      
        
        if self.raw_header["flag_a_4"]:  
            
            # Number of subsets
            self.raw_header["nbSubsets"] , position = self.read_single_word(bb,self.itype, position)
            self.raw_arrays["nbSubsets"]            = {}
            for _nbSubsets in range(0,self.raw_header["nbSubsets"]):

                # Subset name
                subsetText , position = self.read_single_word(bb,[str, 50], position)
                self.raw_arrays["nbSubsets"][subsetText]={}
                
                # Number of parent
                self.raw_arrays["nbSubsets"][subsetText]["numParent"], position         = self.read_single_word(bb,self.itype, position)

                # Number of Sons
                self.raw_arrays["nbSubsets"][subsetText]["nbSubsetSon"], position       = self.read_single_word(bb, self.itype, position)
                
                # Son array
                if self.raw_arrays["nbSubsets"][subsetText]["nbSubsetSon"]:
                    self.raw_arrays["nbSubsets"][subsetText]["subsetSonA"] , position   = self.add_array(self.raw_arrays["nbSubsets"][subsetText]["nbSubsetSon"], 1, position, bb, self.itype)
                
                # Number of 2D subset parts
                self.raw_arrays["nbSubsets"][subsetText]["nbSubPart2D"], position       = self.read_single_word(bb, self.itype, position) 
                
                # 2D subset part array
                if self.raw_arrays["nbSubsets"][subsetText]["nbSubPart2D"]:
                    self.raw_arrays["nbSubsets"][subsetText]["subPart2DA"] , position   = self.add_array(self.raw_arrays["nbSubsets"][subsetText]["nbSubPart2D"], 1, position, bb, self.itype)
                    
                # Number of 3d subsets
                self.raw_arrays["nbSubsets"][subsetText]["nbSubPart3D"], position       = self.read_single_word(bb, self.itype, position) 
                
                # 3D subset part array
                if self.raw_arrays["nbSubsets"][subsetText]["nbSubPart3D"]:
                    self.raw_arrays["nbSubsets"][subsetText]["subPart3DA"] , position   = self.add_array(self.raw_arrays["nbSubsets"][subsetText]["nbSubPart3D"], 1, position, bb, self.itype)
                    
                # Number of 1D subset parts
                self.raw_arrays["nbSubsets"][subsetText]["nbSubPart1D"], position       = self.read_single_word(bb, self.itype, position) 
                
                # 1D subset part array
                if self.raw_arrays["nbSubsets"][subsetText]["nbSubPart1D"]:
                    self.raw_arrays["nbSubsets"][subsetText]["subPart1DA"] , position   = self.add_array(self.raw_arrays["nbSubsets"][subsetText]["nbSubPart1D"], 1, position, bb, self.itype)    
                                        
            # Number of materials
            self.raw_header["nbMaterials"] , position       = self.read_single_word(bb,self.itype, position)
            
            # Number of properties
            self.raw_header["nbProperties"] , position      = self.read_single_word(bb,self.itype, position)
            
            # Names of materials
            self.raw_arrays["materialTextA"] , position     = self.add_array(self.raw_header["nbMaterials"], 50, position, bb, str)
            
            # Material type
            self.raw_arrays["materialTypeA"], position      = self.add_array(self.raw_header["nbMaterials"], 1, position, bb, self.itype)
            
            # Property name
            self.raw_arrays["propertiesTextA"] , position   = self.add_array(self.raw_header["nbProperties"], 50, position, bb, str)
            
            # Property type array
            self.raw_arrays["propertiesTypeA"], position    = self.add_array(self.raw_header["nbProperties"], 1, position, bb, self.itype)
            
        else:
            # Number of materials
            self.raw_header["nbMaterials"]  = 0
            
            # Number of properties
            self.raw_header["nbProperties"] = 0  
        

        # ********************
        # NODES/ELTS FOR Time History ( nodes & elems that are also selected for Time History output)
        # ********************
        
        if self.raw_header["flag_a_5"]:  

            # Number of time history nodes
            self.raw_header["nbNodesTH"] , position     = self.read_single_word(bb,self.itype, position)
            
            # Number of time history shells
            self.raw_header["nbElts2DTH"] , position    = self.read_single_word(bb,self.itype, position)
            
            # Number of time history solids
            self.raw_header["nbElts3DTH"] , position    = self.read_single_word(bb,self.itype, position)
            
            # Number of time history 1D
            self.raw_header["nbElts1DTH"] , position    = self.read_single_word(bb,self.itype, position)
            
            # Node array
            self.raw_arrays["nodes2THA"], position      = self.add_array(self.raw_header["nbNodesTH"], 1, position, bb, self.itype)
            
            # Node names
            self.raw_arrays["n2thTextA"] , position     = self.add_array(self.raw_header["nbNodesTH"], 50, position, bb, str)
            
            # Shell array
            self.raw_arrays["elt2DTHA"], position       = self.add_array(self.raw_header["nbElts2DTH"], 1, position, bb, self.itype)
            
            # Shell names
            self.raw_arrays["elt2DthTextA"] , position  = self.add_array(self.raw_header["nbElts2DTH"], 50, position, bb, str)
            
            # Solid array
            self.raw_arrays["elt3DTHA"], position       = self.add_array(self.raw_header["nbElts3DTH"], 1, position, bb, self.itype)
            
            # Solid names
            self.raw_arrays["elt3DthTextA"] , position  = self.add_array(self.raw_header["nbElts3DTH"], 50, position, bb, str)
            
            # 1D array
            self.raw_arrays["elt1DTHA"], position       = self.add_array(self.raw_header["nbElts1DTH"], 1, position, bb, self.itype)
            
            # 1D names
            self.raw_arrays["elt1DthTextA"] , position  = self.add_array(self.raw_header["nbElts1DTH"], 50, position, bb, str)
            
        else:
            
            # Number of time history nodes
            self.raw_header["nbNodesTH"]    = 0
            
            # Number of time history shells
            self.raw_header["nbElts2DTH"]   = 0
            
            # Number of time history solids
            self.raw_header["nbElts3DTH"]   = 0
            
            # Number of time history 1D
            self.raw_header["nbElts1DTH"]   = 0      
            
            

        # ********************
        # READ SPH PART */
        # ********************
        if self.raw_header["flag_a_7"]:  

            # Number of sph elements
            self.raw_header["nbEltsSPH"] , position     = self.read_single_word(bb,self.itype, position)
            
            # Number of sph parts
            self.raw_header["nbPartsSPH"] , position    = self.read_single_word(bb,self.itype, position)
            
            # Number of sph scalars
            self.raw_header["nbEFuncSPH"] , position    = self.read_single_word(bb,self.itype, position)
            
            # Number of sph tensors
            self.raw_header["nbTensSPH"] , position     = self.read_single_word(bb,self.itype, position)    
                        
            if self.raw_header["nbEltsSPH"]:
                
                # sph connectivity array
                self.raw_arrays["connecSPH"], position      = self.add_array(self.raw_header["nbEltsSPH"], 1, position, bb, self.itype)
                
            if self.raw_header["nbPartsSPH"]:
                
                # Sph part indexes - needs unpacking
                self.raw_arrays["defPartSPH"], position     = self.add_array(self.raw_header["nbPartsSPH"], 1, position, bb, self.itype)
                
                # Sph part names
                self.raw_arrays["pTextSPH"] , position      = self.add_array(self.raw_header["nbPartsSPH"], 50, position, bb, str)

            if self.raw_header["nbEFuncSPH"]:

                # Sph scalar names
                self.raw_arrays["scalTextSPH"] , position   = self.add_array(self.raw_header["nbEFuncSPH"], 81, position, bb, str)
                
                # Sph scalar array
                self.raw_arrays["eFuncSPH"], position       = self.add_array(self.raw_header["nbEltsSPH"] * self.raw_header["nbEFuncSPH"], 1, position, bb, self.itype)
            
            if self.raw_header["nbTensSPH"]:

                # Sph tensor names
                self.raw_arrays["tensTextSPH"] , position   = self.add_array(self.raw_header["nbTensSPH"], 81, position, bb, str)
                
                # Sph tensor array
                self.raw_arrays["tensValSPH"], position     = self.add_array(self.raw_header["nbEltsSPH"] * self.raw_header["nbTensSPH"], 6, position, bb, self.ftype)
                
            if self.raw_header["flag_a_0" == 1]:
                
                # Sph mass
                self.raw_arrays["eMassSPH"], position       = self.add_array(self.raw_header["nbEltsSPH"], 1, position, bb, self.ftype)

            if self.raw_header["flag_a_1" == 1]:  

                # Sph numbering
                self.raw_arrays["nodNumSPH"], position      = self.add_array(self.raw_header["nbEltsSPH"], 1, position, bb, self.ftype)


            if self.raw_header["flag_a_4"]:       
                
                # ********************
                # Parse the Heirarchy 
                # ******************** 
                
                # Sph parent number array
                self.raw_arrays["numParentSPH"], position   = self.add_array(self.raw_header["nbPartsSPH"], 1, position, bb, self.itype)
                
                # Sph material number array
                self.raw_arrays["matPartSPH"], position     = self.add_array(self.raw_header["nbPartsSPH"], 1, position, bb, self.itype)
                
                # Sph property number array
                self.raw_arrays["propPartSPH"], position    = self.add_array(self.raw_header["nbPartsSPH"], 1, position, bb, self.itype)
        else:
            
            # Number of sph elements
            self.raw_header["nbEltsSPH"]    = 0
            
            # Number of sph parts
            self.raw_header["nbPartsSPH"]   = 0
            
            # Number of sph scalars
            self.raw_header["nbEFuncSPH"]   = 0
            
            # Number of sph tensors
            self.raw_header["nbTensSPH"]    = 0             
                
                    
        # *******************************************************
        # Unpack raw_arrays into arrays
        # *******************************************************
        self.arrays={}
                
        self.arrays["timesteps"]                        = self.raw_header["time"]
        
        if "coorA" in self.raw_arrays:
            self.arrays["node_coordinates"]             = np.array(self.raw_arrays["coorA"])
        
        # Element nodes indexes are unchanged
        if "connect1DA" in self.raw_arrays:        
            self.arrays["element_beam_node_indexes"]    = np.array(self.raw_arrays["connect1DA"])
        if "connectA" in self.raw_arrays:
            self.arrays["element_shell_node_indexes"]   = np.array(self.raw_arrays["connectA"])
        if "connect3DA" in self.raw_arrays:
            self.arrays["element_solid_node_indexes"]   = np.array(self.raw_arrays["connect3DA"])
        if "connectSPH" in self.raw_arrays:
            self.arrays["element_sph_node_indexes"]     = np.array(self.raw_arrays["connecSPH"])
                        
        # Node ids are unchanged
        if "nodNumA" in self.raw_arrays:    
            self.arrays["node_ids"] = np.array(self.raw_arrays["nodNumA"])
                  
        # Nodal scalars are unpacked
        if "fTextA" in self.raw_arrays:        
            for ifun in range(0, self.raw_header["nbFunc"]):
                node_scalar     =   "node_" + str(self.raw_arrays["fTextA"][ifun]).lower().replace(" ", "_").strip()
                start           =   ifun*self.raw_header["nbNodes"]
                end             =   (ifun + 1) *self.raw_header["nbNodes"]
                tmp_list        =   self.raw_arrays["funcA"][start : end]   
                self.arrays[node_scalar] = np.array(tmp_list)

        # Nodal vectors are unpacked
        if "vTextA" in self.raw_arrays:        
            for ivect in range(0, self.raw_header["nbVect"]):
                node_vector     = "node_" + str(self.raw_arrays["vTextA"][ivect]).lower().replace(" ", "_").strip()
                start           = ivect * self.raw_header["nbNodes"]
                end             = (ivect+1) * self.raw_header["nbNodes"]
                tmp_list        = self.raw_arrays["vectValA"][start : end]
                self.arrays[node_vector] = np.array(tmp_list)       

        # Element 1D ids are unchanged
        if "elNum1DA" in self.raw_arrays:   
            self.arrays["element_beam_ids"]     = np.array(self.raw_arrays["elNum1DA"])
            
         # Element shell ids are unchanged
        if "elNumA" in self.raw_arrays:   
            self.arrays["element_shell_ids"]    = np.array(self.raw_arrays["elNumA"]) 
            
         # Element solid ids are unchanged
        if "elNum3DA" in self.raw_arrays:  
            self.arrays["element_solid_ids"]    = np.array(self.raw_arrays["elNum3DA"])     
            
         # Element sph ids are unchanged
        if "nodNumSPH" in self.raw_arrays:  
            self.arrays["element_sph_ids"]      = np.array(self.raw_arrays["nodNumSPH"])             

        # Unpack 1D element part indexes

        if "defPart1DA" in self.raw_arrays:        
            if self.raw_arrays["defPart1DA"].any():
                start=0
                tmp_list_i  = np.empty(self.raw_arrays["defPart1DA"][-1], int)
                tmp_list_n  = np.empty(self.raw_arrays["defPart1DA"][-1], str)
                small_dict  = self.raw_arrays["pText1DA"]
                                                
                for _ipart1d, ipart1d in enumerate(self.raw_arrays["defPart1DA"]):
                    end                     =   ipart1d
                    num_el                  =   end + start 
                    _name                   =   small_dict[_ipart1d].strip()
                    _index                  =   _ipart1d
                    tmp_list_i[start:end]   =   _index
                    tmp_list_n[start:end]   =   _name
                    start                   =   end
    
                self.arrays["element_beam_part_indexes"]    = np.array(tmp_list_i)
                self.arrays["element_beam_part_ids"]        = np.array(tmp_list_n)                                                            

        # Unpack 2D element part indexes       

        if self.raw_header["nbParts"] > 0 and "defPartA" in self.raw_arrays:        
            if self.raw_arrays["defPartA"].any():
                start=0
                tmp_list_i  = np.empty(self.raw_arrays["defPartA"][-1], int)
                tmp_list_n  = np.empty(self.raw_arrays["defPartA"][-1], str)
                small_dict  = self.raw_arrays["pTextA"]
                     
                
                for _ipart2d, ipart2d in enumerate(self.raw_arrays["defPartA"]):
                    end                     =   ipart2d
                    num_el                  =   end + start 
                    _name                   =   small_dict[_ipart2d].strip()
                    _index                  =   _ipart2d
                    tmp_list_i[start:end]   =   _index
                    tmp_list_n[start:end]   =   _name
                    start                   =   end
    
                self.arrays["element_shell_part_indexes"]=np.array(tmp_list_i)
                self.arrays["element_shell_part_ids"]=np.array(tmp_list_n)                   
        
        # Unpack solid element part indexes
        if "defPart3DA" in self.raw_arrays:  
            if self.raw_arrays["defPart3DA"].any():
                start=0
                tmp_list_i  = np.empty(self.raw_arrays["defPart3DA"][-1], int)
                tmp_list_n  = np.empty(self.raw_arrays["defPart3DA"][-1], str)
                small_dict  = self.raw_arrays["pText3DA"]
                                                
                for _ipart3d, ipart3d in enumerate(self.raw_arrays["defPart3DA"]):
                    end                     =   ipart3d
                    num_el                  =   end + start 
                    _name                   =   small_dict[_ipart3d].strip()
                    _index                  =   _ipart3d
                    tmp_list_i[start:end]   =   _index
                    tmp_list_n[start:end]   =   _name
                    start                   =   end
    
                self.arrays["element_solid_part_indexes"]   = np.array(tmp_list_i)
                self.arrays["element_solid_part_ids"]       = np.array(tmp_list_n)               
        
        # Unpack sph element part indexes        
        if "defPartSPH" in self.raw_arrays:
            start=0
            tmp_list_i=[]
            tmp_list_n=[]
            for _ipart0d, ipart0d in enumerate(self.raw_arrays["defPartSPH"]):
                end             =   ipart0d
                num_el          =   end - start 
                _tmp_list_i     =   [_ipart0d] * num_el
                _tmp_list_n     =   [self.raw_arrays["pTextSPH"][_ipart0d].strip()] * num_el
                start           =   end
                tmp_list_i.extend(_tmp_list_i)
                tmp_list_n.extend(_tmp_list_n)
            self.arrays["element_sph_part_indexes"]     = np.array(tmp_list_i)
            self.arrays["element_sph_part_names"]       = np.array(tmp_list_n)                     
        
         # Eroded 1D element array is unchanged
        if "delElt1DA" in self.raw_arrays:
            self.arrays["element_beam_is_alive"] = np.array(self.raw_arrays["delElt1DA"])
        
         # Eroded shell element array is unchanged
        if "delEltA" in self.raw_arrays:        
            self.arrays["element_shell_is_alive"] = np.array(self.raw_arrays["delEltA"])
        
         # Eroded solid element array is unchanged
        if "delElt3DA" in self.raw_arrays:
            self.arrays["element_solid_is_alive"] = np.array(self.raw_arrays["delElt3DA"])                     
        
        # Unpack 1D scalar arrays
        if "nbEFunc1D" in self.raw_header:        
            for iefun in range(0, self.raw_header["nbEFunc1D"]):
                scalar              = "element_beam_" + str(self.raw_arrays["fText1DA"][iefun]).lower().replace(" ", "_").strip()
                start               = iefun * self.raw_header["nbElts1D"]
                end                 = (iefun+1) * self.raw_header["nbElts1D"]
                tmp_list            = self.raw_arrays["eFunc1DA"][start : end]
                self.arrays[scalar] = np.array(tmp_list)                     
                
        # Unpack 1D tensor arrays 
        if "nbTors1d" in self.raw_header:        
            for ietens in range(0, self.raw_header["nbTors1d"]):
                tensor                  = "element_beam_" + str(self.raw_arrays["fText1DA"][ietens]).lower().replace(" ", "_").strip()
                start                   = ietens * self.raw_header["nbElts1D"]
                end                     = (ietens+1) * self.raw_header["nbElts1D"]
                _tmp_list               = np.array(self.raw_arrays["tensVal1DA"][start : end])       
                
                self.arrays[tensor] = _tmp_list                                                
                
        # Unpack shell scalar arrays

        if "nbEFunc" in self.raw_header:        
            for iefun in range(0, self.raw_header["nbEFunc"]):
                scalar              = "element_shell_" + str(self.raw_arrays["fTextA"][iefun + self.raw_header["nbFunc"]]).lower().replace(" ", "_").strip()
                start               = iefun * self.raw_header["nbFacets"]
                end                 = (iefun+1) * self.raw_header["nbFacets"]
                tmp_list            = self.raw_arrays["eFuncA"][start : end]
                self.arrays[scalar] = np.array(tmp_list)                                         

        # Unpack shell tensor arrays     
        tTextA                      = np.array(self.raw_arrays["tTextA"])
        nbFacets                    = int(self.raw_header["nbFacets"])
        
        if "nbTens" in self.raw_header:        
            for ietens in range(0, self.raw_header["nbTens"]):
                tensor                  = tTextA[ietens]
                start                   = ietens * nbFacets
                end                     = (ietens+1) * nbFacets
                _tmp_list               = np.array(self.raw_arrays["tensValA"][start : end])
                                   
                self.arrays[tensor] = _tmp_list
         
        # Unpack solid scalar arrays 
        if "nbEFunc3D" in self.raw_header:        
            for iefun in range(0, self.raw_header["nbEFunc3D"]):
                scalar              = "element_solid_" + str(self.raw_arrays["fText3DA"][iefun]).lower().replace(" ", "_").strip()
                start               = iefun * self.raw_header["nbElts3D"]
                end                 = (iefun+1) * self.raw_header["nbElts3D"]
                tmp_list            = self.raw_arrays["eFunc3DA"][start : end]
                self.arrays[scalar] = np.array(tmp_list)                                  

        # Unpack solid tensor arrays 
        if "nbTens3D" in self.raw_header:        
            for ietens in range(0, self.raw_header["nbTens3D"]):
                tensor                  = "element_solid_" + str(self.raw_arrays["tText3DA"][ietens]).lower().replace(" ", "_").strip()
                start                   = ietens * self.raw_header["nbElts3D"]
                end                     = (ietens+1) * self.raw_header["nbElts3D"]
                _tmp_list               = np.array(self.raw_arrays["tensVal3DA"][start : end])                     
                
                self.arrays[tensor] = _tmp_list        

        # Unpack sph scalar arrays         
        if self.raw_header["flag_a_7"]:      
            for iefun in range(0, self.raw_header["nbEFuncSPH"]):
                scalar              = "element_sph_" + str(self.raw_arrays["scalTextSPH"][iefun]).lower().replace(" ", "_").strip()
                start               = iefun * self.raw_header["nbEltsSPH"]
                end                 = (iefun+1) * self.raw_header["nbEltsSPH"]
                tmp_list            = self.raw_arrays["eFuncSPH"][start : end]
                self.arrays[scalar] = np.array(tmp_list)   
                
            print("SPH Tensors")
            print((time.time() - time_start))
            time_start = time.time()                 
                    
            # Unpack sph tensor arrays                       
            if "nbTens" in self.raw_header:        
                for ietens in range(0, self.raw_header["nbTensSPH"]):
                    tensor                  = "element_sph_" + str(self.raw_arrays["tensTextSPH"][ietens]).lower().replace(" ", "_").strip()
                    start                   = ietens * self.raw_header["nbEltsSPH"]
                    end                     = (ietens+1) * self.raw_header["nbEltsSPH"]
                    _tmp_list               = np.array(self.raw_arrays["nbEltsSPH"][start : end])                          
                    
                    self.arrays[tensor] = _tmp_list                      

        # Unpack subsets        
        if self.raw_header["flag_a_4"]: 
            for nbSubset in self.raw_arrays["nbSubsets"]:
                self.arrays[nbSubset.lower().replace(" ", "_")] = self.raw_arrays["nbSubsets"][nbSubset]
        """ Unsure how to unpack these """                

        # Names of materials
        self.arrays["material_text_a"]    =  self.raw_arrays["materialTextA"]
        
        # Material type
        self.arrays["material_type_a"]    =  self.raw_arrays["materialTypeA"]
        
        # Property name
        self.arrays["properties_text_a"]  =  self.raw_arrays["propertiesTextA"]
        
        # Property type array
        self.arrays["properties_type_a"]  =  self.raw_arrays["propertiesTypeA"]

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
            try:
                val = bb.read_text(position,dtype[1]) .replace("\x01", "").replace("\x00", "").strip() 
            except:
                val=""
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
