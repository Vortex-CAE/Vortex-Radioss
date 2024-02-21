![Vortex-Radioss](https://github.com/Vortex-CAE/Vortex-Radioss/blob/main/docs/vortex-logo.png)  
Vortex-Radioss
======================
<img align="right" src="https://github.com/Vortex-CAE/Vortex-Radioss/assets/144378169/e97457cb-2355-4eca-8849-fd4b81b21f4b" alt="My Image">

Vortex-radioss is an open-source Python CAE library that allows for the reading of Radioss Animation files, Radioss Time-History files and for the conversion of Radioss Animation files into LS-Dyna's D3plot format.


Prerequisites
------------


This library is an extention to the open-source lasso-python library which needs installing as a prerequisite.


    python -m pip install lasso-python

Alternatively the library can be found here: https://github.com/open-lasso-python/lasso-python/

Modules
------------

### Radioss Animation Reader
This module will read in the animation data from one Radioss Animation file (One timestep), as the raw arrays are compressed and difficult to use, we reccommend using the unpacked arrays as described below. 

    # Load the module
    from vortex_radioss.animtod3plot.RadiossReader import RadiossReader
    
    # Generate animation object
    rr = RadiossReader("fileA001")

    ### To display high level information about the arrays
    # Return header dictionary
    rr.raw_header
    
    # Print header dictionary keys
    print(rr.raw_header.keys())

    ### To display the arrays as they come out of Radioss
    # Return raw arrays dictionary
    rr.raw_arrays
    
    # Print raw arrays dictionary keys
    print(rr.raw_arrays.keys())

    ### To display the arrays once they are unpacked
    # Return unpacked dictionary arrays
    rr.arrays
    
    # Print unpacked arrays dictionary keys
    print(rr.arrays.keys())    

### Radioss Time-History Reader    
This module will read in the data from one Radioss Time-History file (Multiple timesteps).

    # Load the module
    from vortex_radioss.animtod3plot.RadiossTHReader import RadiossTHReader
    
    # Generate Time-History object
    th = RadiossReader("fileT01")

    ### To display high level information about the arrays
    # Return header dictionary
    th.raw_header
    
    # Print header dictionary keys
    print(th.raw_header.keys())

    ### To display the arrays once they are unpacked
    # Return unpacked dictionary arrays
    th.arrays
    
    # Print unpacked arrays dictionary keys
    print(th.arrays.keys())     

### Radioss Animation to LS-Dyna D3Plot
This module converts Radioss Animation Files to LS-Dyna's D3plot format. Tensor and scalar array conversions are currently limited to thickness, stress, energy density and plastic strain for shell elements.

       # Load the module
      from vortex_radioss.animtod3plot.Anim_to_D3plot import readAndConvert
      
      # Use the file stem e.g for "folder/fileA001", "folder/fileA002" would be:
      readAndConvert("folder/file")
      
      # The D3plot files "folder/file.d3plot*" will be generated

CONTACT
------------
To get in touch please feel free to contact us on LinkedIn

Company:   
https://www.linkedin.com/company/vortex-cae     
      
LICENSING
------------

We chose to open-source license this code-base under Mozilla Public License, version 2.0.
This license offers minimal restrictions on the user, with the condition that modifications or enhancements made to the code should be submitted to this repository.

For more information please see the license file and https://www.mozilla.org/en-US/MPL/2.0/FAQ/

WHO WE ARE
------------

Vortex CAE is a trading name of: 

Vortex Engineering Group Ltd,    
Windsor House,      
Troon Way Business Centre,    
Humberstone Lane,   
Thurmaston,     
Leicestershire,     
England,   
LE4 9HA  

Company number: 11004422
