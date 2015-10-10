# hyperQube
Flat-sky Power spectrum code for polarization and temperature maps.
The code has non standart dependencies (flipper, flipperPol and speck).
In his current state the code use the standard EB formalism (no pure B modes).
Please drop me an email (thbtouis@gmail.com) if you want a standard params file for the code (global.dict).
To run it:
HQcutPatches.py global.dict
HQcomputeMCM.py global.dict
HQcompileSpectra.py global.dict
HQcomputeAnalyticCovariance.py global.dict

Most of the code is written in python, apart from a part of the mode coupling calculation written in fortran: fortran.f90  (and need to be compiled with f2py).
