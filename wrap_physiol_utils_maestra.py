#!/usr/bin/env python
""" Use f2py to wrap some fortran 

*** NOTE if it complains on compile it is beacause
MAESTRA needs to be compiled with same fortran compiler e..g

FC=/opt/local/bin/gfortran-mp-4.4

"""

import os
import sys
import subprocess

__author__  = "Martin De Kauwe"
__version__ = "1.0 (01.07.2011)"
__email__   = "mdekauwe@gmail.com"


arg ="make clean; make"
subprocess.call(arg, shell=True)

arg = "f2py -m maestcom -h maestcom.pyf maestcom.f90 --overwrite-signature"
subprocess.call(arg, shell=True)
arg = "f2py --fcompiler=gnu95 --f90exec=/opt/local/bin/gfortran-mp-4.4 \
        -c maestcom.pyf maestcom.o"
subprocess.call(arg, shell=True)

arg = "f2py -m utils -h utils.pyf utils.f90 metcom.f90 --overwrite-signature"
subprocess.call(arg, shell=True)
arg = "f2py --fcompiler=gnu95 --f90exec=/opt/local/bin/gfortran-mp-4.4 \
        -c utils.pyf utils.o"
subprocess.call(arg, shell=True)

#arg = "f2py -m physiol -h physiol.pyf physiol.f90 maestcom.f90 utils.f90 --overwrite-signature"
#subprocess.call(arg, shell=True)
#arg = "f2py --fcompiler=gnu95 --f90exec=/opt/local/bin/gfortran-mp-4.4 \
#        -c physiol.pyf maestcom.o physiol.o maestcom.o utils.o"
#subprocess.call(arg, shell=True)


arg = "f2py -m physiol -h physiol.pyf physiol.pyf physiol.f90 maestcom.f90 utils.f90 --overwrite-signature"
subprocess.call(arg, shell=True)

arg = "f2py --fcompiler=gnu95 --f90exec=/opt/local/bin/gfortran-mp-4.4 -c physiol.pyf physiol.o maestcom.o utils.o"
subprocess.call(arg, shell=True)