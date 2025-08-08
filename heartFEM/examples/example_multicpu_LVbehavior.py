#set basic case parameters    
stl_filename='t0' #rename stl to 't0.stl'
casePath=""
meshsize=0.6 #smaller = more refinement
clip_ratio=0.95
Laxis=[0,0,1] # set the base plane normal vector in x y z
baseplane_stl=None # you can set this to casePath+'base_cut.stl' instead of adjusting the Laxis (set Laxis to [-1,-1,-1] if to flip the normal of the baseplane)
EDV=1
ESV=2 #set volume

unload_geometry=False
#more setting can be found below

import os
import numpy as np
import shutil
from os import path
import sys
from heartFEM import LVclosed, lcleeHeart, heartModelTranslator, heartParameters

import motionSegmentation.BsplineFourier as BsplineFourier
import trimesh
import vtk_py3
try:
    from dolfin import *
except:
    from fenics import *
import logging
logging.basicConfig(level=logging.INFO,filename=casePath+'FEA_test.log')
#set default parameters with <defaultAge> and chambers configuration with <heartModel>
hr=LVclosed(defaultParameters=None,defaultAge='fetal28',heartModel='lcleeHeart') #the integer after "fetal is the number of weeks, defaultAge can also be set to 'adult'

#adjust the default parameters to preference
hr.defaultParameters.setParameter("FEniCS solver type",0) # can be set to 1 for more delicate solver when dealing with nan
hr.defaultParameters.setParameter('path of folder with case',casePath)
hr.defaultParameters.setParameter('filename of stl',stl_filename)
hr.defaultParameters.setParameter('mesh element size factor',meshsize)
hr.defaultParameters.setParameter('ratio to clip mesh',clip_ratio)
hr.defaultParameters.setParameter('solver volume adjustment',0.05)# make this smaller is you are dealing with smaller units
hr.defaultParameters.setParameter('LV longitudinal axis superior X',Laxis[0]) # longitudinal axis X towards base
hr.defaultParameters.setParameter('LV longitudinal axis superior Y',Laxis[1]) # longitudinal axis Y
hr.defaultParameters.setParameter('LV longitudinal axis superior Z',Laxis[2])
hr.defaultParameters.setParameter('Planar basal stl file for alignment',baseplane_stl)
hr.defaultParameters.setParameter('LV geometrical fiber direction variation file',None)
hr.defaultParameters.setParameter('LV geometrical fiber delay variation file',None)
hr.defaultParameters.setParameter('LV geometrical overall stiffness constant variation file',None)
hr.defaultParameters.setParameter('LV geometrical fiber strength variation file',None)
hr.defaultParameters.setParameter('LV end systolic volume in mL',EDV)
hr.defaultParameters.setParameter('LV end diastolic volume in mL',EDV)

hr.defaultParameters.setParameter('strain energy density function coefficient in Pa',100.)
# mode set to getting LV behaviour table
hr.defaultRunMode='LVbehaviorRun'
hr.(editParameters=None,
    unloadGeo=unload_geometry,
    minESV={"LV":ESV*0.95}, # it is better to set slightly lower than ESV if you are planning to use windkessel  
    maxEDV={"LV":EDV*1.05},
    volstep=50, #discretise results calculation between minESV and maxEDV
    multicpu=1, #set number of cpu to use
    outputResultList=None, #set outputResultList=["displacement"] to get the displacement results
    behavior_only=True) # True if you do not want to run windkessel



