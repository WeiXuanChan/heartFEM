########################################################################

import sys
import vtk
import os
import inspect
from heartFEM.lcleeHeart import vtk_py

########################################################################

def createLVmesh(casename, meshsize, epifilename, endofilename, verbose=True):

    if (verbose): print ('*** createLVmesh ***')

    cur_dir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))) 
    savePath=os.path.dirname(os.path.abspath(casename))
    LVgeofile = cur_dir + "/LV.geo"
    LVtempgeofile = savePath+"/LVtemp.geo"
    mshfilename = casename + ".msh"
    vtkfilename = casename + ".vtk"
    print('cur_dir',cur_dir)
    cmd = "cp " + LVgeofile + " " + LVtempgeofile
    os.system(cmd)
    cmd = "sed -i.bak s/'<<mesh_d>>'/'" + str(meshsize) + "'/g " + LVtempgeofile
    os.system(cmd)
    cmd = "sed -i.bak s/'<<Endofilename>>'/'" + endofilename.replace('/','\/') + "'/g " + LVtempgeofile
    os.system(cmd)
    cmd = "sed -i.bak s/'<<Epifilename>>'/'" + epifilename.replace('/','\/') + "'/g " + LVtempgeofile
    os.system(cmd)
    cmd = "gmsh -3 "+LVtempgeofile+" -o " + mshfilename
    os.system(cmd)
    cmd = "gmsh -3 "+LVtempgeofile+" -o " + vtkfilename
    os.system(cmd)
    cmd = "rm "+LVtempgeofile
    os.system(cmd)



 
