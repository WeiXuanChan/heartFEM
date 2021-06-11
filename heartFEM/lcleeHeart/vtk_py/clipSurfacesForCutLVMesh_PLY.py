########################################################################

import sys
import vtk

from heartFEM.lcleeHeart.vtk_py.mat_vec_tools import *
from heartFEM.lcleeHeart.vtk_py.readPLY       import *
from heartFEM.lcleeHeart.vtk_py.writeSTL      import *

########################################################################

def clipSurfacesForCutLVMesh_PLY(endo, epi, height, verbose=True):

    if (verbose): print ('*** clipSurfacesForCutLVMesh ***')

    plane = vtk.vtkPlane()
    plane.SetNormal(0,0,-1)
    plane.SetOrigin(0,0,height)

    clip = vtk.vtkClipPolyData()
    clip.SetClipFunction(plane)
    if (vtk.vtkVersion.GetVTKMajorVersion() >= 6):
        clip.SetInputData(endo)
    else:
        clip.SetInput(endo)
    clip.Update()
    clipped_endo = clip.GetOutput(0)

    clip = vtk.vtkClipPolyData()
    clip.SetClipFunction(plane)
    if (vtk.vtkVersion.GetVTKMajorVersion() >= 6):
        clip.SetInputData(epi)
    else:
        clip.SetInput(epi)
    clip.Update()
    clipped_epi = clip.GetOutput(0)

    return clipped_endo, clipped_epi

if (__name__ == "__main__"):
    assert (len(sys.argv) in [3,4]), 'Number of arguments must be 2 or 3.'
    if (len(sys.argv) == 3):
        endo_filename = sys.argv[1]
        epi_filename = sys.argv[1]
        clipped_endo_filename = sys.argv[1] + '_CutLV-Endo.ply'
        clipped_epi_filename = sys.argv[1] + '_CutLV-Epi.ply'
        height = float(sys.argv[2])
    elif (len(sys.argv) == 4):
        endo_filename = sys.argv[1]
        epi_filename = sys.argv[2]
        clipped_endo_filename = 'clipped_endo.stl'
        clipped_epi_filename = 'clipped_epi.stl'
        height = float(sys.argv[3])
    endo = readPLY(endo_filename)
    epi = readPLY(epi_filename)
    clipped_endo, clipped_epi = clipSurfacesForCutLVMesh(endo, epi, height)
    writeSTL(clipped_endo, clipped_endo_filename)
    writeSTL(clipped_epi, clipped_epi_filename)

