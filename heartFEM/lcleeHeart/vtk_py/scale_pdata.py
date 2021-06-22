########################################################################

import sys
import numpy
import vtk
from heartFEM.lcleeHeart.vtk_py import *

########################################################################

def scale_pdata(pdata, scale):
		if isinstance(scale,(float,int)):
		    scale=(scale,scale,scale)
		if origin is not None:
		    pdata=transform_pdata(pdata, (-origin[0],-origin[1],-origin[2]), (0.,0.,0.))
		transform = vtk.vtkTransform()
		
		transform.Scale(scale)
		transform.Update()

		transformfilter = vtk.vtkTransformPolyDataFilter()
		if (vtk.vtkVersion.GetVTKMajorVersion() >= 6):
			transformfilter.SetInputData(pdata)
		else:
			transformfilter.SetInput(pdata)
			transformfilter.SetTransform(transform)
			transformfilter.Update()
		pdata=transformfilter.GetOutput()
		if origin is not None:
		    pdata=transform_pdata(pdata, (-origin[0],-origin[1],-origin[2]), (0.,0.,0.))
		return pdata





