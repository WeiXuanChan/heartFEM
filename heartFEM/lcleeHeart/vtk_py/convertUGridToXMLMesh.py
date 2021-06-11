########################################################################

import vtk
import dolfin as dolfin
import numpy as np
from vtk.util import numpy_support

########################################################################

def convertUGridToXMLMesh(ugrid):


	num_pts = ugrid.GetNumberOfPoints()
	num_cells =  ugrid.GetNumberOfCells()
	
	celltypes = numpy_support.vtk_to_numpy(ugrid.GetCellTypesArray())
	num_tetra = np.count_nonzero(celltypes == 10)

	print ("Number of points  = ", num_pts)
	print ("Number of tetra  = ", num_tetra)

	mesh = dolfin.Mesh()
	editor = dolfin.MeshEditor()
#	c_type = cell_type()
#	c_str = c_type.type2string(c_type.cell_type())
	editor.open(mesh,'tetrahedron',3, 3)#for python3
	#editor.open(mesh,4, 3, 3)  # top. and geom. dimension are both 2
# The following argument types                     are supported:
 #   1. (self: dolfin.cpp.mesh.MeshEditor, mesh: dolfin.cpp.mesh.Mesh, type: str,tdim: int, gdim: int, degree: int=1) -> None

	editor.init_vertices(num_pts)  # number of vertices
	editor.init_cells(num_tetra)     # number of cells

	for p in range(0, num_pts):
		pt = ugrid.GetPoints().GetPoint(p)
		editor.add_vertex(p, [pt[0], pt[1], pt[2]])



	cnt =  0
	for p in range(0, num_cells):
		pts = vtk.vtkIdList()
		ugrid.GetCellPoints(p, pts)
		if(pts.GetNumberOfIds() == 4):
			editor.add_cell(cnt, [pts.GetId(0),  pts.GetId(1), pts.GetId(2), pts.GetId(3)])
			cnt = cnt + 1
		
	editor.close()

	return mesh
	



