########################################################################

import numpy
import vtk
from pyquaternion import Quaternion

########################################################################

def rotateUGrid_w_axis(old_ugrid, angle, axis, verbose=True):

    if (verbose): print('*** rotatePData_w_axis ***')

    nb_points = old_ugrid.GetNumberOfPoints()

    new_points = vtk.vtkPoints()
    new_points.SetNumberOfPoints(nb_points)

    my_quaternion = Quaternion(axis=axis, angle=angle)

    old_point = numpy.array([0.]*3)

    for num_point in range(nb_points):
        old_ugrid.GetPoint(num_point, old_point)
        #print old_point

        new_point = my_quaternion.rotate(old_point)
        #print new_point

        new_points.InsertPoint(num_point, new_point)

    #new_ugrid = vtk.vtkUnstructuredGrid()
    old_ugrid.SetPoints(new_points)
    print(old_ugrid)
    #new_ugrid.SetCells(old_ugrid.GetCells())
    
    return old_ugrid 
 
