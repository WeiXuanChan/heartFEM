########################################################################

import sys
import math
import numpy
import vtk
from random import uniform as randuniform
from heartFEM.lcleeHeart.vtk_py.addLocalFiberOrientation            import *
from heartFEM.lcleeHeart.vtk_py.addLocalProlateSpheroidalDirections import *
from heartFEM.lcleeHeart.vtk_py.createFloatArray                    import *
from heartFEM.lcleeHeart.vtk_py.getABPointsFromBoundsAndCenter      import *
from heartFEM.lcleeHeart.vtk_py.readSTL                             import *
from heartFEM.lcleeHeart.vtk_py.readUGrid                           import *
from heartFEM.lcleeHeart.vtk_py.writeUGrid                          import *

########################################################################
def normalDistribution(x,mu,sigma):
    if sigma==0.:
        if x==mu:
            return 1.
        else:
            return 0.
    return numpy.exp(-((x-mu)/sigma)**2./2.)/(sigma* 2.5066282746310002)
def randomWithNearby(randFunc,x,mu,sigma,numtry=10):
    if isinstance(randFunc,(int,float)):
        return randFunc
    randVal=numpy.zeros(10)
    probVal=numpy.zeros(10)
    randVal[0]=randFunc()
    probVal[0]=normalDistribution(randVal[0],mu,sigma)
    for n in range(1,numtry):
        randVal[n]=randFunc()
        probVal[n]=normalDistribution(randVal[n],mu,sigma)+probVal[n-1]
    return randVal[numpy.argmax(randuniform(0,probVal[-1])<=probVal)]
def addLocalFiberOrientation(ugrid_wall,
                             fiber_angle_end,
                             fiber_angle_epi,
                             points_AB=None,
                             fiberSheetletAngle=0., #can be a function(coord=), include **kwargs in input
                             fiberSheetletWidth=0., #can be a class with __call__(coord=,sheetlet_angle=), include **kwargs in input and self.max=max_width
                             radialFiberAngle=0., #can be a function(coord=), include **kwargs in input
                             fiberLength=0., #can be a class with __call__(coord=,fiber_radial_angle=), include **kwargs in input and self.max=max_length
                             verbose=True):

    if (verbose): print ('*** addLocalFiberOrientation ***')

    if (points_AB == None):
        points_AB = getABPointsFromBoundsAndCenter(ugrid_wall, verbose)
    assert (points_AB.GetNumberOfPoints() >= 2), "\"points_AB\" must have at least two points. Aborting."
    point_A = numpy.array([0.]*3)
    point_B = numpy.array([0.]*3)
    points_AB.GetPoint(                              0, point_A)
    points_AB.GetPoint(points_AB.GetNumberOfPoints()-1, point_B)
    distAB  = point_B - point_A
    longitudinalLength=numpy.linalg.norm(distAB)
    eL = distAB/longitudinalLength

    if (verbose): print ("Computing local fiber orientation...")

    farray_norm_dist_end = ugrid_wall.GetCellData().GetArray("norm_dist_end")
    farray_norm_dist_epi = ugrid_wall.GetCellData().GetArray("norm_dist_epi")
    farray_eRR = ugrid_wall.GetCellData().GetArray("eRR")
    farray_eCC = ugrid_wall.GetCellData().GetArray("eCC")
    farray_eLL = ugrid_wall.GetCellData().GetArray("eLL")

    nb_cells = ugrid_wall.GetNumberOfCells()

    farray_fiber_angle = createFloatArray("fiber_angle", 1, nb_cells)
    farray_fiber_radial_angle=createFloatArray("fiber_radial_angle", 1, nb_cells)
    farray_radialPos=createFloatArray("fiber_radial_position", 1, nb_cells)
    farray_sheetlet_angle=createFloatArray("fiber_sheetlet_angle", 1, nb_cells)
    farray_sheeletPos=createFloatArray("fiber_sheetlet_position", 1, nb_cells)
    farray_fiber_shift_towards_endo=createFloatArray("fiber_shift_towards_endo", 1, nb_cells)
    farray_fiber_normDensity = createFloatArray("fiber_normDensity", 1, nb_cells)
    farray_eF = createFloatArray("fiber vectors", 3, nb_cells)
    farray_eS = createFloatArray("sheet vectors", 3, nb_cells)
    farray_eN = createFloatArray("sheet normal vectors", 3, nb_cells)
    
    
    # LCL hack to have homogeneous fibers near apex
    #bds = ugrid_wall.GetBounds()
    #center = vtk.vtkCellCenters()
    #center.SetInputData(ugrid_wall)
    #center.Update()
    ###############################################
    cell_centers=[]
    for num_cell in range(nb_cells):
        cell_centers.append(numpy.array(ugrid_wall.GetPoints().GetPoint(num_cell)))
    cell_centers=numpy.array(cell_centers)
        
    
    for num_cell in range(nb_cells):
        norm_dist_end = farray_norm_dist_end.GetTuple(num_cell)[0]
        norm_dist_epi = farray_norm_dist_epi.GetTuple(num_cell)[0]
        
        eRR = numpy.array(farray_eRR.GetTuple(num_cell))
        eCC = numpy.array(farray_eCC.GetTuple(num_cell))
        eLL = numpy.array(farray_eLL.GetTuple(num_cell))
        if fiberSheetletAngle==0. and radialFiberAngle==0.:
            fiber_angle_in_degrees = (1.-norm_dist_end) * fiber_angle_end + (1.-norm_dist_epi) * fiber_angle_epi
            farray_fiber_angle.InsertTuple(num_cell, [fiber_angle_in_degrees])
            
            fiber_angle_in_radians = math.pi*fiber_angle_in_degrees/180
            eF = math.cos(fiber_angle_in_radians) * eCC + math.sin(fiber_angle_in_radians) * eLL
            eS = eRR
            eN = numpy.cross(eF, eS)
            fiber_radial_angle=0.
            radialPos=0.
            sheetlet_angle=0.
            sheeletPos=0.
            sheeletShifttowardsEndo=0.
            fdensity=1.
        else:
            #get nearby
            if isinstance(fiberSheetletWidth,(int,float)):
                temp_fiberSheetletWidth=fiberSheetletWidth
            else:
                temp_fiberSheetletWidth=fiberSheetletWidth.max
            if isinstance(fiberLength,(int,float)):
                temp_fiberLength=fiberLength
            else:
                temp_fiberLength=fiberLength.max
            maxNearbyRadius=max(temp_fiberSheetletWidth,temp_fiberLength)
            if maxNearbyRadius>0 and fiberSheetletWidth!=0 and fiberLength!=0 and not(isinstance(fiberSheetletAngle,(int,float))):
                if num_cell>0:
                    temp_cell_centers=cell_centers[:num_cell].copy()
                    temp_cell_centers_Index=numpy.arange(num_cell)
                    celldist=numpy.linalg.norm(numpy.cross(temp_cell_centers-cell_centers[num_cell],eRR),axis=1)
                    temp_cell_centers=temp_cell_centers[celldist<maxNearbyRadius].copy()
                    if len(temp_cell_centers)>0:
                        temp_cell_centers_Index=temp_cell_centers_Index[celldist<maxNearbyRadius].copy()
                        
                        cir_dist=numpy.dot(temp_cell_centers-cell_centers[num_cell],eCC)
                        norm_cir_dist=np.zeros(len(temp_cell_centers_Index))
                        for n in range(len(temp_cell_centers_Index)):
                            temp_radialPos=farray_radialPos.GetTuple(temp_cell_centers_Index[n])[0]
                            temp_fiber_radial_angle=farray_fiber_radial_angle.GetTuple(temp_cell_centers_Index[n])[0]
                            if isinstance(fiberLength,(int,float)):
                                temp_fiberLength=fiberLength
                            else:
                                temp_fiberLength=fiberLength(coord=cell_centers_Index[n],fiber_radial_angle=temp_fiber_radial_angle)
                            if cir_dist[n]>0.:
                                norm_cir_dist[n]=cir_dist[n]/numpy.abs(0.5*temp_fiberLength*math.cos(math.pi*temp_fiber_radial_angle/180)+temp_radialPos)
                            elif cir_dist[n]<0.:
                                norm_cir_dist[n]=cir_dist[n]/numpy.abs(0.5*temp_fiberLength*math.cos(math.pi*temp_fiber_radial_angle/180)-temp_radialPos)
                        
                        longi_dist=numpy.dot(temp_cell_centers-cell_centers[num_cell],eLL)
                        norm_longi_dist=np.zeros(len(temp_cell_centers_Index))
                        for n in range(len(temp_cell_centers_Index)):
                            temp_sheeletPos=farray_sheeletPos.GetTuple(temp_cell_centers_Index[n])[0]
                            temp_sheetlet_angle=farray_sheetlet_angle.GetTuple(temp_cell_centers_Index[n])[0]
                            if isinstance(fiberSheetletWidth,(int,float)):
                                temp_fiberSheetletWidth=fiberSheetletWidth
                            else:
                                temp_fiberSheetletWidth=fiberSheetletWidth(coord=cell_centers_Index[n],sheetlet_angle=temp_sheetlet_angle)
                            if longi_dist[n]>0.:
                                norm_longi_dist[n]=longi_dist[n]/numpy.abs(0.5*temp_fiberSheetletWidth*math.cos(math.pi*temp_sheetlet_angle/180)+temp_sheeletPos)
                            elif longi_dist[n]<0.:
                                norm_longi_dist[n]=longi_dist[n]/numpy.abs(0.5*temp_fiberSheetletWidth*math.cos(math.pi*temp_sheetlet_angle/180)-temp_sheeletPos)
                            
                        norm_celldist=numpy.sqrt(norm_cir_dist**2.+norm_longi_dist**2.)
                        followrandset=randuniform(0,1)
                        temp_cell_centers=temp_cell_centers[norm_celldist<=followrandset]
                        temp_cell_centers_Index=temp_cell_centers_Index[norm_celldist<=followrandset]
                        norm_celldist=norm_celldist[norm_celldist<=followrandset]
                if len(temp_cell_centers)>0:
                    minInd=numpy.argmin(norm_celldist)
                    fiber_radial_angle=farray_fiber_radial_angle.GetTuple(temp_cell_centers_Index[minInd])[0]
                    ref_radialPos=farray_radialPos.GetTuple(temp_cell_centers_Index[minInd])[0]
                    radialPos=ref_radialPos-numpy.dot(temp_cell_centers[minInd]-cell_centers[num_cell],eCC)
                    sheetlet_angle=farray_sheetlet_angle.GetTuple(temp_cell_centers_Index[minInd])[0]
                    ref_sheeletPos=farray_sheeletPos.GetTuple(temp_cell_centers_Index[minInd])[0]
                    sheeletPos=ref_sheeletPos-numpy.dot(temp_cell_centers[minInd]-cell_centers[num_cell],eLL)
                    sheeletShifttowardsEndo=sheeletPos*math.tan(math.pi*sheetlet_angle/180)-radialPos*math.tan(math.pi*fiber_radial_angle/180)
                else:
                    if isinstance(radialFiberAngle,(int,float)):
                        fiber_radial_angle=radialFiberAngle
                    else:
                        fiber_radial_angle=radialFiberAngle(coord=cell_centers[num_cell])
                    if isinstance(fiberLength,(int,float)):
                        temp_fiberLength=fiberLength
                    else:
                        temp_fiberLength=fiberLength(coord=cell_centers[num_cell],fiber_radial_angle=temp_fiber_radial_angle)
                    temp_radialFiberCirHalfLength=0.5*temp_fiberLength*math.cos(math.pi*fiber_radial_angle/180)
                    radialPos=randuniform(-temp_radialFiberCirHalfLength,temp_radialFiberCirHalfLength)
                    if isinstance(fiberSheetletAngle,(int,float)):
                        sheetlet_angle=fiberSheetletAngle
                    else:
                        sheetlet_angle=fiberSheetletAngle(coord=cell_centers[num_cell])
                    if isinstance(fiberSheetletWidth,(int,float)):
                        temp_fiberSheetletWidth=fiberSheetletWidth
                    else:
                        temp_fiberSheetletWidth=fiberSheetletWidth(coord=cell_centers[num_cell],sheetlet_angle=temp_sheetlet_angle)
                    temp_sheetLongiHalfLength=0.5*temp_fiberSheetletWidth*math.cos(math.pi*sheetlet_angle/180)
                    sheeletPos=randuniform(-temp_sheetLongiHalfLength,temp_sheetLongiHalfLength)
                sheeletShifttowardsEndo=sheeletPos*math.tan(math.pi*sheetlet_angle/180)-radialPos*math.tan(math.pi*fiber_radial_angle/180)
            else:
                if isinstance(radialFiberAngle,(int,float)):
                    fiber_radial_angle=radialFiberAngle
                else:
                    fiber_radial_angle=radialFiberAngle(coord=cell_centers[num_cell])
                if isinstance(fiberLength,(int,float)):
                    temp_fiberLength=fiberLength
                else:
                    temp_fiberLength=fiberLength(coord=cell_centers[num_cell],fiber_radial_angle=temp_fiber_radial_angle)
                temp_radialFiberCirHalfLength=0.5*temp_fiberLength*math.cos(math.pi*fiber_radial_angle/180)
                radialPos=randuniform(-temp_radialFiberCirHalfLength,temp_radialFiberCirHalfLength)
                if isinstance(fiberSheetletAngle,(int,float)):
                    sheetlet_angle=fiberSheetletAngle
                else:
                    sheetlet_angle=fiberSheetletAngle(coord=cell_centers[num_cell])
                if isinstance(fiberSheetletWidth,(int,float)):
                    temp_fiberSheetletWidth=fiberSheetletWidth
                else:
                    temp_fiberSheetletWidth=fiberSheetletWidth(coord=cell_centers[num_cell],sheetlet_angle=temp_sheetlet_angle)
                temp_sheetLongiHalfLength=0.5*temp_fiberSheetletWidth*math.cos(math.pi*sheetlet_angle/180)
                sheeletPos=randuniform(-temp_sheetLongiHalfLength,temp_sheetLongiHalfLength)
                sheeletShifttowardsEndo=sheeletPos*math.tan(math.pi*sheetlet_angle/180)-radialPos*math.tan(math.pi*fiber_radial_angle/180)
            fiberRaidalAngle_in_radians=math.pi*fiber_radial_angle/180
            fiberSheetletAngle_in_radians=math.pi*sheetlet_angle/180
            adjusted_norm_dist_end=norm_dist_end+sheeletShifttowardsEndo
            adjusted_norm_dist_epi=norm_dist_epi-sheeletShifttowardsEndo
            fiber_angle_in_degrees = max(0,min(1,1.-adjusted_norm_dist_end)) * fiber_angle_end + max(0,min(1,1.-adjusted_norm_dist_epi)) * fiber_angle_epi
            if adjusted_norm_dist_end>1 or adjusted_norm_dist_end<0:
                fdensity=0.
            else:
                fdensity=1.
            
            farray_fiber_angle.InsertTuple(num_cell, [fiber_angle_in_degrees])
            

            fiber_angle_in_radians = math.pi*fiber_angle_in_degrees/180
            temp_eF = math.cos(fiber_angle_in_radians) * eCC + math.sin(fiber_angle_in_radians) * eLL
            temp_eS = eRR
            temp_eN = numpy.cross(temp_eF, temp_eS)
            
            eF=math.cos(fiberRaidalAngle_in_radians) * temp_eF + math.sin(fiberRaidalAngle_in_radians) * eRR
            temp_eS = math.cos(fiberRaidalAngle_in_radians) * eRR - math.sin(fiberRaidalAngle_in_radians) * temp_eF
            
            eN=math.cos(fiberSheetletAngle_in_radians) * temp_eN + math.sin(fiberSheetletAngle_in_radians) * temp_eS
            eS=math.cos(fiberSheetletAngle_in_radians) * temp_eS - math.sin(fiberSheetletAngle_in_radians) * temp_eN
            
        farray_eF.InsertTuple(num_cell, eF)
        farray_eS.InsertTuple(num_cell, eS)
        farray_eN.InsertTuple(num_cell, eN)
        farray_fiber_radial_angle.InsertTuple(num_cell, [fiber_radial_angle])
        farray_radialPos.InsertTuple(num_cell, [radialPos])
        farray_sheetlet_angle.InsertTuple(num_cell, [sheetlet_angle])
        farray_sheeletPos.InsertTuple(num_cell, [sheeletPos])
        farray_fiber_shift_towards_endo.InsertTuple(num_cell, [sheeletShifttowardsEndo])
        farray_fiber_normDensity.InsertTuple(num_cell, [fdensity])
        

    if (verbose): print ("Filling mesh...")

    ugrid_wall.GetCellData().AddArray(farray_fiber_angle)
    ugrid_wall.GetCellData().AddArray(farray_fiber_radial_angle)
    ugrid_wall.GetCellData().AddArray(farray_radialPos)
    ugrid_wall.GetCellData().AddArray(farray_sheetlet_angle)
    ugrid_wall.GetCellData().AddArray(farray_sheeletPos)
    ugrid_wall.GetCellData().AddArray(farray_fiber_shift_towards_endo)
    ugrid_wall.GetCellData().AddArray(farray_fiber_normDensity)
    ugrid_wall.GetCellData().AddArray(farray_eF)
    ugrid_wall.GetCellData().AddArray(farray_eS)
    ugrid_wall.GetCellData().AddArray(farray_eN)

if (__name__ == "__main__"):
    assert (len(sys.argv) in [4]), "Number of arguments must be 3. Aborting."
    basename = sys.argv[1]
    ugrid_wall = readUGrid(basename + "-Mesh.vtk")
    pdata_end = readSTL(basename + "-End.stl")
    pdata_epi = readSTL(basename + "-Epi.stl")
    angle_end = float(sys.argv[2])
    angle_epi = float(sys.argv[3])
    addLocalProlateSpheroidalDirections(ugrid_wall, pdata_end, pdata_epi)
    addLocalFiberOrientation(ugrid_wall, angle_end, angle_epi)
    writeUGrid(ugrid_wall, basename + "-Mesh.vtk")
