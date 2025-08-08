'''
################################################
MIT License
Copyright (c) 2021 L. C. Lee
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
################################################
File: __init__.py
Description: Create mesh and setup object for solving Heart dynamics using FEniCS
             Contains externally usable class
             Note: Requires manual removal of exception in FEniCS (newer versions) during initial run
History:
    Date    Programmer SAR# - Description
    ---------- ---------- ----------------------------
  Author: L.C. LEE                   08MAR2019           - Created
  Author: w.x.chan@gmail.com         08MAR2021           - v2.0.0
                                                            -collated and turn into heart object
  Author: w.x.chan@gmail.com         10APR2021           - v2.1.0
                                                            -added function changeFiberAngles
  Author: w.x.chan@gmail.com         25MAY2021           - v2.2.0
                                                            -added inverse to inverse the energy function for calculating unloaded geometry
  Author: w.x.chan@gmail.com         25MAY2021           - v2.3.0
                                                            -debug inverse to solve for Finv = deformation gradient - I 
  Author: w.x.chan@gmail.com         09JUN2021           - v3.0.0
                                                            -integrate into heartFEM
  Author: w.x.chan@gmail.com         08JUL2021           - v3.1.0
                                                            -added outOfplaneDeg for meshing
  Author: w.x.chan@gmail.com         14JUL2021           - v3.3.0
                                                            -removed outOfplaneDeg and added fiberSheetletAngle,fiberSheetletWidth,radialFiberAngle,fiberLength for meshing
                                                            -added fiberSheetletAngle,fiberSheetletWidth,radialFiberAngle,fiberLength for changeFiberAngles
  Author: w.x.chan@gmail.com         28JUL2021           - v3.4.3
                                                            - added strainEnergy in ulforms "StrainEnergyDensityFunction_Coef"
  Author: w.x.chan@gmail.com         05OCT2021           - v3.6.0
                                                            - added  "StrainEnergyDensityFunction_Coef"
                                                            - added trackphase and updatephase
  Author: w.x.chan@gmail.com         20OCT2021           - v3.6.1
                                                            - added  input of endo and epi angles as file with columns x,y,z,angle
  Author: w.x.chan@gmail.com         25Dec2021           - v4.0.0
                                                            - added addBiVfiber into vtk
  Author: w.x.chan@gmail.com         06Dec2022           - v4.1.0
                                                            - added activeFiberStrength
                                                            - added activeFiberDelay
  Author: w.x.chan@gmail.com         01Mar2023           - v4.2.0
                                                            - added cstrainfactor
  Author: w.x.chan@gmail.com         08Nov2023           - v4.3.2
                                                            - using abs: on Laxis for calculative raxis
  Author: w.x.chan@gmail.com         07Jun2024           - v4.4.0
                                                            - added geometry: AVsplit
'''
_version='4.4.0'

import sys

import os as os
#sys.path.append(os.path.dirname(os.path.realpath(__file__)))
#import shutil
try:
    from dolfin import *
except:
    from fenics import *
import numpy as np
from .forms_shavik import Forms , Forms_AV
from .simpleActiveMaterial import SimpleActiveMaterial as Active
from .nsolver import NSolver as NSolver 
from .addfiber_matid import *

import vtk
from vtk.util import numpy_support
#from . import vtk_py as vtk_py
import vtk_py3 as vtk_py
#from mpi4py import MPI as pyMPI
from .rotateUGrid_w_axis import rotateUGrid_w_axis
from pyquaternion import Quaternion
import shutil

def generateMesh(casePath,stlname,Laxis,endo_angle='',epi_angle='',gmsh_surf_id=None,fiberSheetletAngle=0.,fiberSheetletWidth=0.,radialFiberAngle=0.,fiberLength=0.,cstrainFactor=1.,activeFiberStrength=1.,activeFiberDelay=0.,fiberDirectionVec=None,clipratio=0.75,meshsize=0.6,meshname=None,saveSubFolder='',RV_vector=None,fixed_transform=None,split_avj_extend=None,split_avj_refinement=None,model=''):
    if meshname is None:
        meshname=stlname
    if split_avj_refinement is None:
        split_avj_refinement=[40.,200.,10.,300.]
    meshfilename = casePath+ '/'+stlname+ '.stl'
    pdata = vtk_py.readSTL(meshfilename)
    if fixed_transform is not None:
        reverse_temp=False
        try:
            if np.all(Laxis<-0.99):
                reverse_temp=True
        except:
            print('Unable to detect if reversal of fixed_transform is needed. Laxis',repr(Laxis))
        import trimesh
        fixed_transform = trimesh.load(fixed_transform)
        fixed_transform_center=np.mean(fixed_transform.vertices,axis=0)
        fixed_transform_normal=np.mean(fixed_transform.face_normals,axis=0)
        Laxis=fixed_transform_normal/np.linalg.norm(fixed_transform_normal)
        if reverse_temp:
            Laxis=Laxis*-1.
            print('Laxis reversed.')
    else:
        Laxis = Laxis / np.linalg.norm(Laxis)
        print(pdata)
    C = vtk_py.getcentroid(pdata)
    angle = np.arccos(Laxis[2])
    if abs(Laxis[2])>0.99:
        raxis = np.cross(Laxis, np.array([1,0,0]))
        raxis = raxis / np.linalg.norm(raxis)
    else:
        raxis = np.cross(Laxis, np.array([0,0,1]))
        raxis = raxis / np.linalg.norm(raxis)
            
    ztop = pdata.GetBounds()[5]
    C = [C[0], C[1], 0.98*ztop+0.02*pdata.GetBounds()[4]]
    print(pdata.GetBounds())
    rotatedpdata = vtk_py.rotatePData_w_axis(pdata, angle, raxis)
    
    my_quaternion = Quaternion(axis=raxis, angle=angle)
    transformation_matrix_R=np.zeros((4,4))
    transformation_matrix_R[3,3]=1.
    transformation_matrix_R[:3,0]=my_quaternion.rotate(np.array([1,0,0]))
    transformation_matrix_R[:3,1]=my_quaternion.rotate(np.array([0,1,0]))
    transformation_matrix_R[:3,2]=my_quaternion.rotate(np.array([0,0,1]))
    if RV_vector is not None:
        RV_vector2=my_quaternion.rotate(np.array(RV_vector))
        angle2 = -np.arcsin(RV_vector2[1]/np.linalg.norm(RV_vector2[:2]))
        if np.arccos(RV_vector2[1]/np.linalg.norm(RV_vector2[:2]))<0:
            angle2=np.pi-angle2
        rotatedpdata = vtk_py.rotatePData_w_axis(rotatedpdata, angle2, np.array([0,0,1]))
        my_quaternion2 = Quaternion(axis=np.array([0,0,1]), angle=angle2)
        transformation_matrix_R2=np.zeros((4,4))
        transformation_matrix_R2[3,3]=1.
        transformation_matrix_R2[:3,0]=my_quaternion2.rotate(np.array([1,0,0]))
        transformation_matrix_R2[:3,1]=my_quaternion2.rotate(np.array([0,1,0]))
        transformation_matrix_R2[:3,2]=my_quaternion2.rotate(np.array([0,0,1]))
        transformation_matrix_R=transformation_matrix_R2.dot(transformation_matrix_R)
    
    if not(isinstance(clipratio,float)):
        clipratio=(transformation_matrix_R.dot(clipratio.reshape((-1,1))).reshape(-1)[2]-rotatedpdata.GetBounds()[4])/(rotatedpdata.GetBounds()[5]-rotatedpdata.GetBounds()[4])
    print(rotatedpdata.GetBounds())
    vtk_py.writeSTL(rotatedpdata, casePath+saveSubFolder+"/TEST.stl")
    isinsideout = 1
    C = vtk_py.getcentroid(rotatedpdata)
    C = [C[0], C[1], rotatedpdata.GetBounds()[4]+(rotatedpdata.GetBounds()[5]-rotatedpdata.GetBounds()[4])*clipratio] #you have xmin, xmax, ymin,ymax,zmin,zmax (is just to draw a box to clip the geometry, then you can get the edge of boundary)
	
    
    if fixed_transform is None:
        print(C)
        clippedheart = vtk_py.clipheart(rotatedpdata, C, [0,0,1],isinsideout, True) 
    else:
        clippedheart=rotatedpdata
    vtk_py.writeSTL(clippedheart, casePath+saveSubFolder+"/TEST2.stl")
    if model != 'AVsplit':
        if RV_vector is not None:
            epi, endo2 , endo= vtk_py.splitDomainBetweenEndoAndEpi(clippedheart,RV=True)
        else:
            epi, endo= vtk_py.splitDomainBetweenEndoAndEpi(clippedheart)
    
        cleanepipdata = vtk.vtkCleanPolyData()
        if (vtk.vtkVersion.GetVTKMajorVersion() >= 6):
            cleanepipdata.SetInputData(epi)
        else:
            cleanepipdata.SetInput(epi)
        cleanepipdata.Update()
        cleanepi = cleanepipdata.GetOutput()	

        cleanendopdata = vtk.vtkCleanPolyData()
        if (vtk.vtkVersion.GetVTKMajorVersion() >= 6):
            cleanendopdata.SetInputData(endo)
        else:
            cleanendopdata.SetInput(endo)
        cleanendopdata.Update()
        cleanendo = cleanendopdata.GetOutput()	
    
    #points_AB = getABPointsFromBoundsAndCenter(cleanepi, verbose)
    
    endofile = casePath+saveSubFolder+"/endo_lv" +'.stl'
    epifile = casePath+saveSubFolder+"/epi_lv" +'.stl'
   
    if model != 'AVsplit' and (abs(cleanendo.GetBounds()[5] - cleanendo.GetBounds()[4]) > abs(cleanepi.GetBounds()[5] - cleanepi.GetBounds()[4])):
        temp = cleanendo
        cleanendo = cleanepi
        cleanepi = temp
    
    
    if model != 'AVsplit' and RV_vector is not None:
        cleanendo2pdata = vtk.vtkCleanPolyData()
        if (vtk.vtkVersion.GetVTKMajorVersion() >= 6):
            cleanendo2pdata.SetInputData(endo2)
        else:
            cleanendo2pdata.SetInput(endo2)
        cleanendo2pdata.Update()
        cleanendo2 = cleanendo2pdata.GetOutput()	
        
        if (abs(cleanendo2.GetBounds()[5] - cleanendo2.GetBounds()[4]) > abs(cleanepi.GetBounds()[5] - cleanepi.GetBounds()[4])):
            temp = cleanendo2
            cleanendo2 = cleanepi
            cleanepi = temp
        if (cleanendo2.GetBounds()[1] < cleanendo.GetBounds()[1]):
            temp = cleanendo
            cleanendo = cleanendo2
            cleanendo2 = temp
        #points_AB = getABPointsFromBoundsAndCenter(cleanepi, verbose)
        
        endo2file = casePath+saveSubFolder+"/endo_rv" +'.stl'
        vtk_py.writeSTL(cleanendo2, endo2file)

    if model != 'AVsplit':
        vtk_py.writeSTL(cleanepi, epifile)
        vtk_py.writeSTL(cleanendo, endofile)
    
    filename = casePath+saveSubFolder+'/New_mesh' 
    if model == 'AVsplit':
        shutil.copy(casePath+saveSubFolder+"/TEST2.stl",casePath+saveSubFolder+"/gmsh_use.stl")
        if gmsh_surf_id is None:
            vtk_py.createAVmesh(filename, meshsize, casePath+saveSubFolder+"/gmsh_use.stl",split_avj_refinement)
        else:
            vtk_py.createAVmesh(filename, meshsize, casePath+saveSubFolder+"/gmsh_use.stl",split_avj_refinement,endo_atr_id=gmsh_surf_id[0],endo_ven_id=gmsh_surf_id[1],epi_atr_id=gmsh_surf_id[2],epi_ven_id=gmsh_surf_id[3],avj_id=gmsh_surf_id[4])
        
    elif RV_vector is not None:
        vtk_py.create_BiVmesh(casePath+saveSubFolder+"/epi_lv.stl", casePath+saveSubFolder+"/endo_lv.stl",casePath+saveSubFolder+"/endo_rv.stl", casename=filename, meshsize=meshsize)
    else:
        vtk_py.createLVmesh(filename, meshsize, casePath+saveSubFolder+"/epi_lv.stl", casePath+saveSubFolder+"/endo_lv.stl")
   
   
    ugridfilename = filename + '.vtk'
    pdatafilename = filename + '.vtp'
    ugrid = vtk_py.readUGrid(ugridfilename)
   
    mesh = vtk_py.convertUGridToXMLMesh(ugrid)
   	
    print (mesh)
   
    #print (mesh)
    #comm2 = pyMPI.COMM_WORLD
   
    # Rotate mesh (LCL include)
    rugrid = vtk_py.rotateUGrid(ugrid, rx=360)#original is rx=180, so it is not working, so I change it to rx=360, it works. 
    rugrid = vtk_py.rotateUGrid(rugrid, sx=0.1, sy=0.1, sz=0.1) #add the scale, the code did not working, so I remove it it works 
    transformation_matrix_S=np.array([[0.1,0,0,0],[0,0.1,0,0],[0,0,0.1,0],[0,0,0,1.]])
    
    vtk_py.writeUGrid(rugrid,casePath+saveSubFolder+"/rtemp.vtk")
    if model == 'AVsplit':
        fenics_mesh_ref, fenics_facet_ref, fenics_edge_ref, endoAtr_ugrid, endoVen_ugrid, epiAtr_ugrid, epiVen_ugrid, inlet_ugrid, avj_ugrid, outlet_ugrid = vtk_py.extractFeNiCsBiVFacet_AVsplit(rugrid,savePath=casePath+saveSubFolder+'/',split_avj_extend=split_avj_extend)
    elif RV_vector is not None:
        fenics_mesh_ref, fenics_facet_ref, fenics_edge_ref = vtk_py.extractFeNiCsBiVFacet(rugrid,savePath=casePath+saveSubFolder+'/', geometry = "BiV")
    else:
        fenics_mesh_ref, fenics_facet_ref, fenics_edge_ref = vtk_py.extractFeNiCsBiVFacet(rugrid,savePath=casePath+saveSubFolder+'/', geometry = "LV")

    matid = MeshFunction('size_t',fenics_mesh_ref, 3, mesh.domains())
    
    transformation_matrix=transformation_matrix_S.dot(transformation_matrix_R)
    if fixed_transform is not None:
        center_to_shift=transformation_matrix_R.dot(np.array([[fixed_transform_center[0]],[fixed_transform_center[1]],[fixed_transform_center[2]],[1.]])).reshape(-1)
        print("DISPLAY center_to_shift",center_to_shift)
        transformation_matrix[0,3]=-center_to_shift[0]/10.
        transformation_matrix[1,3]=-center_to_shift[1]/10.
        transformation_matrix[2,3]=-center_to_shift[2]/10.
        
        ztrans = Expression((str(-center_to_shift[0]/10.), str(-center_to_shift[1]/10.), str(-center_to_shift[2]/10.)), degree = 1)
    else:
        ztop =  max(fenics_mesh_ref.coordinates()[:,2])
        center_to_shift=np.array(cleanepi.GetCenter())
        print("DISPLAY center_to_shift",center_to_shift)
        
        transformation_matrix[0,3]=-center_to_shift[0]/10.
        transformation_matrix[1,3]=-center_to_shift[1]/10.
        transformation_matrix[2,3]=-ztop
        center_to_shift[2]=ztop*10.
        
        ztrans = Expression((str(-center_to_shift[0]/10.), str(-center_to_shift[1]/10.), str(-center_to_shift[2]/10.)), degree = 1)
    
    np.savetxt(casePath+saveSubFolder+'/'+meshname+'_Tmatrix_fromSTLtoFEniCS.txt',transformation_matrix)
    invtransformation_matrix=transformation_matrix_R.T.dot(np.array([[10.,0,0,0],[0,10.,0,0],[0,0,10.,0],[0,0,0,1.]]))
    inv_translate=invtransformation_matrix.dot(transformation_matrix[:,3:])
    invtransformation_matrix[:3,3]=-inv_translate[:3,0]
    np.savetxt(casePath+saveSubFolder+'/'+meshname+'_Tmatrix_fromFEniCStoSTL.txt',invtransformation_matrix)
    
    
    #ztrans = Expression(("0.0", "0.0", str(-ztop)), degree = 1)
   
   
    ALE.move(fenics_mesh_ref,ztrans) #if (dolfin.dolfin_version() != '1.6.0'): fenics_mesh_ref.move(ztrans)
    #shift and scale rotatedpdata
    if model != 'AVsplit':
        scaledrotatedpdata=vtk_py.scale_pdata(rotatedpdata,0.1)
        shiftscaledrotatedpdata=vtk_py.transform_pdata(scaledrotatedpdata,[-center_to_shift[0]/10.,-center_to_shift[1]/10.,-ztop],[0.,0.,0.])
        vtk_py.writeSTL(shiftscaledrotatedpdata, casePath+saveSubFolder+"/UnclipedFinal.stl")
        scaledrotatedpdata=vtk_py.scale_pdata(clippedheart,0.1)
        shiftscaledrotatedpdata=vtk_py.transform_pdata(scaledrotatedpdata,[-center_to_shift[0]/10.,-center_to_shift[1]/10.,-ztop],[0.,0.,0.])
        vtk_py.writeSTL(shiftscaledrotatedpdata, casePath+saveSubFolder+"/ClipedFinal.stl")
    else:
        endoAtr_pdata=vtk_py.convertUGridtoPdata(endoAtr_ugrid)
        endoAtr_pdata=vtk_py.transform_pdata(endoAtr_pdata,[-center_to_shift[0]/10.,-center_to_shift[1]/10.,-center_to_shift[2]/10.],[0.,0.,0.])
        endoVen_pdata=vtk_py.convertUGridtoPdata(endoVen_ugrid)
        endoVen_pdata=vtk_py.transform_pdata(endoVen_pdata,[-center_to_shift[0]/10.,-center_to_shift[1]/10.,-center_to_shift[2]/10.],[0.,0.,0.])
        epiAtr_pdata=vtk_py.convertUGridtoPdata(epiAtr_ugrid)
        epiAtr_pdata=vtk_py.transform_pdata(epiAtr_pdata,[-center_to_shift[0]/10.,-center_to_shift[1]/10.,-center_to_shift[2]/10.],[0.,0.,0.])
        epiVen_pdata=vtk_py.convertUGridtoPdata(epiVen_ugrid)
        epiVen_pdata=vtk_py.transform_pdata(epiVen_pdata,[-center_to_shift[0]/10.,-center_to_shift[1]/10.,-center_to_shift[2]/10.],[0.,0.,0.])
        
        inlet_pdata=vtk_py.convertUGridtoPdata(inlet_ugrid)
        inlet_pdata=vtk_py.transform_pdata(inlet_pdata,[-center_to_shift[0]/10.,-center_to_shift[1]/10.,-center_to_shift[2]/10.],[0.,0.,0.])
        avj_pdata=vtk_py.convertUGridtoPdata(avj_ugrid)
        avj_pdata=vtk_py.transform_pdata(avj_pdata,[-center_to_shift[0]/10.,-center_to_shift[1]/10.,-center_to_shift[2]/10.],[0.,0.,0.])
        outlet_pdata=vtk_py.convertUGridtoPdata(outlet_ugrid)
        outlet_pdata=vtk_py.transform_pdata(outlet_pdata,[-center_to_shift[0]/10.,-center_to_shift[1]/10.,-center_to_shift[2]/10.],[0.,0.,0.])
        
    mesh = fenics_mesh_ref
   
    gdim = mesh.geometry().dim()
   
    outdir = casePath+saveSubFolder+"/"
   
    quad_deg = 4
    VQuadelem = VectorElement("Quadrature", mesh.ufl_cell(), degree=quad_deg, quad_scheme="default")
    VQuadelem._quad_scheme = 'default'
    fiberFS = FunctionSpace(mesh, VQuadelem)
    SQuadelem = FiniteElement("Quadrature", mesh.ufl_cell(), degree=quad_deg, quad_scheme="default")
    SQuadelem._quad_scheme = 'default'
    scalarFS = FunctionSpace(mesh, SQuadelem)
    isepiflip = True #False
    isendoflip = True #True
    casedir=casePath+saveSubFolder+"/"+meshname+"_";  
   
    matid_filename = casePath+saveSubFolder+'/'+meshname + "_matid.pvd"
    File(matid_filename) << matid
    
    File(casePath+saveSubFolder+"/facetboundaries"+".pvd") << fenics_facet_ref
    File(casePath+saveSubFolder+"/edgeboundaries"+".pvd") << fenics_edge_ref
    File(casePath+saveSubFolder+"/mesh" + ".pvd") << mesh
    File(casePath+saveSubFolder+"/matid" +".pvd") << matid
   	
    File(casePath+saveSubFolder+"/gen_"+meshname+"_facetboundaries"+".pvd") << fenics_facet_ref
    File(casePath+saveSubFolder+"/gen_"+meshname+"_edgeboundaries"+".pvd") << fenics_edge_ref
    File(casePath+saveSubFolder+"/gen_"+meshname+"_mesh" + ".pvd") << mesh
    File(casePath+saveSubFolder+"/gen_"+meshname+"_matid" +".pvd") << matid
    if model == 'AVsplit':
        ef, es, en, eC, eL, eR,cstrainFactor,activefiberStrength,activefiberDelay = vtk_py.addLVfiber_AVsplit(mesh, scalarFS,fiberFS, "lv", endo_angle, epi_angle, casedir, endoAtr_pdata,endoVen_pdata, epiAtr_pdata, epiVen_pdata,inlet_pdata, avj_pdata, outlet_pdata,fiberSheetletAngle=fiberSheetletAngle,fiberSheetletWidth=fiberSheetletWidth,radialFiberAngle=radialFiberAngle,fiberLength=fiberLength,cstrainFactor=cstrainFactor,activeFiberStrength=activeFiberStrength,activeFiberDelay=activeFiberDelay,fiberDirectionVec=fiberDirectionVec)
    elif RV_vector is not None:
        eC, eL, eR=vtk_py.addBiVfiber(mesh,fenics_facet_ref, 0, 0,fiberSheetletAngle=0, saveDir=casedir)
        ef, es, en=vtk_py.addBiVfiber(mesh,fenics_facet_ref, endo_angle, epi_angle,fiberSheetletAngle=fiberSheetletAngle, saveDir=casedir)
    else:
        ef, es, en, eC, eL, eR,cstrainFactor,activefiberStrength,activefiberDelay = vtk_py.addLVfiber(mesh, scalarFS,fiberFS, "lv", endo_angle, epi_angle, casedir,fiberSheetletAngle=fiberSheetletAngle,fiberSheetletWidth=fiberSheetletWidth,radialFiberAngle=radialFiberAngle,fiberLength=fiberLength,cstrainFactor=cstrainFactor,activeFiberStrength=activeFiberStrength,activeFiberDelay=activeFiberDelay,fiberDirectionVec=fiberDirectionVec)
    
    f = HDF5File(mesh.mpi_comm(), casePath+saveSubFolder+'/'+meshname+".hdf5", 'w')
    f.write(mesh, meshname)
    f.close()
   	
    f = HDF5File(mesh.mpi_comm(), casePath+saveSubFolder+'/'+meshname+".hdf5", 'a') 
    f.write(fenics_facet_ref, meshname+"/"+"facetboundaries") 
    f.write(fenics_edge_ref, meshname+"/"+"edgeboundaries") 
    f.write(matid, meshname+"/"+"matid") 
    f.write(ef, meshname+"/"+"eF") 
    f.write(es, meshname+"/"+"eS") 
    f.write(en, meshname+"/"+"eN")
    f.write(eC, meshname+"/"+"eC") 
    f.write(eL, meshname+"/"+"eL") 
    f.write(eR, meshname+"/"+"eR")
    if RV_vector is None:
        f.write(cstrainFactor, meshname+"/"+"cstrainfactor")
        f.write(activefiberStrength, meshname+"/"+"activefiberstrength")
        f.write(activefiberDelay, meshname+"/"+"activefiberdelay")
    f.close()
   
   
    
   
   
   	# Rotate back fiber.vtu to check
    fiberugrid = vtk_py.readXMLUGrid(casePath+saveSubFolder+"/mesh000000.vtu")
    transfiberugrid = vtk_py.translate_mesh(fiberugrid, [center_to_shift[0],center_to_shift[1],ztop])
    Laxis_temp = np.array([0, 0, 1])
    angle_temp = np.arccos(Laxis[2])
    raxis_temp = np.cross(Laxis_temp, Laxis)
    raxis_temp = raxis_temp / np.linalg.norm(raxis_temp)
    rotatedugrid = rotateUGrid_w_axis(transfiberugrid, angle_temp, raxis_temp)
    rotatedugrid= vtk_py.rotateUGrid(rotatedugrid, sx=10, sy=10, sz=10) #add the scale, the code did not working, so I remove it it works 
    vtk_py.writeXMLUGrid(rotatedugrid, casePath+saveSubFolder+"/rotated_mesh.vtu")
   
    #src_dir="casePath/fiber.vtu"
    #dst_dir="casePath/fiber"+saveaddstr+".vtu"
    #shutil.copy(src_dir,dst_dir)
   
    
    #comm2 = pyMPI.COMM_WORLD
    #meshname = filename
   
    #f = HDF5File(mesh.mpi_comm(),  meshname+".hdf5", 'w')
    #f.write(mesh, meshname)
    #f.close()
    
   
    return 0
def changeFiberAngles(casePath,meshname,endo_angle,epi_angle,fiberSheetletAngle=0.,fiberSheetletWidth=0.,radialFiberAngle=0.,fiberLength=0.,cstrainFactor=1.,activeFiberStrength=1.,activeFiberDelay=0.,RV_vector=None,fiberDirectionVec=None,saveSubFolder='',loadSubFolder='',model='',surfID=None):
    os.makedirs(casePath+saveSubFolder,exist_ok=True)
    mesh = Mesh()
    f = HDF5File(MPI.comm_world, casePath+loadSubFolder+'/'+meshname+".hdf5", 'r') 
    f.read(mesh, meshname, False)
    
    facetboundaries = MeshFunction("size_t", mesh, 2)#!! when changing fiber directions, need to extract endoAtr_pdata,endoVen_pdata, epiAtr_pdata, epiVen_pdata,inlet_pdata, avj_pdata, outlet_pdata
    f.read(facetboundaries, meshname+"/"+"facetboundaries")
    
    edgeboundaries = MeshFunction('size_t', mesh,1)
    f.read(edgeboundaries, meshname+"/"+"edgeboundaries") 
    
    f.close()
    if model=='AVsplit':
        facet_filepath=casePath+saveSubFolder+"/"+meshname+"_facetboundaries000000.vtu"
        facets_labels=vtk_py.readXMLUGrid(facet_filepath)
        print(facet_filepath)
        print('surfID',surfID)
        for p in range(1,9):
            #thresd = vtk.vtkPolyData()
            pts_vtkThreshold=vtk.vtkThreshold()
            pts_vtkThreshold.SetInputData(facets_labels)
            pts_vtkThreshold.ThresholdBetween(p-0.5,p+0.5)
            pts_vtkThreshold.SetInputArrayToProcess(0, 0, 0, 1, 'f')
            pts_vtkThreshold.Update()
            temp_ugrid=pts_vtkThreshold.GetOutput()
            temp_pdata=vtk_py.convertUGridtoPdata(temp_ugrid)
            if p==surfID[0]:
                epiAtr_ugrid=pts_vtkThreshold.GetOutput()
                epiAtr_pdata=vtk_py.convertUGridtoPdata(epiAtr_ugrid)
            elif p==surfID[1]:
                epiVen_ugrid=pts_vtkThreshold.GetOutput()
                epiVen_pdata=vtk_py.convertUGridtoPdata(epiVen_ugrid)
            elif p==surfID[2]:
                endoAtr_ugrid=pts_vtkThreshold.GetOutput()
                endoAtr_pdata=vtk_py.convertUGridtoPdata(endoAtr_ugrid)
            elif p==surfID[3]:
                endoVen_ugrid=pts_vtkThreshold.GetOutput()
                endoVen_pdata=vtk_py.convertUGridtoPdata(endoVen_ugrid)
            elif p==surfID[4]:
                inlet_ugrid=pts_vtkThreshold.GetOutput()
                inlet_pdata=vtk_py.convertUGridtoPdata(inlet_ugrid)
            elif p==surfID[5]:
                outlet_ugrid=pts_vtkThreshold.GetOutput()
                outlet_pdata=vtk_py.convertUGridtoPdata(outlet_ugrid)
            elif p==surfID[6]:
                avj_ugrid=pts_vtkThreshold.GetOutput()
                avj_pdata=vtk_py.convertUGridtoPdata(avj_ugrid)
    gdim = mesh.geometry().dim()
   
    quad_deg = 4
    VQuadelem = VectorElement("Quadrature", mesh.ufl_cell(), degree=quad_deg, quad_scheme="default")
    VQuadelem._quad_scheme = 'default'
    fiberFS = FunctionSpace(mesh, VQuadelem)
    SQuadelem = FiniteElement("Quadrature", mesh.ufl_cell(), degree=quad_deg, quad_scheme="default")
    SQuadelem._quad_scheme = 'default'
    scalarFS = FunctionSpace(mesh, SQuadelem)
    casedir=casePath+saveSubFolder+"/"+meshname+"_"
    matid = MeshFunction('size_t',mesh, 3, mesh.domains())
    matid_filename = casePath+saveSubFolder+'/'+meshname + "_matid.pvd"
    if os.path.isfile(matid_filename):
        os.remove(matid_filename)
    File(matid_filename) << matid
   	
    File(casePath+saveSubFolder+"/facetboundaries"+".pvd") << facetboundaries
    File(casePath+saveSubFolder+"/edgeboundaries"+".pvd") << edgeboundaries
    File(casePath+saveSubFolder+"/mesh" + ".pvd") << mesh
    File(casePath+saveSubFolder+"/matid" +".pvd") << matid
    if model == 'AVsplit':
        ef, es, en, eC, eL, eR,cstrainFactor,activefiberStrength,activefiberDelay = vtk_py.addLVfiber_AVsplit(mesh, scalarFS,fiberFS, "lv", endo_angle, epi_angle, casedir, endoAtr_pdata,endoVen_pdata, epiAtr_pdata, epiVen_pdata,inlet_pdata, avj_pdata, outlet_pdata,fiberSheetletAngle=fiberSheetletAngle,fiberSheetletWidth=fiberSheetletWidth,radialFiberAngle=radialFiberAngle,fiberLength=fiberLength,cstrainFactor=cstrainFactor,activeFiberStrength=activeFiberStrength,activeFiberDelay=activeFiberDelay,fiberDirectionVec=fiberDirectionVec)
    elif RV_vector is not None:
        eC, eL, eR=vtk_py.addBiVfiber(mesh,facetboundaries, 0, 0,fiberSheetletAngle=0, saveDir=casedir)
        ef, es, en=vtk_py.addBiVfiber(mesh,facetboundaries, endo_angle, epi_angle,fiberSheetletAngle=fiberSheetletAngle, saveDir=casedir)
    else:
        ef, es, en, eC, eL, eR,cstrainFactor,activefiberStrength,activefiberDelay = vtk_py.addLVfiber(mesh, scalarFS,fiberFS, "lv", endo_angle, epi_angle, casedir,fiberSheetletAngle=fiberSheetletAngle,fiberSheetletWidth=fiberSheetletWidth,radialFiberAngle=radialFiberAngle,fiberLength=fiberLength,cstrainFactor=cstrainFactor,activeFiberStrength=activeFiberStrength,activeFiberDelay=activeFiberDelay,fiberDirectionVec=fiberDirectionVec)
                                    
    if os.path.isfile(casePath+saveSubFolder+'/'+meshname+".hdf5"):
        os.remove(casePath+saveSubFolder+'/'+meshname+".hdf5")
    f = HDF5File(mesh.mpi_comm(), casePath+saveSubFolder+'/'+meshname+".hdf5", 'w')
    f.write(mesh, meshname)
    f.close()
   	
    f = HDF5File(mesh.mpi_comm(), casePath+saveSubFolder+'/'+meshname+".hdf5", 'a') 
    f.write(facetboundaries, meshname+"/"+"facetboundaries") 
    f.write(edgeboundaries, meshname+"/"+"edgeboundaries") 
    f.write(matid, meshname+"/"+"matid") 
    f.write(ef, meshname+"/"+"eF") 
    f.write(es, meshname+"/"+"eS") 
    f.write(en, meshname+"/"+"eN")
    f.write(eC, meshname+"/"+"eC") 
    f.write(eL, meshname+"/"+"eL") 
    f.write(eR, meshname+"/"+"eR")
    if RV_vector is None:
        f.write(cstrainFactor, meshname+"/"+"cstrainfactor")
        f.write(activefiberStrength, meshname+"/"+"activefiberstrength")
        f.write(activefiberDelay, meshname+"/"+"activefiberdelay")
   
    f.close()
   
   
   
    return 0
def print_hdf5_para(read_path,meshname=None,print_path=None):
    #from fenics import *
    mesh = Mesh()
    if meshname is None:
        meshname='t0_unloadedmesh'
    f = HDF5File(MPI.comm_world, read_path, 'r') #"/media/yap/drive3/WeiXuan/zebrafish/Renee_1/fish_test3/t0_unloadedmesh_w_fiberstrength.hdf5"
    f.read(mesh, meshname, False)
    
    facetboundaries = MeshFunction("size_t", mesh, 2)
    f.read(facetboundaries, meshname+"/"+"facetboundaries")
    edgeboundaries = MeshFunction('size_t', mesh,1)
    f.read(edgeboundaries, meshname+"/"+"edgeboundaries") 
    fiberFS = VectorFunctionSpace(mesh, 'DG', 0)
    VQuadelem = VectorElement("Quadrature", mesh.ufl_cell(), degree=4, quad_scheme="default") 
    VQuadelem._quad_scheme = 'default'
    for e in VQuadelem.sub_elements():
    	e._quad_scheme = 'default'
    if f.has_dataset(meshname+"/"+"matid"):
        matid = MeshFunction('size_t',mesh, 3, mesh.domains())
        f.read(matid, meshname+"/"+"matid")
    else:
        print('matid does not exist in '+read_path)
        matid=None
    fiberFS = FunctionSpace(mesh, VQuadelem)
    f0 = Function(fiberFS)
    s0 = Function(fiberFS)
    n0 = Function(fiberFS)
    
    f.read(f0, meshname+"/"+"eF")
    f.read(s0, meshname+"/"+"eS")
    f.read(n0, meshname+"/"+"eN")
    
    long0 = Function(fiberFS)
    cir0 = Function(fiberFS)
    rad0 = Function(fiberFS)
    
    f.read(long0, meshname+"/"+"eL")
    f.read(cir0, meshname+"/"+"eC")
    f.read(rad0, meshname+"/"+"eR")
    if print_path is not None:
        os.makedirs(print_path,exist_ok=True)
    if f.has_dataset(meshname+"/"+"cstrainfactor"):
        cstrainFactor = Function(fiberFS)
        f.read(cstrainFactor, meshname+"/"+"cstrainfactor")
        print("cstrainFactor")
        print(cstrainFactor.sub(0).vector()[:].min(),' to ',cstrainFactor.sub(0).vector()[:].max())
        if print_path is not None:
            file = fenics.File(print_path+"/cstrainFactor.pvd")
            result=fenics.project(cstrainFactor,fenics.VectorFunctionSpace(mesh, "DG",0))
            result.rename("cstrainFactor","cstrainFactor")
            file << result
    if f.has_dataset(meshname+"/"+"activefiberstrength"):
        activefiberstrength = Function(fiberFS)
        f.read(activefiberstrength, meshname+"/"+"activefiberstrength")
        print("activefiberstrength")
        print(activefiberstrength.sub(0).vector()[:].min(),' to ',activefiberstrength.sub(0).vector()[:].max())
        if print_path is not None:
            file = fenics.File(print_path+"/activefiberstrength.pvd")
            result=fenics.project(activefiberstrength,fenics.VectorFunctionSpace(mesh, "DG", 0))
            result.rename("activefiberstrength","activefiberstrength")
            file << result
    if f.has_dataset(meshname+"/"+"activefiberdelay"):
        activefiberdelay = Function(fiberFS)
        f.read(activefiberdelay, meshname+"/"+"activefiberdelay")
        print("activefiberdelay")
        print(activefiberdelay.sub(0).vector()[:].min(),' to ',activefiberdelay.sub(0).vector()[:].max())
        if print_path is not None:
            file = fenics.File(print_path+"/activefiberdelay.pvd")
            result=fenics.project(activefiberdelay,fenics.VectorFunctionSpace(mesh, "DG", 0))
            result.rename("activefiberdelay","activefiberdelay")
            file << result
    f.close()
    
    
def adjust_hdf5_para(read_path,save_path,meshname=None,new_meshname=None,cstrainfactor_assign=None,activefiberstrength_assign=None,activefiberdelay_assign=None):
    #from fenics import *
    mesh = Mesh()
    if meshname is None:
        meshname='t0_unloadedmesh'
    if new_meshname is None:
        new_meshname=meshname
    f = HDF5File(MPI.comm_world, read_path, 'r') #"/media/yap/drive3/WeiXuan/zebrafish/Renee_1/fish_test3/t0_unloadedmesh_w_fiberstrength.hdf5"
    f.read(mesh, meshname, False)
    
    facetboundaries = MeshFunction("size_t", mesh, 2)
    f.read(facetboundaries, meshname+"/"+"facetboundaries")
    edgeboundaries = MeshFunction('size_t', mesh,1)
    f.read(edgeboundaries, meshname+"/"+"edgeboundaries") 
    fiberFS = VectorFunctionSpace(mesh, 'DG', 0)
    VQuadelem = VectorElement("Quadrature", mesh.ufl_cell(), degree=4, quad_scheme="default") 
    VQuadelem._quad_scheme = 'default'
    for e in VQuadelem.sub_elements():
    	e._quad_scheme = 'default'
    if f.has_dataset(meshname+"/"+"matid"):
        matid = MeshFunction('size_t',mesh, 3, mesh.domains())
        f.read(matid, meshname+"/"+"matid")
    else:
        print('matid does not exist in '+read_path)
        matid=None
    fiberFS = FunctionSpace(mesh, VQuadelem)
    f0 = Function(fiberFS)
    s0 = Function(fiberFS)
    n0 = Function(fiberFS)
    
    f.read(f0, meshname+"/"+"eF")
    f.read(s0, meshname+"/"+"eS")
    f.read(n0, meshname+"/"+"eN")
    
    long0 = Function(fiberFS)
    cir0 = Function(fiberFS)
    rad0 = Function(fiberFS)
    
    f.read(long0, meshname+"/"+"eL")
    f.read(cir0, meshname+"/"+"eC")
    f.read(rad0, meshname+"/"+"eR")
    
    if f.has_dataset(meshname+"/"+"cstrainfactor"):
        cstrainFactor = Function(fiberFS)
        f.read(cstrainFactor, meshname+"/"+"cstrainfactor")
    else:
        print('cstrainFactor does not exist in '+read_path)
        if cstrainfactor_assign is None:
            cstrainfactor_assign=False
    if f.has_dataset(meshname+"/"+"activefiberstrength"):
        activefiberstrength = Function(fiberFS)
        f.read(activefiberstrength, meshname+"/"+"activefiberstrength")
    else:
        print('activefiberstrength does not exist in '+read_path)
        if activefiberstrength_assign is None:
            activefiberstrength_assign=False
    if f.has_dataset(meshname+"/"+"activefiberdelay"):
        activefiberdelay = Function(fiberFS)
        f.read(activefiberdelay, meshname+"/"+"activefiberdelay")
    else:
        print('activefiberdelay does not exist in '+read_path)
        if activefiberdelay_assign is None:
            activefiberdelay_assign=False
    f.close()
    
    f = HDF5File(MPI.comm_world, save_path, 'w')#"/media/yap/drive3/WeiXuan/zebrafish/Renee_1/fish_test3/t0_unloadedmesh.hdf5"
    f.write(mesh, new_meshname)
    f.close()
    
    if cstrainfactor_assign is not None:
        if cstrainfactor_assign is not False:
            if isinstance(cstrainfactor_assign,(int,float,list,np.ndarray)):
                cstrainFactor.sub(0).vector()[:]=cstrainfactor_assign
                cstrainFactor.sub(1).vector()[:]=cstrainfactor_assign
                cstrainFactor.sub(2).vector()[:]=cstrainfactor_assign
            else:
                temp=np.stack((cstrainFactor.sub(0).vector()[:],cstrainFactor.sub(1).vector()[:],cstrainFactor.sub(2).vector()[:]),axis=0)
                cstrainFactor.sub(0).vector()[:]=cstrainfactor_assign(temp,axis=0)
                cstrainFactor.sub(1).vector()[:]=cstrainfactor_assign(temp,axis=1)
                cstrainFactor.sub(2).vector()[:]=cstrainfactor_assign(temp,axis=2)
    if activefiberstrength_assign is not None:
        if activefiberstrength_assign is not False:
            if isinstance(activefiberstrength_assign,(int,float,list,np.ndarray)):
                activefiberstrength.sub(0).vector()[:]=activefiberstrength_assign
                activefiberstrength.sub(1).vector()[:]=activefiberstrength_assign
                activefiberstrength.sub(2).vector()[:]=activefiberstrength_assign
            else:
                temp=np.stack((activefiberstrength.sub(0).vector()[:],activefiberstrength.sub(1).vector()[:],activefiberstrength.sub(2).vector()[:]),axis=0)
                activefiberstrength.sub(0).vector()[:]=activefiberstrength_assign(temp,axis=0)
                activefiberstrength.sub(1).vector()[:]=activefiberstrength_assign(temp,axis=1)
                activefiberstrength.sub(2).vector()[:]=activefiberstrength_assign(temp,axis=2)
    if activefiberdelay_assign is not None:
        if activefiberdelay_assign is not False:
            if isinstance(activefiberdelay_assign,(int,float,list,np.ndarray)):
                activefiberdelay.sub(0).vector()[:]=activefiberdelay_assign
                activefiberdelay.sub(1).vector()[:]=activefiberdelay_assign
                activefiberdelay.sub(2).vector()[:]=activefiberdelay_assign
            else:
                temp=np.stack((activefiberdelay.sub(0).vector()[:],activefiberdelay.sub(1).vector()[:],activefiberdelay.sub(2).vector()[:]),axis=0)
                activefiberdelay.sub(0).vector()[:]=activefiberdelay_assign(temp,axis=0)
                activefiberdelay.sub(1).vector()[:]=activefiberdelay_assign(temp,axis=1)
                activefiberdelay.sub(2).vector()[:]=activefiberdelay_assign(temp,axis=2)
        
    f = HDF5File(MPI.comm_world, save_path, 'a') 
    f.write(facetboundaries, new_meshname+"/"+"facetboundaries") 
    f.write(edgeboundaries, new_meshname+"/"+"edgeboundaries") 
    if matid is not None:
        f.write(matid, new_meshname+"/"+"matid") 
    f.write(f0, new_meshname+"/"+"eF") 
    f.write(s0, new_meshname+"/"+"eS") 
    f.write(n0, new_meshname+"/"+"eN")
    f.write(cir0, new_meshname+"/"+"eC") 
    f.write(long0, new_meshname+"/"+"eL") 
    f.write(rad0, new_meshname+"/"+"eR")
    if cstrainfactor_assign is not False:
        print("writing cstrainfactor")
        f.write(cstrainFactor, new_meshname+"/"+"cstrainfactor")
    if activefiberstrength_assign is not False:
        print("writing activefiberstrength")
        f.write(activefiberstrength, new_meshname+"/"+"activefiberstrength")
    if activefiberdelay_assign is not False:
        print("writing activefiberdelay")
        f.write(activefiberdelay, new_meshname+"/"+"activefiberdelay")
    
    f.close()
    
def transfer_vtu_data_to_mesh(read_data,mesh,save_file,clean=True):
    from scipy.spatial.distance import cdist
    
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(read_data)
    reader.Update()
    reader_data = reader.GetOutput()
    
    mesher = vtk.vtkXMLUnstructuredGridReader()
    mesher.SetFileName(mesh)
    mesher.Update()
    mesh_data = mesher.GetOutput()
    
    meshpointcoord=numpy_support.vtk_to_numpy(mesh_data.GetPoints().GetData())
    readpointcoord=numpy_support.vtk_to_numpy(reader_data.GetPoints().GetData())
    
    centroider=vtk.vtkCellCenters()
    centroider.SetInputData(reader_data)
    centroider.Update()
    centroid_data=centroider.GetOutput()
    
    readcellcoord=numpy_support.vtk_to_numpy(centroid_data.GetPoints().GetData())
    if clean:
        while mesh_data.GetPointData().GetArrayName(0) is not None:
            mesh_data.GetPointData().RemoveArray(mesh_data.GetPointData().GetArrayName(0))
        while mesh_data.GetCellData().GetArrayName(0) is not None:
            mesh_data.GetCellData().RemoveArray(mesh_data.GetCellData().GetArrayName(0))
    for pointdataN in range(reader_data.GetPointData().GetNumberOfArrays()):
        if pointdataN==0:
            nearest_data_assignment=cdist(meshpointcoord,readpointcoord)
            nearest_data_assignment=np.argmin(nearest_data_assignment,axis=-1)
        temp_name=reader_data.GetPointData().GetArrayName(pointdataN)
        print("Copying point data",temp_name)
        temp_array=numpy_support.vtk_to_numpy(reader_data.GetPointData().GetArray(temp_name))
        temp_field = numpy_support.numpy_to_vtk(temp_array[nearest_data_assignment])
        temp_field.SetName(temp_name)
        mesh_data.GetPointData().AddArray(temp_field)
        mesher.Update()
        
    for celldataN in range(reader_data.GetCellData().GetNumberOfArrays()):
        if celldataN==0:
            nearest_data_assignment=cdist(meshpointcoord,readcellcoord)
            nearest_data_assignment=np.argmin(nearest_data_assignment,axis=-1)
        temp_name=reader_data.GetCellData().GetArrayName(celldataN)
        print("Copying cell data",temp_name)
        temp_array=numpy_support.vtk_to_numpy(reader_data.GetCellData().GetArray(temp_name))
        temp_field = numpy_support.numpy_to_vtk(temp_array[nearest_data_assignment])
        temp_field.SetName(temp_name)
        mesh_data.GetPointData().AddArray(temp_field)
        mesher.Update()
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetInputData(mesh_data)
    writer.SetFileName(save_file)
    writer.Update()
    return 0
    
class heart:
    def __new__(cls, *args, **kwargs):
        return super().__new__(cls)
    def __init__(self,casename,meshname,runParameters,inverseHeart=False,trackphase=False):
        self.debug_call=[]
        self.inverseHeart=inverseHeart
        self.trackphase=trackphase
        self.backup=(casename,meshname,dict(runParameters))
        runParameters=dict(runParameters)
        parameters["form_compiler"]["quadrature_degree"]=4
        parameters["form_compiler"]["representation"] = "uflacs"
        parameters['ghost_mode'] = 'shared_facet'
        
        self.mesh = Mesh()
        f = HDF5File(MPI.comm_world, casename+'/'+meshname+".hdf5", 'r') 
        f.read(self.mesh, meshname, False)
        File(casename+"/facetboundaries.pvd.pvd")<<self.mesh
        
        self.facetboundaries = MeshFunction("size_t", self.mesh, 2)
        f.read(self.facetboundaries, meshname+"/"+"facetboundaries")
        
        #f.read(self.facetboundaries, meshname+"/"+"facetboundaries")
        self.ds = dolfin.ds(subdomain_data = self.facetboundaries)
        File(casename+"/"+"facetboundaries.pvd") << self.facetboundaries
        
        self.fiberFS = VectorFunctionSpace(self.mesh, 'DG', 0)
        VQuadelem = VectorElement("Quadrature", self.mesh.ufl_cell(), degree=4, quad_scheme="default") #segmentation fault happened due to degree = 2. degree = 4 should be implied. I changed it. (Joy) 
        VQuadelem._quad_scheme = 'default'
        for e in VQuadelem.sub_elements():
        	e._quad_scheme = 'default'
        SQuadelem = FiniteElement("Quadrature", self.mesh.ufl_cell(), degree=4, quad_scheme="default")
        SQuadelem._quad_scheme = 'default'
        scalarFS = FunctionSpace(mesh, SQuadelem)
        
        self.fiberFS = FunctionSpace(self.mesh, VQuadelem)
        self.f0 = Function(self.fiberFS)
        self.s0 = Function(self.fiberFS)
        self.n0 = Function(self.fiberFS)
        
        f.read(self.f0, meshname+"/"+"eF")
        f.read(self.s0, meshname+"/"+"eS")
        f.read(self.n0, meshname+"/"+"eN")
        
        self.long0 = Function(self.fiberFS)
        self.cir0 = Function(self.fiberFS)
        self.rad0 = Function(self.fiberFS)
        
        f.read(self.long0, meshname+"/"+"eL")
        f.read(self.cir0, meshname+"/"+"eC")
        f.read(self.rad0, meshname+"/"+"eR")
        
        try:
            self.cstrainFactor = Function(self.fiberFS)
            f.read(self.cstrainFactor, meshname+"/"+"cstrainfactor")
            print("cstrainFactor read.")
            print("cstrainFactor min-max:", self.cstrainFactor.sub(0).vector()[:].min(), self.cstrainFactor.sub(0).vector()[:].max())
        except Exception as e:
            self.cstrainFactor =None
            print("cstrainFactor NOT read.")
            print(repr(e))
        try:
            self.activefiberstrength = Function(self.fiberFS)
            f.read(self.activefiberstrength, meshname+"/"+"activefiberstrength")
            print("activefiberstrength read.")
            print("activefiberstrength min-max:", self.activefiberstrength.sub(0).vector()[:].min(), self.activefiberstrength.sub(0).vector()[:].max())
        except:
            self.activefiberstrength =None
            print("activefiberstrength NOT read.")
        try:
            self.activefiberdelay = Function(self.fiberFS)
            f.read(self.activefiberdelay, meshname+"/"+"activefiberdelay")
            print("activefiberdelay read.")
            print("activefiberdelay min-max:", self.activefiberdelay.sub(0).vector()[:].min(), self.activefiberdelay.sub(0).vector()[:].max())
        except:
            self.activefiberdelay = None
            print("activefiberdelay NOT read.")
        
        f.close()
        
        topid = runParameters['topid']
        endoid = runParameters['endoid']
        epiid = runParameters['epiid']#####################################################################
        
        self.comm = self.mesh.mpi_comm()
        
        
        N = FacetNormal (self.mesh)
        X = SpatialCoordinate (self.mesh)
        Kspring = Constant(runParameters['Kspring_constant'])
        self.Press = Expression(("P"), P=0.0, degree=2)
        self.Cavityvol = Expression(("vol"), vol=0.0, degree=2)
        
        
        V = VectorFunctionSpace(self.mesh, 'CG', 2)
        TF = TensorFunctionSpace(self.mesh, 'DG', 1)
        self.Q = FunctionSpace(self.mesh,'CG',1)
        R = FunctionSpace(self.mesh,'Real',0)
        
        Velem = VectorElement("CG", self.mesh.ufl_cell(), 2, quad_scheme="default")
        Velem._quad_scheme = 'default'
        for e in Velem.sub_elements():
        	e._quad_scheme = 'default'
        Qelem = FiniteElement("CG", self.mesh.ufl_cell(), 1, quad_scheme="default")
        Qelem._quad_scheme = 'default'
        for e in Qelem.sub_elements():
        	e._quad_scheme = 'default'
        Relem = FiniteElement("Real", self.mesh.ufl_cell(), 0, quad_scheme="default")
        Relem._quad_scheme = 'default'
        for e in Relem.sub_elements():
        	e._quad_scheme = 'default'
        Quadelem = FiniteElement("Quadrature", self.mesh.ufl_cell(), degree=4, quad_scheme="default")
        Quadelem._quad_scheme = 'default'
        for e in Quadelem.sub_elements():
        	e._quad_scheme = 'default'
        
        Telem2 = TensorElement("Quadrature", self.mesh.ufl_cell(), degree=4, shape=2*(3,), quad_scheme='default')
        Telem2._quad_scheme = 'default'
        for e in Telem2.sub_elements():
        	e._quad_scheme = 'default'
        Telem3 = TensorElement("Lagrange", self.mesh.ufl_cell(), degree=2, shape=(3,3), quad_scheme='default')
        Telem3._quad_scheme = 'default'
        for e in Telem3.sub_elements():
        	e._quad_scheme = 'default'
        Telem4 = TensorElement("Quadrature", self.mesh.ufl_cell(), degree=4, shape=4*(3,), quad_scheme='default')
        Telem4._quad_scheme = 'default'
        for e in Telem4.sub_elements():
        	e._quad_scheme = 'default'
        ####### Mixed element for rigid body motion #####################################
        VRelem = MixedElement([Relem, Relem, Relem, Relem, Relem])
        for e in VRelem.sub_elements():
        	e._quad_scheme = 'default'
        #################################################################################
        #W = FunctionSpace(self.mesh, MixedElement([Velem,Qelem,Relem]))
        #W = FunctionSpace(self.mesh, MixedElement([Velem,Relem,VRelem]))
        if inverseHeart:
            Mixelem=MixedElement([Velem,VRelem])
        else:
            Mixelem=MixedElement([Velem,Qelem,Relem,VRelem])
        Mixelem._quad_scheme = 'default'
        for e in Mixelem.sub_elements():
        	e._quad_scheme = 'default'
        W = FunctionSpace(self.mesh, Mixelem)
        self.invFspace=FunctionSpace(self.mesh, Telem3)
        self.displacementSpace=FunctionSpace(self.mesh, Velem)
        self.phaseSpace=FunctionSpace(self.mesh, Qelem)
        self.dphaseSpace=FunctionSpace(self.mesh, Qelem)
        self.trSpace=FunctionSpace(self.mesh, Qelem)
        self.tr_prevSpace=FunctionSpace(self.mesh, Qelem)
        Quad = FunctionSpace(self.mesh, Quadelem)
        
        bctop = DirichletBC(W.sub(0).sub(2), Expression(("0.0"), degree = 4), self.facetboundaries, topid)
        bcs = [bctop]
        self.invF= Function(self.invFspace)
        self.invF_prev= Function(self.invFspace)
        self.phase= Function(self.phaseSpace)
        self.dphase= Function(self.dphaseSpace)
        self.tr= Function(self.trSpace)
        self.tr_prev= Function(self.tr_prevSpace)
        self.w = Function(W)
        dw = TrialFunction(W)
        wtest = TestFunction(W)
        
        if inverseHeart:
            du,dc11 = TrialFunctions(W)
            (self.u,c11) = split(self.w)
            (v,v11) = TestFunctions(W)
            self.p=None
            pendo=None
        else:
            du,dp, dpendo,dc11 = TrialFunctions(W)
            (self.u,self.p, pendo,c11) = split(self.w)
            (v,q, qendo,v11) = TestFunctions(W)
        
        PK1activestress = Function(self.Q)
        PK1activestress.rename("active stress", "active stress")
        
        self.t_a = Expression(("t_a"), t_a=0.0, degree=1)
        self.dt = Expression(("dt"), dt=0.0, degree=1)
        T0_LV = runParameters['T0_LV']
        
        params= {"mesh": self.mesh,
                 "facetboundaries": self.facetboundaries,
                 "facet_normal": N,
        	 "mixedfunctionspace": W,
        	 "mixedfunction": self.w,
                 "displacement_variable": self.u, 
                 "pressure_variable": self.p,
        	 "volconst_variable": pendo,
        	 "constrained_vol":self.Cavityvol,
                 "endoid": endoid,
        	 "fiber": self.f0,
                 "sheet": self.s0,
                 "sheet-normal": self.n0,
                 "invF_variable":self.invF,
                 "cstrainfactor": self.cstrainFactor,
                 "StrainEnergyDensityFunction_Coef":runParameters["StrainEnergyDensityFunction_Coef"],
                 "StrainEnergyDensityFunction_Cff":runParameters["StrainEnergyDensityFunction_Cff"],
                 "StrainEnergyDensityFunction_Css":runParameters["StrainEnergyDensityFunction_Css"],
                 "StrainEnergyDensityFunction_Cnn":runParameters["StrainEnergyDensityFunction_Cnn"],
                 "StrainEnergyDensityFunction_Cns":runParameters["StrainEnergyDensityFunction_Cns"],
                 "StrainEnergyDensityFunction_Cfs":runParameters["StrainEnergyDensityFunction_Cfs"],
                 "StrainEnergyDensityFunction_Cfn":runParameters["StrainEnergyDensityFunction_Cfn"]}
        
        activeparams = {"mesh": self.mesh,
                        "facetboundaries": self.facetboundaries,
                        "facet_normal": N,
                        "displacement_variable": self.u, 
                        "pressure_variable": self.p,
                        "relaxphase":self.phase,
                        "drelaxphase":self.dphase,
                        "tr": self.tr,
                        "prev_tr": self.tr_prev,
                        "fiber": self.f0,
                        "sheet": self.s0,
                        "sheet-normal": self.n0,
                "BCL": runParameters['BCL'],
        		"t_a": self.t_a,
        		"dt": self.dt,
                "activefiberstrength": self.activefiberstrength,
                "activefiberdelay": self.activefiberdelay,
        		"Tact": T0_LV,
        		"T0": T0_LV}
        
         
        self.uflforms = Forms(params)
        self.activeforms = Active(activeparams,trackphase=self.trackphase)
        self.activeforms.set_default_parameters(runParameters)
        
        self.Fmat = self.uflforms.Fmat()
        self.Cmat = (self.Fmat.T*self.Fmat)
        self.Emat = self.uflforms.Emat()
        J = self.uflforms.J()
        
        n = J*inv(self.Fmat.T)*N
        dx = dolfin.dx(self.mesh,metadata = {"integration_order":2})
        
        Ematrix = project(self.Emat, TF)
        
        
        
        # Automatic differentiation  #####################################################################################################
        if inverseHeart:
            Fmat_inv = self.uflforms.Fmat(inverse=True)


            F4 = inner(Fmat_inv, grad(v))*dx
            L5 = inner(runParameters['pinned_movement']*as_vector([c11[0], c11[1], 0.0]), self.u)*dx
            L6 = inner(runParameters['pinned_movement']*as_vector([0.0, 0.0, c11[2]]), cross(X, self.u))*dx + \
            	 inner(runParameters['pinned_movement']*as_vector([c11[3], 0.0, 0.0]), cross(X, self.u))*dx + \
            	 inner(runParameters['pinned_movement']*as_vector([0.0, c11[4], 0.0]), cross(X, self.u))*dx
            F5 = derivative(L5, self.w, wtest)
            F6 = derivative(L6, self.w, wtest)
            Ftotal = F4 + F5 + F6
            Jac4 = derivative(F4, self.w, dw) 
            Jac5 = derivative(F5, self.w, dw)
            Jac6 = derivative(F6, self.w, dw)
            Jac = Jac4 + Jac5 + Jac6
        else:
            Wp = self.uflforms.PassiveMatSEF()
            Wvol = self.uflforms.V0constrainedE()
            Pactive = self.activeforms.PK1StressTensor()
        
            F1 = derivative(Wp, self.w, wtest)
            F2 = derivative(Wvol, self.w, wtest)
            F3 = Kspring*inner(dot(self.u,n)*n,v)/Constant(runParameters['length_scale']**2.)*self.ds(epiid)
            
            F4 = inner(self.Fmat*Pactive, grad(v))*dx
            L5 = inner(runParameters['pinned_movement']*as_vector([c11[0], c11[1], 0.0])*Constant(runParameters['pressure_scale']), self.u)*dx
            L6 = inner(runParameters['pinned_movement']*as_vector([0.0, 0.0, c11[2]])*Constant(runParameters['pressure_scale']/runParameters['length_scale']), cross(X, self.u))*dx + \
            	 inner(runParameters['pinned_movement']*as_vector([c11[3], 0.0, 0.0])*Constant(runParameters['pressure_scale']/runParameters['length_scale']), cross(X, self.u))*dx + \
            	 inner(runParameters['pinned_movement']*as_vector([0.0, c11[4], 0.0])*Constant(runParameters['pressure_scale']/runParameters['length_scale']), cross(X, self.u))*dx
            F5 = derivative(L5, self.w, wtest)
            F6 = derivative(L6, self.w, wtest)
            Ftotal = F1 + F2 + F3 + F4 + F5 +F6
            Jac1 = derivative(F1, self.w, dw) 
            Jac2 = derivative(F2, self.w, dw) 
            Jac3 = derivative(F3, self.w, dw) 
            Jac4 = derivative(F4, self.w, dw) 
            Jac5 = derivative(F5, self.w, dw)
            Jac6 = derivative(F6, self.w, dw)
            Jac = Jac1 + Jac2 + Jac3 + Jac4 + Jac5 + Jac6
        ##################################################################################################################################
        
        solverparams = {"Jacobian": Jac,
                        "F": Ftotal,
                        "w": self.w,
                        "boundary_conditions": bcs,
        		"Type": 0,
        		"mesh": self.mesh,
                "debug_call": self.debug_call
        		}
        
        
        
        self.solver= NSolver(solverparams)
        self.mesh_volume = assemble(Constant(1.0)*dx)
        return
    def updatePhase(self):
        self.phase.vector()[:]=self.phase.vector()[:]+self.dphase.vector()[:]
        for n in range(len(self.phase.vector()[:])):
            if self.dphase.vector()[n]<0:
                self.phase.vector()[n]=0
            elif (self.phase.vector()[n]+self.dphase.vector()[n])>np.pi:
                self.phase.vector()[n]=np.pi
            else:
                self.phase.vector()[n]=self.phase.vector()[n]+self.dphase.vector()[n]
        self.tr_prev.vector()[:]=self.tr.vector()[:].copy()
    def copy(self):
        return heart(*self.backup)
    
class heart_AV:
    def __new__(cls, *args, **kwargs):
        return super().__new__(cls)
    def __init__(self,casename,meshname,runParameters,inverseHeart=False,trackphase=False,AVJ_control=0):
        self.debug_call=[]
        self.inverseHeart=inverseHeart
        self.trackphase=trackphase
        self.backup=(casename,meshname,dict(runParameters))
        runParameters=dict(runParameters)
        parameters["form_compiler"]["quadrature_degree"]=4
        parameters["form_compiler"]["representation"] = "uflacs"
        parameters['ghost_mode'] = 'shared_facet'
        
        self.mesh = Mesh()
        print("reading mesh",casename+'/'+meshname+".hdf5")
        f = HDF5File(MPI.comm_world, casename+'/'+meshname+".hdf5", 'r') 
        f.read(self.mesh, meshname, False)
        File(casename+"/facetboundaries.pvd.pvd")<<self.mesh
        
        self.facetboundaries = MeshFunction("size_t", self.mesh, 2)
        f.read(self.facetboundaries, meshname+"/"+"facetboundaries")
        
        #f.read(self.facetboundaries, meshname+"/"+"facetboundaries")
        self.ds = dolfin.ds(subdomain_data = self.facetboundaries)
        File(casename+"/"+"facetboundaries.pvd") << self.facetboundaries
        
        self.fiberFS = VectorFunctionSpace(self.mesh, 'DG', 0)
        VQuadelem = VectorElement("Quadrature", self.mesh.ufl_cell(), degree=4, quad_scheme="default") #segmentation fault happened due to degree = 2. degree = 4 should be implied. I changed it. (Joy) 
        VQuadelem._quad_scheme = 'default'
        for e in VQuadelem.sub_elements():
        	e._quad_scheme = 'default'
        SQuadelem = FiniteElement("Quadrature", self.mesh.ufl_cell(), degree=4, quad_scheme="default")
        SQuadelem._quad_scheme = 'default'
        scalarFS = FunctionSpace(mesh, SQuadelem)
        
        self.fiberFS = FunctionSpace(self.mesh, VQuadelem)
        self.f0 = Function(self.fiberFS)
        self.s0 = Function(self.fiberFS)
        self.n0 = Function(self.fiberFS)
        
        f.read(self.f0, meshname+"/"+"eF")
        f.read(self.s0, meshname+"/"+"eS")
        f.read(self.n0, meshname+"/"+"eN")
        
        self.long0 = Function(self.fiberFS)
        self.cir0 = Function(self.fiberFS)
        self.rad0 = Function(self.fiberFS)
        
        f.read(self.long0, meshname+"/"+"eL")
        f.read(self.cir0, meshname+"/"+"eC")
        f.read(self.rad0, meshname+"/"+"eR")
        
        try:
            self.cstrainFactor = Function(self.fiberFS)
            f.read(self.cstrainFactor, meshname+"/"+"cstrainfactor")
            print("cstrainFactor read.")
            print("cstrainFactor min-max:", self.cstrainFactor.sub(0).vector()[:].min(), self.cstrainFactor.sub(0).vector()[:].max())
        except Exception as e:
            self.cstrainFactor =None
            print("cstrainFactor NOT read.")
            print(repr(e))
        try:
            self.activefiberstrength = Function(self.fiberFS)
            f.read(self.activefiberstrength, meshname+"/"+"activefiberstrength")
            print("activefiberstrength read.")
            print("activefiberstrength min-max:", self.activefiberstrength.sub(0).vector()[:].min(), self.activefiberstrength.sub(0).vector()[:].max())
        except:
            self.activefiberstrength =None
            print("activefiberstrength NOT read.")
        try:
            self.activefiberdelay = Function(self.fiberFS)
            f.read(self.activefiberdelay, meshname+"/"+"activefiberdelay")
            print("activefiberdelay read.")
            print("activefiberdelay min-max:", self.activefiberdelay.sub(0).vector()[:].min(), self.activefiberdelay.sub(0).vector()[:].max())
        except:
            self.activefiberdelay = None
            print("activefiberdelay NOT read.")
        
        f.close()
        
        avjid = runParameters['avjid']
        endoAtrid = runParameters['endoAtrid']
        epiAtrid = runParameters['epiAtrid']
        endoVenid = runParameters['endoVenid']
        epiVenid = runParameters['epiVenid']
        endoAtrAVJid = runParameters['endoAtrAVJid']
        endoVenAVJid = runParameters['endoVenAVJid']
        #####################################################################
        self.comm = self.mesh.mpi_comm()
        
        
        N = FacetNormal (self.mesh)
        X = SpatialCoordinate (self.mesh)
        Kspring = Constant(runParameters['Kspring_constant'])
        self.PressAtr = Expression(("patr"), patr=0.0, degree=2)
        self.CavityvolAtr = Expression(("vol"), vol=0.0, degree=2)
        self.PressVen = Expression(("pven"), pven=0.0, degree=2)
        self.CavityvolVen = Expression(("vol"), vol=0.0, degree=2)
        if AVJ_control:
            self.PressAVJ = Expression(("pavj"), pavj=0.0, degree=2)
            self.CavityvolAVJ = Expression(("vol"), vol=0.0, degree=2)
        V = VectorFunctionSpace(self.mesh, 'CG', 2)
        TF = TensorFunctionSpace(self.mesh, 'DG', 1)
        self.Q = FunctionSpace(self.mesh,'CG',1)
        R = FunctionSpace(self.mesh,'Real',0)
        
        Velem = VectorElement("CG", self.mesh.ufl_cell(), 2, quad_scheme="default")
        Velem._quad_scheme = 'default'
        for e in Velem.sub_elements():
        	e._quad_scheme = 'default'
        Qelem = FiniteElement("CG", self.mesh.ufl_cell(), 1, quad_scheme="default")
        Qelem._quad_scheme = 'default'
        for e in Qelem.sub_elements():
        	e._quad_scheme = 'default'
        Relem = FiniteElement("Real", self.mesh.ufl_cell(), 0, quad_scheme="default")
        Relem._quad_scheme = 'default'
        for e in Relem.sub_elements():
        	e._quad_scheme = 'default'
        Quadelem = FiniteElement("Quadrature", self.mesh.ufl_cell(), degree=4, quad_scheme="default")
        Quadelem._quad_scheme = 'default'
        for e in Quadelem.sub_elements():
        	e._quad_scheme = 'default'
        
        Telem2 = TensorElement("Quadrature", self.mesh.ufl_cell(), degree=4, shape=2*(3,), quad_scheme='default')
        Telem2._quad_scheme = 'default'
        for e in Telem2.sub_elements():
        	e._quad_scheme = 'default'
        Telem3 = TensorElement("Lagrange", self.mesh.ufl_cell(), degree=2, shape=(3,3), quad_scheme='default')
        Telem3._quad_scheme = 'default'
        for e in Telem3.sub_elements():
        	e._quad_scheme = 'default'
        Telem4 = TensorElement("Quadrature", self.mesh.ufl_cell(), degree=4, shape=4*(3,), quad_scheme='default')
        Telem4._quad_scheme = 'default'
        for e in Telem4.sub_elements():
        	e._quad_scheme = 'default'
        ####### Mixed element for rigid body motion #####################################
        if AVJ_control==3:
            VRelem = MixedElement([Relem, Relem, Relem, Relem, Relem, Relem, Relem, Relem, Relem, Relem, Relem])
        else:
            VRelem = MixedElement([Relem, Relem, Relem, Relem, Relem])
        for e in VRelem.sub_elements():
        	e._quad_scheme = 'default'
        #################################################################################
        #W = FunctionSpace(self.mesh, MixedElement([Velem,Qelem,Relem]))
        #W = FunctionSpace(self.mesh, MixedElement([Velem,Relem,VRelem]))
        if inverseHeart:
            Mixelem=MixedElement([Velem,VRelem])
        elif AVJ_control:
            Mixelem=MixedElement([Velem,Qelem,Relem,Relem,Relem,VRelem])
        else:
            Mixelem=MixedElement([Velem,Qelem,Relem,Relem,VRelem])
        Mixelem._quad_scheme = 'default'
        for e in Mixelem.sub_elements():
        	e._quad_scheme = 'default'
        W = FunctionSpace(self.mesh, Mixelem)
        self.invFspace=FunctionSpace(self.mesh, Telem3)
        self.displacementSpace=FunctionSpace(self.mesh, Velem)
        self.phaseSpace=FunctionSpace(self.mesh, Qelem)
        self.dphaseSpace=FunctionSpace(self.mesh, Qelem)
        self.trSpace=FunctionSpace(self.mesh, Qelem)
        self.tr_prevSpace=FunctionSpace(self.mesh, Qelem)
        Quad = FunctionSpace(self.mesh, Quadelem)
        
        bcavj = DirichletBC(W.sub(0).sub(2), Expression(("0.0"), degree = 4), self.facetboundaries, avjid)
        bcs = [bcavj]
        self.invF= Function(self.invFspace)
        self.invF_prev= Function(self.invFspace)
        self.phase= Function(self.phaseSpace)
        self.dphase= Function(self.dphaseSpace)
        self.tr= Function(self.trSpace)
        self.tr_prev= Function(self.tr_prevSpace)
        self.w = Function(W)
        dw = TrialFunction(W)
        wtest = TestFunction(W)
        
        if inverseHeart:
            du,dc11 = TrialFunctions(W)
            (self.u,c11) = split(self.w)
            (v,v11) = TestFunctions(W)
            self.p=None
            pendoAtr=None
            pendoVen=None
            if AVJ_control:
                pendoAVJ=None
        elif AVJ_control:
            du,dp, dpendoAtr,dpendoVen,dpendoAVJ,dc11 = TrialFunctions(W)
            (self.u,self.p, pendoAtr,pendoVen,pendoAVJ,c11) = split(self.w)
            (v,q, qendoAtr,qendoVen,qendoAVJ,v11) = TestFunctions(W)
        else:
            du,dp, dpendoAtr,dpendoVen,dc11 = TrialFunctions(W)
            (self.u,self.p, pendoAtr,pendoVen,c11) = split(self.w)
            (v,q, qendoAtr,qendoVen,v11) = TestFunctions(W)
        
        PK1activestress = Function(self.Q)
        PK1activestress.rename("active stress", "active stress")
        
        self.t_a = Expression(("t_a"), t_a=0.0, degree=1)
        self.dt = Expression(("dt"), dt=0.0, degree=1)
        T0_LV = runParameters['T0_LV']
        
        params= {"mesh": self.mesh,
                 "facetboundaries": self.facetboundaries,
                 "facet_normal": N,
        	 "mixedfunctionspace": W,
        	 "mixedfunction": self.w,
                 "displacement_variable": self.u, 
                 "pressure_variable": self.p,
        	 "volAtrconst_variable": pendoAtr,
        	 "volVenconst_variable": pendoVen,
        	 "constrained_volAtr":self.CavityvolAtr,
        	 "constrained_volVen":self.CavityvolVen,
                 "endoAtrid": endoAtrid,
                  "endoVenid": endoVenid,
                  "endoAtrAVJid": endoAtrAVJid,
                  "endoVenAVJid": endoVenAVJid,
        	 "fiber": self.f0,
                 "sheet": self.s0,
                 "sheet-normal": self.n0,
                 "invF_variable":self.invF,
                 "cstrainfactor": self.cstrainFactor,
                 "AVJ_control": 0,
                 "StrainEnergyDensityFunction_Coef":runParameters["StrainEnergyDensityFunction_Coef"],
                 "StrainEnergyDensityFunction_Cff":runParameters["StrainEnergyDensityFunction_Cff"],
                 "StrainEnergyDensityFunction_Css":runParameters["StrainEnergyDensityFunction_Css"],
                 "StrainEnergyDensityFunction_Cnn":runParameters["StrainEnergyDensityFunction_Cnn"],
                 "StrainEnergyDensityFunction_Cns":runParameters["StrainEnergyDensityFunction_Cns"],
                 "StrainEnergyDensityFunction_Cfs":runParameters["StrainEnergyDensityFunction_Cfs"],
                 "StrainEnergyDensityFunction_Cfn":runParameters["StrainEnergyDensityFunction_Cfn"],}
        if AVJ_control:
            params["volAVJconst_variable"]=pendoAVJ
            params["constrained_volAVJ"]=self.CavityvolAVJ
            params["AVJ_control"]=AVJ_control
            if AVJ_control==3:
                params["c11"]=c11
        activeparams = {"mesh": self.mesh,
                        "facetboundaries": self.facetboundaries,
                        "facet_normal": N,
                        "displacement_variable": self.u, 
                        "pressure_variable": self.p,
                        "relaxphase":self.phase,
                        "drelaxphase":self.dphase,
                        "tr": self.tr,
                        "prev_tr": self.tr_prev,
                        "fiber": self.f0,
                        "sheet": self.s0,
                        "sheet-normal": self.n0,
                "BCL": runParameters['BCL'],
        		"t_a": self.t_a,
        		"dt": self.dt,
                "activefiberstrength": self.activefiberstrength,
                "activefiberdelay": self.activefiberdelay,
        		"Tact": T0_LV,
        		"T0": T0_LV}
        
         
        self.uflforms = Forms_AV(params)            
        self.activeforms = Active(activeparams,trackphase=self.trackphase)
        self.activeforms.set_default_parameters(runParameters)
        
        self.Fmat = self.uflforms.Fmat()
        self.Cmat = (self.Fmat.T*self.Fmat)
        self.Emat = self.uflforms.Emat()
        J = self.uflforms.J()
        
        n = J*inv(self.Fmat.T)*N
        dx = dolfin.dx(self.mesh,metadata = {"integration_order":2})
        
        Ematrix = project(self.Emat, TF)
        
        
        
        # Automatic differentiation  #####################################################################################################
        dsavj=ds(avjid, domain = self.mesh, subdomain_data = self.facetboundaries)
        if inverseHeart:
            Fmat_inv = self.uflforms.Fmat(inverse=True)


            F4 = inner(Fmat_inv, grad(v))*dx
            L5 = inner(runParameters['pinned_movement']*as_vector([c11[0], c11[1], 0.0]), self.u)*dx
            L6 = inner(runParameters['pinned_movement']*as_vector([0.0, 0.0, c11[2]]), cross(X, self.u))*dx + \
            	 inner(runParameters['pinned_movement']*as_vector([c11[3], 0.0, 0.0]), cross(X, self.u))*dx + \
            	 inner(runParameters['pinned_movement']*as_vector([0.0, c11[4], 0.0]), cross(X, self.u))*dx
            F5 = derivative(L5, self.w, wtest)
            F6 = derivative(L6, self.w, wtest)
            Ftotal = F4 + F5 + F6
            if AVJ_control==3:
                L7=inner( X+self.u-Constant(0.5)*as_vector([c11[5], c11[6],  c11[7]]),as_vector([c11[5], c11[6],  c11[7]]))*self.ds(endoAtrid)
                L8=inner( X+self.u-Constant(0.5)*as_vector([c11[8], c11[9],  c11[10]]), as_vector([c11[8], c11[9],  c11[10]]))*self.ds(endoVenid)
                F7=derivative(L7, self.w, wtest)
                F8=derivative(L8, self.w, wtest)
                Ftotal=Ftotal+F7+F8
            Jac4 = derivative(F4, self.w, dw) 
            Jac5 = derivative(F5, self.w, dw)
            Jac6 = derivative(F6, self.w, dw)
            Jac = Jac4 + Jac5 + Jac6
            if AVJ_control==3:
                Jac=Jac+derivative(F7, self.w, dw)+derivative(F8, self.w, dw)
            F1=None
            F2=None
            F3_Atr=None
            F3_Ven=None
            Jac1=None
            Jac2=None
            Jac3_Atr=None
            Jac3_Ven=None
        else:
            Wp = self.uflforms.PassiveMatSEF()
            #Wvol_Atr = self.uflforms.V0constrainedE(addstr="Atr")
            #Wvol_Ven = self.uflforms.V0constrainedE(addstr="Ven")
            Wvol = self.uflforms.V0constrainedE()
            Pactive = self.activeforms.PK1StressTensor()
        
            F1 = derivative(Wp, self.w, wtest)
            #F2_Atr = derivative(Wvol_Atr, self.w, wtest)
            #F2_Ven = derivative(Wvol_Ven, self.w, wtest)
            F2 = derivative(Wvol, self.w, wtest)
            F3_Atr = Kspring*inner(dot(self.u,n)*n,v)/Constant(runParameters['length_scale']**2.)*self.ds(epiAtrid)
            F3_Ven = Kspring*inner(dot(self.u,n)*n,v)/Constant(runParameters['length_scale']**2.)*self.ds(epiVenid)
            F4 = inner(self.Fmat*Pactive, grad(v))*dx
            L5 = inner(runParameters['pinned_movement']*as_vector([c11[0], c11[1], 0.0])*Constant(runParameters['pressure_scale']), self.u)*self.ds#dsavj
            L6 = inner(runParameters['pinned_movement']*as_vector([0.0, 0.0, c11[2]])*Constant(runParameters['pressure_scale']/runParameters['length_scale']), cross(X, self.u))*self.ds + \
            	 inner(runParameters['pinned_movement']*as_vector([c11[3], 0.0, 0.0])*Constant(runParameters['pressure_scale']/runParameters['length_scale']), cross(X, self.u))*self.ds + \
            	 inner(runParameters['pinned_movement']*as_vector([0.0, c11[4], 0.0])*Constant(runParameters['pressure_scale']/runParameters['length_scale']), cross(X, self.u))*self.ds
            F5 = derivative(L5, self.w, wtest)
            F6 = derivative(L6, self.w, wtest)
            #Ftotal = F1 + F4 + F5 +F6
            #Ftotal = F1 + F2_Atr + F2_Ven + F3_Atr + F3_Ven + F4 + F5 +F6
            Ftotal = F1 + F2 + F3_Atr + F3_Ven + F4 + F5 +F6###original
            #Ftotal = F1 + F2 + F3_Atr + F3_Ven + F5 +F6
            if AVJ_control==3:
                L7=inner( X+self.u-as_vector([c11[5], c11[6],  c11[7]]), X+self.u-as_vector([c11[5], c11[6],  c11[7]]))*self.ds(endoAtrid)
                L8=inner( X+self.u-as_vector([c11[8], c11[9],  c11[10]]), X+self.u-as_vector([c11[8], c11[9],  c11[10]]))*self.ds(endoVenid)
                F7=derivative(L7, self.w, wtest)
                F8=derivative(L8, self.w, wtest)
                Ftotal=Ftotal+F7+F8
            Jac1 = derivative(F1, self.w, dw) 
            #Jac2_Atr = derivative(F2_Atr, self.w, dw) 
            #Jac2_Ven = derivative(F2_Ven, self.w, dw)
            Jac2 = derivative(F2, self.w, dw)
            Jac3_Atr = derivative(F3_Atr, self.w, dw) 
            Jac3_Ven = derivative(F3_Ven, self.w, dw) 
            Jac4 = derivative(F4, self.w, dw) 
            Jac5 = derivative(F5, self.w, dw)
            Jac6 = derivative(F6, self.w, dw)
            #Jac = Jac1 + Jac4 + Jac5 + Jac6
            Jac = Jac1 + Jac2 + Jac3_Atr + Jac3_Ven + Jac4 + Jac5 + Jac6#original
            #Jac = Jac1 + Jac2 + Jac3_Atr + Jac3_Ven + Jac5 + Jac6
            if AVJ_control==3:
                Jac=Jac+derivative(F7, self.w, dw)+derivative(F8, self.w, dw)
        ##################################################################################################################################
        
        solverparams = {"Jacobian": Jac,
                        "F": Ftotal,
                        "F1": F1,
                        #"F2_Atr": F2_Atr,
                        #"F2_Ven": F2_Ven,
                        "F2": F2,
                        "F3_Atr": F3_Atr,
                        "F3_Ven": F3_Ven,
                        "F4": F4,
                        "F5": F5,
                        "F6": F6,
                        "Jac1": Jac1,
                        #"Jac2_Atr": Jac2_Atr,
                        #"Jac2_Ven": Jac2_Ven,
                        "Jac2": Jac2,
                        "Jac3_Atr": Jac3_Atr,
                        "Jac3_Ven": Jac3_Ven,
                        "Jac4": Jac4,
                        "Jac5": Jac5,
                        "Jac6": Jac6,
                        "w": self.w,
                        "boundary_conditions": bcs,
                        'lr':1.,
                        "max_iter":1000,
        		"Type": runParameters["solver type"],
        		"mesh": self.mesh,
                "debug_call": self.debug_call
        		}
        
        
        
        self.solver= NSolver(solverparams)
        self.mesh_volume = assemble(Constant(1.0)*dx)
        return
    def updatePhase(self):
        self.phase.vector()[:]=self.phase.vector()[:]+self.dphase.vector()[:]
        for n in range(len(self.phase.vector()[:])):
            if self.dphase.vector()[n]<0:
                self.phase.vector()[n]=0
            elif (self.phase.vector()[n]+self.dphase.vector()[n])>np.pi:
                self.phase.vector()[n]=np.pi
            else:
                self.phase.vector()[n]=self.phase.vector()[n]+self.dphase.vector()[n]
        self.tr_prev.vector()[:]=self.tr.vector()[:].copy()
    def copy(self):
        return heart(*self.backup)