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
                                                            - added strainEnergy in ulforms 
'''
_version='3.4.3'

import sys

import os as os
#sys.path.append(os.path.dirname(os.path.realpath(__file__)))
#import shutil
try:
    from dolfin import *
except:
    from fenics import *
import numpy as np
from heartFEM.lcleeHeart.forms_shavik import Forms
from heartFEM.lcleeHeart.simpleActiveMaterial import SimpleActiveMaterial as Active
from heartFEM.lcleeHeart.nsolver import NSolver as NSolver 
from heartFEM.lcleeHeart.addfiber_matid import *

import vtk
from heartFEM.lcleeHeart import vtk_py as vtk_py
from mpi4py import MPI as pyMPI
from heartFEM.lcleeHeart.rotateUGrid_w_axis import rotateUGrid_w_axis
from pyquaternion import Quaternion

def generateMesh(casePath,stlname,Laxis,endo_angle='',epi_angle='',fiberSheetletAngle=0.,fiberSheetletWidth=0.,radialFiberAngle=0.,fiberLength=0.,clipratio=0.75,meshsize=0.6,meshname=None,saveSubFolder=''):
    if meshname is None:
        meshname=stlname
    Laxis = Laxis / np.linalg.norm(Laxis)
    meshfilename = casePath+ '/'+stlname+ '.stl'
    pdata = vtk_py.readSTL(meshfilename)
        
    print(pdata)

    C = vtk_py.getcentroid(pdata)

    angle = np.arccos(Laxis[2])
    if Laxis[2]>0.99:
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
    
    if not(isinstance(clipratio,float)):
        clipratio=(transformation_matrix_R.dot(clipratio.reshape((-1,1))).reshape(-1)[2]-rotatedpdata.GetBounds()[4])/(rotatedpdata.GetBounds()[5]-rotatedpdata.GetBounds()[4])
    print(rotatedpdata.GetBounds())
    vtk_py.writeSTL(rotatedpdata, casePath+saveSubFolder+"/TEST.stl")
    isinsideout = 1
    C = vtk_py.getcentroid(rotatedpdata)
    C = [C[0], C[1], rotatedpdata.GetBounds()[4]+(rotatedpdata.GetBounds()[5]-rotatedpdata.GetBounds()[4])*clipratio] #you have xmin, xmax, ymin,ymax,zmin,zmax (is just to draw a box to clip the geometry, then you can get the edge of boundary)
	
    print(C)
    clippedheart = vtk_py.clipheart(rotatedpdata, C, [0,0,1],isinsideout, True) 

    vtk_py.writeSTL(clippedheart, casePath+saveSubFolder+"/TEST2.stl")
	
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
   
    vtk_py.writeSTL(cleanepi, epifile)
    vtk_py.writeSTL(cleanendo, endofile)
   
    filename = casePath+saveSubFolder+'/New_mesh' 
   
    vtk_py.createLVmesh(filename, meshsize, casePath+saveSubFolder+"/epi_lv.stl", casePath+saveSubFolder+"/endo_lv.stl")
   
   
    ugridfilename = filename + '.vtk'
    pdatafilename = filename + '.vtp'
    ugrid = vtk_py.readUGrid(ugridfilename)
   
    mesh = vtk_py.convertUGridToXMLMesh(ugrid)
   	
    print (mesh)
   
    #print (mesh)
    comm2 = pyMPI.COMM_WORLD
   
    # Rotate mesh (LCL include)
    rugrid = vtk_py.rotateUGrid(ugrid, rx=360)#original is rx=180, so it is not working, so I change it to rx=360, it works. 
    rugrid = vtk_py.rotateUGrid(rugrid, sx=0.1, sy=0.1, sz=0.1) #add the scale, the code did not working, so I remove it it works 
    transformation_matrix_S=np.array([[0.1,0,0,0],[0,0.1,0,0],[0,0,0.1,0],[0,0,0,1.]])
    
    vtk_py.writeUGrid(rugrid,casePath+saveSubFolder+"/rtemp.vtk")
   
   	
    fenics_mesh_ref, fenics_facet_ref, fenics_edge_ref = vtk_py.extractFeNiCsBiVFacet(rugrid,savePath=casePath+saveSubFolder+'/', geometry = "LV")
   
   	
    matid = MeshFunction('size_t',fenics_mesh_ref, 3, mesh.domains())
   
   	
    ztop =  max(fenics_mesh_ref.coordinates()[:,2])
    center_to_shift=cleanepi.GetCenter()
    print("DISPLAY center_to_shift",center_to_shift)
    
    transformation_matrix=transformation_matrix_S.dot(transformation_matrix_R)
    transformation_matrix[0,3]=-center_to_shift[0]/10.
    transformation_matrix[1,3]=-center_to_shift[1]/10.
    transformation_matrix[2,3]=-ztop
    
    np.savetxt(casePath+saveSubFolder+'/'+meshname+'_Tmatrix_fromSTLtoFEniCS.txt',transformation_matrix)
    invtransformation_matrix=transformation_matrix_R.T.dot(np.array([[10.,0,0,0],[0,10.,0,0],[0,0,10.,0],[0,0,0,1.]]))
    inv_translate=invtransformation_matrix.dot(transformation_matrix[:,3:])
    invtransformation_matrix[:3,3]=-inv_translate[:3,0]
    np.savetxt(casePath+saveSubFolder+'/'+meshname+'_Tmatrix_fromFEniCStoSTL.txt',invtransformation_matrix)
    
    ztrans = Expression((str(-center_to_shift[0]/10.), str(-center_to_shift[1]/10.), str(-ztop)), degree = 1)
    #ztrans = Expression(("0.0", "0.0", str(-ztop)), degree = 1)
   
   
    ALE.move(fenics_mesh_ref,ztrans) #if (dolfin.dolfin_version() != '1.6.0'): fenics_mesh_ref.move(ztrans)
   
   		
   
    mesh = fenics_mesh_ref
   
    gdim = mesh.geometry().dim()
   
    outdir = casePath+saveSubFolder+"/"
   
    quad_deg = 4
    VQuadelem = VectorElement("Quadrature", mesh.ufl_cell(), degree=quad_deg, quad_scheme="default")
    VQuadelem._quad_scheme = 'default'
    fiberFS = FunctionSpace(mesh, VQuadelem)
    isepiflip = True #False
    isendoflip = True #True
    casedir=casePath+saveSubFolder+"/"+meshname+"_";  
   
    ef, es, en, eC, eL, eR = vtk_py.addLVfiber(mesh, fiberFS, "lv", endo_angle, epi_angle, casedir,fiberSheetletAngle=fiberSheetletAngle,fiberSheetletWidth=fiberSheetletWidth,radialFiberAngle=radialFiberAngle,fiberLength=fiberLength)
   
    matid_filename = casePath+saveSubFolder+'/'+meshname + "_matid.pvd"
    File(matid_filename) << matid
   	
   
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
   
    f.close()
   
   
    File(casePath+saveSubFolder+"/facetboundaries"+".pvd") << fenics_facet_ref
    File(casePath+saveSubFolder+"/edgeboundaries"+".pvd") << fenics_edge_ref
    File(casePath+saveSubFolder+"/mesh" + ".pvd") << mesh
    File(casePath+saveSubFolder+"/matid" +".pvd") << matid
   
   
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
def changeFiberAngles(casePath,meshname,endo_angle,epi_angle,fiberSheetletAngle=0.,fiberSheetletWidth=0.,radialFiberAngle=0.,fiberLength=0.,saveSubFolder='',loadSubFolder=''):
    os.makedirs(casePath+saveSubFolder,exist_ok=True)
    mesh = Mesh()
    f = HDF5File(MPI.comm_world, casePath+loadSubFolder+'/'+meshname+".hdf5", 'r') 
    f.read(mesh, meshname, False)
    
    facetboundaries = MeshFunction("size_t", mesh, 2)
    f.read(facetboundaries, meshname+"/"+"facetboundaries")
    
    edgeboundaries = MeshFunction('size_t', mesh,1)
    f.read(edgeboundaries, meshname+"/"+"edgeboundaries") 
    
    f.close()
    
    gdim = mesh.geometry().dim()
   
    quad_deg = 4
    VQuadelem = VectorElement("Quadrature", mesh.ufl_cell(), degree=quad_deg, quad_scheme="default")
    VQuadelem._quad_scheme = 'default'
    fiberFS = FunctionSpace(mesh, VQuadelem)
    casedir=casePath+saveSubFolder+"/"+meshname+"_"
    ef, es, en, eC, eL, eR = vtk_py.addLVfiber(mesh, fiberFS, "lv", endo_angle, epi_angle, casedir,fiberSheetletAngle=fiberSheetletAngle,fiberSheetletWidth=fiberSheetletWidth,radialFiberAngle=radialFiberAngle,fiberLength=fiberLength)
    matid = MeshFunction('size_t',mesh, 3, mesh.domains())
    matid_filename = casePath+saveSubFolder+'/'+meshname + "_matid.pvd"
    if os.path.isfile(matid_filename):
        os.remove(matid_filename)
    File(matid_filename) << matid
   	
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
   
    f.close()
   
   
    File(casePath+saveSubFolder+"/facetboundaries"+".pvd") << facetboundaries
    File(casePath+saveSubFolder+"/edgeboundaries"+".pvd") << edgeboundaries
    File(casePath+saveSubFolder+"/mesh" + ".pvd") << mesh
    File(casePath+saveSubFolder+"/matid" +".pvd") << matid
   
    return 0
class heart:
    def __new__(cls, *args, **kwargs):
        return super().__new__(cls)
    def __init__(self,casename,meshname,runParameters,PVinput='volume',inverseHeart=False):
        self.inverseHeart=inverseHeart
        self.backup=(casename,meshname,dict(runParameters),PVinput)
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
        
        fiberFS = VectorFunctionSpace(self.mesh, 'DG', 0)
        VQuadelem = VectorElement("Quadrature", self.mesh.ufl_cell(), degree=4, quad_scheme="default") #segmentation fault happened due to degree = 2. degree = 4 should be implied. I changed it. (Joy) 
        VQuadelem._quad_scheme = 'default'
        for e in VQuadelem.sub_elements():
        	e._quad_scheme = 'default'
        
        fiberFS = FunctionSpace(self.mesh, VQuadelem)
        self.f0 = Function(fiberFS)
        self.s0 = Function(fiberFS)
        self.n0 = Function(fiberFS)
        
        f.read(self.f0, meshname+"/"+"eF")
        f.read(self.s0, meshname+"/"+"eS")
        f.read(self.n0, meshname+"/"+"eN")
        
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
        Quad = FunctionSpace(self.mesh, Quadelem)
        
        bctop = DirichletBC(W.sub(0).sub(2), Expression(("0.0"), degree = 4), self.facetboundaries, topid)
        bcs = [bctop]
        self.invF= Function(self.invFspace)
        self.invF_prev= Function(self.invFspace)
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
        Tact = Constant(runParameters['Tact_constant'])
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
                        "fiber": self.f0,
                        "sheet": self.s0,
                        "sheet-normal": self.n0,
        		"t_a": self.t_a,
        		"dt": self.dt,
        		"Tact": T0_LV,
        		"T0": T0_LV}
        
         
        self.uflforms = Forms(params)
        self.activeforms = Active(activeparams)
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
            L5 = inner(as_vector([c11[0], c11[1], 0.0]), self.u)*dx
            L6 = inner(as_vector([0.0, 0.0, c11[2]]), cross(X, self.u))*dx + \
            	 inner(as_vector([c11[3], 0.0, 0.0]), cross(X, self.u))*dx + \
            	 inner(as_vector([0.0, c11[4], 0.0]), cross(X, self.u))*dx
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
            F3 = Kspring*inner(dot(self.u,n)*n,v)*self.ds(epiid)
            
            F4 = inner(self.Fmat*Pactive, grad(v))*dx
            L5 = inner(as_vector([c11[0], c11[1], 0.0]), self.u)*dx
            L6 = inner(as_vector([0.0, 0.0, c11[2]]), cross(X, self.u))*dx + \
            	 inner(as_vector([c11[3], 0.0, 0.0]), cross(X, self.u))*dx + \
            	 inner(as_vector([0.0, c11[4], 0.0]), cross(X, self.u))*dx
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
        		"mesh": self.mesh
        		}
        
        
        
        self.solver= NSolver(solverparams)
        self.mesh_volume = assemble(Constant(1.0)*dx)
        return
    def copy(self):
        return heart(*self.backup)
    
