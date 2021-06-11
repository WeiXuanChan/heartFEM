import vtk
from heartFEM.lcleeHeart import vtk_py as vtk_py
import glob, os
import numpy as np
import operator
from vtk.util.numpy_support import vtk_to_numpy
from vtk.util.numpy_support import numpy_to_vtk
from matplotlib import pylab as plt

def removetop(pdata):
	
	np_normal_vec = vtk_to_numpy(pdata.GetCellData().GetArray("Normals"))
	np_pts_vec = vtk_to_numpy(vtk_py.getCellCenters(pdata).GetPoints().GetData())

	cnt = 0
	markcell = []
	for vec, pt in zip(np_normal_vec, np_pts_vec):
		if(pt[2] > pdata.GetBounds()[5] - 1 and abs(np.dot(vec, np.array([0,0,1])) - 1) < 1e-2):
			pdata.DeleteCell(cnt)
			markcell.append(cnt)
			
		cnt += 1

  	
	pdata.RemoveDeletedCells()
	return markcell,pdata




def computeDeformationGradient_n_strain(pdata):

	n_cells = pdata.GetNumberOfCells()

	pdata.GetPointData().SetActiveVectors("disp")
	cell_derivatives = vtk.vtkCellDerivatives()
	cell_derivatives.SetVectorModeToPassVectors()
	cell_derivatives.SetTensorModeToComputeGradient()
	cell_derivatives.SetInputData(pdata)
	cell_derivatives.Update()

	farray_gu = cell_derivatives.GetOutput().GetCellData().GetArray("VectorGradient")
	eC = pdata.GetCellData().GetArray("eC")
	eR = pdata.GetCellData().GetArray("eR")
	eL = pdata.GetCellData().GetArray("eL")
	
	farray_defo_grad = vtk.vtkFloatArray()
	farray_defo_grad.SetNumberOfComponents(9)
	farray_defo_grad.SetNumberOfTuples(n_cells)
	farray_defo_grad.SetName("DeformationGradient")
	pdata.GetCellData().AddArray(farray_defo_grad)

	farray_strain = vtk.vtkFloatArray()
	farray_strain.SetNumberOfComponents(6)
	farray_strain.SetNumberOfTuples(n_cells)
	farray_strain.SetName("Cartesian_strain")
	pdata.GetCellData().AddArray(farray_strain)

	farray_LVstrain = vtk.vtkFloatArray()
	farray_LVstrain.SetNumberOfComponents(6)
	farray_LVstrain.SetNumberOfTuples(n_cells)
	farray_LVstrain.SetName("LV_strain")
	pdata.GetCellData().AddArray(farray_LVstrain)

	I = np.eye(3)
	for k_cell in range(n_cells):
	    GU = np.reshape(farray_gu.GetTuple(k_cell), (3,3), order='C')
	    F = I + GU
	    farray_defo_grad.SetTuple(k_cell, np.reshape(F, 9, order='C'))

	    C = np.dot(np.transpose(F), F)
	    E = (C - I)/2
	    E_vec = np.array([E[0,0], E[1,1], E[2,2], E[0,1], E[0,2], E[1,2]])

	    eC_ = np.reshape(eC.GetTuple(k_cell), (3,1), order='C')
	    eL_ = np.reshape(eL.GetTuple(k_cell), (3,1), order='C')
	    eR_ = np.reshape(eR.GetTuple(k_cell), (3,1), order='C')

	    Ecc = np.dot(np.transpose(eC_), np.dot(E, eC_))[0,0]
	    Ell = np.dot(np.transpose(eL_), np.dot(E, eL_))[0,0]
	    Err = np.dot(np.transpose(eR_), np.dot(E, eR_))[0,0]
	    Ecr = np.dot(np.transpose(eC_), np.dot(E, eR_))[0,0]
	    Ecl = np.dot(np.transpose(eC_), np.dot(E, eL_))[0,0]
	    Erl = np.dot(np.transpose(eR_), np.dot(E, eL_))[0,0]
	    E_LV_vec = np.array([Ecc, Ell, Err, Ecr, Ecl, Erl])

	    farray_strain.SetTuple(k_cell, E_vec)
	    farray_LVstrain.SetTuple(k_cell, E_LV_vec)
	

def GetCentremass(pdata):

	centerOfMassFilter = vtk.vtkCenterOfMass()
  	centerOfMassFilter.SetInputData(pdata);
  	centerOfMassFilter.SetUseScalarsAsWeights(False);
  	centerOfMassFilter.Update();

  	return centerOfMassFilter.GetCenter();

def flip_vector(center, vec_array, pts_array):

	new_vec_array = np.copy(vec_array)
	cnt = 0
	for vec, pt in zip(vec_array, pts_array):
		if np.dot(pt - center, vec) < 0:
			new_vec_array[cnt,:] = -1.0*vec
		cnt += 1

	return new_vec_array

def get_ortho_basis(normal, apex_base_dir):

	eR = np.copy(normal)
	eC = np.copy(normal)
	eL = np.copy(normal)
	cnt = 0
	for vec in normal:
		eC[cnt,:] = np.cross(apex_base_dir, vec)
		if(np.sqrt(eC[cnt,0]**2 + eC[cnt,1]**2 + eC[cnt,2]**2) > 1e-10):
			eC[cnt,:] = eC[cnt,:] / np.sqrt(eC[cnt,0]**2 + eC[cnt,1]**2 + eC[cnt,2]**2)
		else:
			eC[cnt,:] = np.array([1,0,0])

		eL[cnt,:] = np.cross(eR[cnt,:], eC[cnt,:])
		if(np.sqrt(eL[cnt,0]**2 + eL[cnt,1]**2 + eL[cnt,2]**2) > 1e-10):
			eL[cnt,:] = eL[cnt,:] / np.sqrt(eL[cnt,0]**2 + eL[cnt,1]**2 + eL[cnt,2]**2)
		else:
			eL[cnt,:] = np.array([0,1,0])
		cnt += 1

	return eR, eC, eL



filedir = "../deformation/"

outdir = "../newdeformation/"


ntpt = 401 # Total time point
reftpt = 4421 # Reference time point

coord_array = []
disp_array = []

# Get all the stl files in the directory and coordinates
for tpt in np.arange(4421, ntpt+4422):
	print tpt
	filename = filedir + "geometry_"+str(tpt)+".stl"
	pdata = vtk_py.readSTL(filename)
	coord_array.append(vtk_to_numpy(pdata.GetPoints().GetData()))
	if(tpt == reftpt):
		pdata_ref = vtk_py.readSTL(filename)


# Compute displacement using reference time
#print coord_array[615]
[disp_array.append(coord_array[(tpt-4421)] - coord_array[(reftpt-4421)]) for tpt in np.arange(4421, ntpt+4422)]


# Get eC, eR, eL direction #################################################################
# Get Normal vector to surface
pdata_ref = vtk_py.getPDataNormals(pdata_ref)
# Get Centroid of the STL file
centroid = GetCentremass(pdata)
# Convert cell normal vector array to numpy array
np_normal_vec = vtk_to_numpy(pdata_ref.GetCellData().GetArray("Normals"))
# Convert cell centroid array to numpy array
np_pts_vec = vtk_to_numpy(vtk_py.getCellCenters(pdata_ref).GetPoints().GetData())
# Flip normal vector
new_np_normal_vec = flip_vector(centroid, np_normal_vec, np_pts_vec)
# Get eC, eR and eL 
np_eR, np_eC, np_eL = get_ortho_basis(new_np_normal_vec, np.array([0,0,1]))
# Convert numpy array to vtk
eR = numpy_to_vtk(num_array=np_eR, deep=True, array_type=vtk.VTK_FLOAT)
eC = numpy_to_vtk(num_array=np_eC, deep=True, array_type=vtk.VTK_FLOAT)
eL = numpy_to_vtk(num_array=np_eL, deep=True, array_type=vtk.VTK_FLOAT)
eR.SetName("eR")
eL.SetName("eL")
eC.SetName("eC")
pdata_ref.GetCellData().AddArray(eR)
pdata_ref.GetCellData().AddArray(eL)
pdata_ref.GetCellData().AddArray(eC)
# Remove top facet
markcell, pdata_ref = removetop(pdata_ref)
#vtk_py.writeXMLPData(pdata_ref, "temp.vtp")

mean_Ecc = []
mean_Ell = []
mean_Err = []

# Compute displacement using reference time
for tpt in np.arange(4421, ntpt+4421):
	newfilename = outdir + "geometry_"+str(tpt)+".vtp"
	sourcefilename = filedir + "geometry_"+str(tpt)+".stl"
	pdata = vtk_py.readSTL(sourcefilename)

	# Set up displacement field
	disp = numpy_to_vtk(num_array=disp_array[tpt-4421], deep=True, array_type=vtk.VTK_FLOAT)
	disp.SetName("disp")

	# Insert displacement field into array
	pdata.GetPointData().AddArray(disp)

	# Insert LV coordinates
	pdata.GetCellData().AddArray(eR)
	pdata.GetCellData().AddArray(eL)
	pdata.GetCellData().AddArray(eC)

	# Compute strain
	computeDeformationGradient_n_strain(pdata)

	# Remove top facet
	[pdata.DeleteCell(cellid) for cellid in markcell]
  	pdata.RemoveDeletedCells()

	# Compute average strain
	Ecc = np.array([pdata.GetCellData().GetArray("LV_strain").GetTuple(k_cell)[0] for k_cell in range(pdata.GetNumberOfCells())])
	mean_Ecc.append(np.mean(Ecc))
	Ell = np.array([pdata.GetCellData().GetArray("LV_strain").GetTuple(k_cell)[1] for k_cell in range(pdata.GetNumberOfCells())])
	mean_Ell.append(np.mean(Ell))
	Err = np.array([pdata.GetCellData().GetArray("LV_strain").GetTuple(k_cell)[2] for k_cell in range(pdata.GetNumberOfCells())])
	mean_Err.append(np.mean(Err))
			

	# Write Ugrid
	vtk_py.writeXMLPData(pdata, newfilename)


plt.plot(np.arange(0, ntpt), mean_Ecc, label="Ecc")
plt.plot(np.arange(0, ntpt), mean_Ell, label="Ell")
np.savetxt("strainv85.txt", np.c_[mean_Ecc, mean_Ell])
#plt.plot(np.arange(0, ntpt), mean_Err, label="Err")
plt.xlabel("Time point")
plt.ylabel("Strain")
plt.legend()
plt.savefig("tempv85.png")

#plt

	
	

