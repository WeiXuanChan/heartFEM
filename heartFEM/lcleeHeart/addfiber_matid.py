from dolfin import *
import vtk as vtk
from heartFEM.lcleeHeart.vtk_py import *

def addfiber_matid(mesh, V, casename, endo_angle, epi_angle,  infarct_endo_angle, infarct_epi_angle, casedir, matid, isepiflip, isendoflip):


		fiberV = Function(V)
		sheetV = Function(V)
		sheetnormV = Function(V)

	#endofilename = casedir + 'endocardium.stl'
	#epifilename = casedir + 'epicardium.stl'

	#endofilename = casename + '_endo.stl'
	#epifilename = casename + '_epi.stl'

	#pdata_endo = readSTL(endofilename)
	#pdata_epi = readSTL(epifilename)

		ugrid=vtk_py.convertXMLMeshToUGrid(mesh)

		pdata = vtk_py.convertUGridtoPdata(ugrid)
		C = vtk_py.getcentroid(pdata)
		ztop = pdata.GetBounds()[5]
		C = [C[0], C[1], ztop-0.05]
		clippedheart = vtk_py.clipheart(pdata, C, [0,0,1], True)
		epi, endo= vtk_py.splitDomainBetweenEndoAndEpi(clippedheart)

		cleanepipdata = vtk.vtkCleanPolyData()
		if (vtk.vtkVersion.GetVTKMajorVersion() >= 6):
				cleanepipdata.SetInputData(epi)
		else:
				cleanepipdata.SetInput(epi)
		cleanepipdata.Update()
		pdata_epi = cleanepipdata.GetOutput()

		cleanendopdata = vtk.vtkCleanPolyData()
		if (vtk.vtkVersion.GetVTKMajorVersion() >= 6):
				cleanendopdata.SetInputData(endo)
		else:
				cleanendopdata.SetInput(endo)
		cleanendopdata.Update()
		pdata_endo = cleanendopdata.GetOutput()

		L_epi = pdata_epi.GetBounds()[5]  -  pdata_epi.GetBounds()[4]
		L_endo = pdata_endo.GetBounds()[5] - pdata_endo.GetBounds()[4]

		if(L_endo > L_epi):
			pdata_epi = temp
			pdata_epi = pdata_endo
			pdata_endo = temp
		

	# Quad points
		gdim = mesh.geometry().dim()
		xdofmap = V.sub(0).dofmap().dofs()
		ydofmap = V.sub(1).dofmap().dofs()
		zdofmap = V.sub(2).dofmap().dofs()

		if(dolfin.dolfin_version() != '1.6.0'):
			xq = V.tabulate_dof_coordinates().reshape((-1, gdim))
			xq0 = xq[xdofmap]  
		else:
			xq = V.dofmap().tabulate_all_coordinates(mesh).reshape((-1, gdim))
			xq0 = xq[xdofmap]  

	# Create an unstructured grid of Gauss Points
			points = vtk.vtkPoints()
			ugrid = vtk.vtkUnstructuredGrid()
			cnt = 0;
		for pt in xq0:
			points.InsertNextPoint([pt[0], pt[1], pt[2]])
			cnt += 1

			ugrid.SetPoints(points)

			CreateVertexFromPoint(ugrid)
			addLocalProlateSpheroidalDirections(ugrid, pdata_endo, pdata_epi, type_of_support="cell", epiflip=isepiflip, endoflip=isendoflip)
			addLocalFiberOrientation_infarct(ugrid, endo_angle, epi_angle, infarct_endo_angle, infarct_epi_angle, matid)

			fiber_vector =  ugrid.GetCellData().GetArray("fiber vectors")
			sheet_vector =  ugrid.GetCellData().GetArray("sheet vectors")
			sheetnorm_vector =  ugrid.GetCellData().GetArray("sheet normal vectors")

			cnt = 0
		for pt in xq0:

			fvec = fiber_vector.GetTuple(cnt)
			svec = sheet_vector.GetTuple(cnt)
			nvec = sheetnorm_vector.GetTuple(cnt)

			fvecnorm = sqrt(fvec[0]**2 + fvec[1]**2 + fvec[2]**2)
			svecnorm = sqrt(svec[0]**2 + svec[1]**2 + svec[2]**2)
			nvecnorm = sqrt(nvec[0]**2 + nvec[1]**2 + nvec[2]**2)

		if(abs(fvecnorm - 1.0) > 1e-7 or  abs(svecnorm - 1.0) > 1e-6 or abs(nvecnorm - 1.0) > 1e-7):
			print (fvecnorm)
			print (svecnorm)
			print (nvecnorm)

		#print(xdofmap[cnt], ydofmap[cnt], zdofmap[cnt])
		fiberV.vector()[xdofmap[cnt]] = fvec[0]; fiberV.vector()[ydofmap[cnt]] = fvec[1]; fiberV.vector()[zdofmap[cnt]] = fvec[2];
		sheetV.vector()[xdofmap[cnt]] = svec[0]; sheetV.vector()[ydofmap[cnt]] = svec[1]; sheetV.vector()[zdofmap[cnt]] = svec[2];
		sheetnormV.vector()[xdofmap[cnt]] = nvec[0]; sheetnormV.vector()[ydofmap[cnt]] = nvec[1]; sheetnormV.vector()[zdofmap[cnt]] = nvec[2];

		cnt += 1


		writeXMLUGrid(ugrid, "test.vtu")

		return fiberV, sheetV, sheetnormV


