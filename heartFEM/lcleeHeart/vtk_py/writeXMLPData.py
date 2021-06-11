########################################################################

import vtk

########################################################################

def writeXMLPData(pdata, pdata_file_name, verbose=True):
    if (verbose): print ('*** writeXMLPData ***')

    pdata_writer = vtk.vtkXMLPolyDataWriter()
    pdata_writer.SetFileName(pdata_file_name)
    if (vtk.vtkVersion.GetVTKMajorVersion() >= 6):
        pdata_writer.SetInputData(pdata)
    else:
        pdata_writer.SetInput(pdata)
    pdata_writer.Update()
    pdata_writer.Write()
