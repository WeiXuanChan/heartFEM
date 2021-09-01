'''
File: generateLVtable.py
Description: creates a table to be used as LV source 2dtable
History:
    Date    Programmer SAR# - Description
    ---------- ---------- ----------------------------
  Author: w.x.chan@gmail.com         08MAR2021           - Created
  Author: w.x.chan@gmail.com         08MAR2021           - v1.0.0
  Author: w.x.chan@gmail.com         13APR2021           - v2.0.0
  Author: w.x.chan@gmail.com         21APR2021           - v2.1.0
  Author: w.x.chan@gmail.com         21APR2021           - v2.1.1
                                                            -debug when stackaddstr is not sorted and/or not unique
  Author: w.x.chan@gmail.com         12MAY2021           - v2.3.4
                                                            -debug read LVtablefile before savetxt
  Author: w.x.chan@gmail.com         12Aug2021           - v3.5.0
                                                            -added 'Windkessel_scale_T0_LV'
'''
########################################################################
_version='3.5.0'
import logging
logger = logging.getLogger(__name__)


import sys
import vtk
import os
import inspect
from heartFEM import ngspice_py
import numpy as np
from scipy import interpolate
########################################################################
Components=['lv','la','rv','ra','aa','ao1','ao2','ao3','ao4','br','ca','ub','he','inte','ivc','kid','leg','lung','pa1','pa2','plac','svc','uv']
linkComponents=['aaao1','ao1ao2','ao2ao3','ao3ao4','pa1pa2','pa2lung','da','ao1ca','cabr','brsvc','ao1ub','ubsvc','ao3he','ao3inte','intehe','ao3kid','kidivc','ao4plac','placuv','ao4leg','legivc','uvhe','heivc','dv','svcra','ivcra','lungla','fo','raravalv','rvrvvalv','lvlvvalv','lalavalv']
keyToFill=[]
for comp in Components:
    keyToFill.append(comp+'c')

for comp in linkComponents:
    keyToFill.append(comp+'r')
    keyToFill.append(comp+'l')
    keyToFill.append(comp+'k')
    keyToFill.append(comp+'b')
    
suffixDict={4:'T  ',3:'g  ',2:'meg',1:'k  ',0:' ',-1:'m  ',-2:'u  ',-3:'n  ',-4:'p  ',-5:'f  '}

def generateLVtable(casename,period,timetopeak=None,verbose=True,stackaddstr=None,loading_casename=None,scale_T0_LV=1.):
    #period same units as timeSpace
    logger.info('*** generateLVtable ***')
    if loading_casename is None:
        loading_casename=casename
    period=period/1000.
    if stackaddstr is None:
        stackaddstr=['']
    LVtablefile = casename + "_lvcirtable.txt"
    
    
    ytime=np.loadtxt(loading_casename+"_Press_timeSpace.txt")/1000.
    
    ytime=np.round(ytime,decimals=4)
    
    xvol=[]
    for addstr in stackaddstr:
        xvol.append(np.loadtxt(loading_casename+"_Press_volumeSpace"+addstr+".txt"))
    if len(xvol)>1:
        xvol=np.concatenate(xvol,axis=0)
    else:
        xvol=np.array(xvol[0])
    xvol, sortVolume=np.unique(xvol, return_index=True)
    
    datatable=[]
    for addstr in stackaddstr:
        datatable.append(np.loadtxt(loading_casename+"_Press_VolTime"+addstr+".txt"))
    if len(datatable)>1:
        datatable=np.concatenate(datatable,axis=0)
    else:
        datatable=np.array(datatable[0])
    datatable=datatable[sortVolume]
    
    datatable=datatable.T*scale_T0_LV
    np.savetxt(LVtablefile,datatable,fmt='%.9e')

    with open(LVtablefile,'r') as f:
        lines=f.readlines()
    
    newline=['*table source\n',
             '*number of columns (x)\n',
             str(len(xvol))+'\n',
             '*number of rows (y)\n',
             str(len(ytime))+'\n',
             '*x horizontal (column) address values (real numbers)\n',
             ' '.join(xvol.astype(str))+'\n',
             '*y vertical (row) address values (real numbers)\n',
             ' '.join(ytime.astype(str))+'\n',
             '*table with output data (horizontally addressed by x, vertically by y)\n']
    lines=newline+lines
    with open(LVtablefile,'w') as f:
        f.writelines(lines)
    
    if timetopeak is not None:
        timetopeak=timetopeak/1000.
        LVtabletrfile = casename + "_lvcirtabletr.txt"
        trdata=[]
        for Nvol in range(len(xvol)):
            spl = interpolate.splrep(ytime,datatable[:,Nvol])
            temp_maxpress = interpolate.splev(np.array([timetopeak]), spl)[0]
            temp_datatable=datatable[:,Nvol][ytime>timetopeak]
            temp_ytime=ytime[ytime>timetopeak]
            spl = interpolate.splrep(temp_ytime,temp_datatable)
            tryTimeInd=np.argwhere(temp_datatable<temp_maxpress/2.)[0,0]
            tryTime=temp_ytime[tryTimeInd-1]
            tryTimemin=temp_ytime[tryTimeInd-1]
            tryTimemax=temp_ytime[tryTimeInd]
            
            while abs(tryTimemax-tryTimemin)>10**-6:
                tryTime_pressure=interpolate.splev(np.array([tryTime]), spl)[0]
                if tryTime_pressure>temp_maxpress/2.:
                    tryTimemin=tryTime
                elif tryTime_pressure<temp_maxpress/2.:
                    tryTimemax=tryTime
                else:
                    tryTimemin=tryTime
                    tryTimemax=tryTime
                    break
                tryTime=(tryTimemin+tryTimemax)/2.
            trdata.append((tryTime-timetopeak)*2)
        trdata=np.array(trdata)
        trdata=np.concatenate((trdata.reshape((-1,1)),trdata.reshape((-1,1)),trdata.reshape((-1,1))),axis=1)
        np.savetxt(LVtabletrfile,trdata,fmt='%.9e')
        with open(LVtabletrfile,'r') as f:
            lines=f.readlines()
        newline=['*table source\n',
             '*number of columns (x)\n',
             str(3)+'\n',
             '*number of rows (y)\n',
             str(len(xvol))+'\n',
             '*x horizontal (column) address values (real numbers)\n',
             ' '.join(np.array([-1,0,1]).astype(str))+'\n',
             '*y vertical (row) address values (real numbers)\n',
             ' '.join(xvol.astype(str))+'\n',
             '*table with output data (horizontally addressed by x, vertically by y)\n']
        lines=newline+lines
        with open(LVtabletrfile,'w') as f:
            f.writelines(lines)
            
    LVtablebasefile = casename + "_lvcirtablebase.txt"
    datatable=[]
    for addstr in stackaddstr:
        datatable.append(np.loadtxt(loading_casename+"_Press_VolTime_base"+addstr+".txt"))
    if len(datatable)>1:
        datatable=np.concatenate(datatable,axis=0)
    else:
        datatable=np.array(datatable[0])
    datatable=datatable[::-1]
    datatable=np.concatenate((datatable.reshape((-1,1)),datatable.reshape((-1,1)),datatable.reshape((-1,1))),axis=1)
    np.savetxt(LVtablebasefile,datatable,fmt='%.9e')
    with open(LVtablebasefile,'r') as f:
        lines=f.readlines()
    newline=['*table source\n',
             '*number of columns (x)\n',
             str(3)+'\n',
             '*number of rows (y)\n',
             str(len(xvol))+'\n',
             '*x horizontal (column) address values (real numbers)\n',
             ' '.join(np.array([-1,0,1]).astype(str))+'\n',
             '*y vertical (row) address values (real numbers)\n',
             ' '.join(xvol.astype(str))+'\n',
             '*table with output data (horizontally addressed by x, vertically by y)\n']
    lines=newline+lines
    with open(LVtablebasefile,'w') as f:
        f.writelines(lines)
        
