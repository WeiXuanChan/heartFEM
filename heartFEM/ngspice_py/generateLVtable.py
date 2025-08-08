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
  Author: w.x.chan@gmail.com         09Nov2021           - v4.0.0
                                                              -change time unit to ms
'''
########################################################################
_version='4.0.0'
import logging
logger = logging.getLogger(__name__)


import sys
import vtk
import os
import inspect
import numpy as np
from scipy import interpolate
import scipy.optimize as opt
########################################################################
def exponential2D_log(xdata_tuple, amplitude_x,amplitude_y, sigma_x, sigma_y, offset):
    (x, y) = xdata_tuple 
    ex =   amplitude_x**2.*np.exp(sigma_x**2.*x)
    ey =   amplitude_y**2.*np.exp(sigma_y**2.*y)   
    g = np.log(offset**2. + ex+ey)                                  
    return g
def exponential1D_log(x, amplitude_x, sigma_x, offset):
    ex =   amplitude_x**2.*np.exp(sigma_x**2.*x)
    g = np.log(offset**2. + ex)                           
    return g
def generatetable(cavity,casename,period,timetopeak=None,verbose=True,stackaddstr=None,loading_casename=None,scale_maximum_fiber_tension=1.,adjust_relaxedPressure=True,extendvoldata=1.):
    #period same units as timeSpace
    logger.info('*** generate'+cavity.upper()+'table ***')
    extend_num_of_points=10
    if loading_casename is None:
        loading_casename=casename
    period=period
    if stackaddstr is None:
        stackaddstr=['']
    LVtablefile = casename + "_"+cavity.lower()+"cirtable.txt"
    
    ytime=np.loadtxt(loading_casename+"_Press_timeSpace.txt")
    
    ytime=np.round(ytime,decimals=4)
    
    xvol=[]
    zvol=[]
    for addstr in stackaddstr:
        xvol_temp=np.loadtxt(loading_casename+"_Press_volumeSpace"+addstr+".txt")
        if len(xvol_temp.shape)>1:
            xvol.append(xvol_temp[0])
            zvol.append(xvol_temp[1])
        else:
            xvol.append(xvol_temp)
    if len(xvol)>1:
        xvol=np.concatenate(xvol,axis=0)
        if len(zvol)>0:
            zvol=np.concatenate(zvol,axis=0)
    else:
        xvol=np.array(xvol[0])
        if len(zvol)>0:
            zvol=np.array(zvol[0])
    xvol, sortVolume=np.unique(xvol, return_index=True)
    if len(zvol)>0:
        zvol, zsortVolume=np.unique(zvol, return_index=True)
    if extendvoldata>0:
        xvol_extended=np.concatenate((xvol[0:1]*(1.-extendvoldata*np.linspace(1,0,num=extend_num_of_points,endpoint=False)),xvol,xvol[-1:]*(1.+extendvoldata*np.linspace(1,0,num=extend_num_of_points,endpoint=False)[::-1])),axis=0)
        if len(zvol)>0:
            zvol_extended=np.concatenate((zvol[0:1]*(1.-extendvoldata*np.linspace(1,0,num=extend_num_of_points,endpoint=False)),zvol,zvol[-1:]*(1.+extendvoldata*np.linspace(1,0,num=extend_num_of_points,endpoint=False)[::-1])),axis=0)
    else:
        xvol_extended=xvol
        if len(zvol)>0:
            zvol_extended=zvol
    datatable=[]
    for addstr in stackaddstr:
        if len(zvol)>0:
            datatable.append(np.loadtxt(loading_casename+"_Press_VolTime_"+cavity.upper()+addstr+".txt").reshape((len(zvol),len(xvol),-1)))
        else:
            datatable.append(np.loadtxt(loading_casename+"_Press_VolTime_"+cavity.upper()+addstr+".txt"))
    if len(datatable)>1:
        datatable=np.concatenate(datatable,axis=-2)
    else:
        datatable=np.array(datatable[0])
    if adjust_relaxedPressure:
        if len(zvol)>0:
            datatable=datatable.reshape((len(zvol)*len(xvol),-1))
        for volN in range(datatable.shape[0]):
            maxInd=np.argmax(datatable[volN])
            lastZeroInd=np.nonzero(datatable[volN][maxInd:]==0)[0][0]-1+maxInd
            if (datatable[volN][maxInd]-datatable[volN][lastZeroInd])==0:
                maxInd=np.argmin(datatable[volN])
                lastZeroInd=np.nonzero(datatable[volN][maxInd:]==0)[0][0]-1+maxInd
            datatable[volN][maxInd:(lastZeroInd+1)]=datatable[volN][maxInd]*(datatable[volN][maxInd:(lastZeroInd+1)]-datatable[volN][lastZeroInd])/(datatable[volN][maxInd]-datatable[volN][lastZeroInd])
        if len(zvol)>0:
            datatable=datatable.reshape((len(zvol),len(xvol),-1))
    if len(zvol)>0:
        new_datatable=[]
        for n in range(len(zvol)):
            datatable[n]=datatable[n][sortVolume]
            new_datatable.append(datatable[n].T*scale_maximum_fiber_tension)
        new_datatable=np.array(new_datatable)[zsortVolume]
        if extendvoldata>0:
            new_datatable=np.pad(new_datatable,((extend_num_of_points,extend_num_of_points),(0,0),(extend_num_of_points,extend_num_of_points)),mode='edge')
        new_datatable=new_datatable.reshape((-1,len(xvol_extended)))
        datatable=datatable[zsortVolume]
        datatable=datatable.reshape((len(zvol)*len(xvol),-1))
    else:
        datatable=datatable[sortVolume]
        new_datatable=datatable.T*scale_maximum_fiber_tension 
        if extendvoldata>0:
            new_datatable=np.pad(new_datatable,((0,0),(extend_num_of_points,extend_num_of_points)),mode='edge')
    np.savetxt(LVtablefile,new_datatable,fmt='%.9e')
    
    with open(LVtablefile,'r') as f:
        lines=f.readlines()
    if len(zvol)>0:
        newline=['*table source\n',
                 '*number of columns (x)\n',
                 str(len(xvol_extended))+'\n',
                 '*number of rows (y)\n',
                 str(len(ytime))+'\n',
                 '*number of blocks (z)\n',
                 str(len(zvol_extended))+'\n',
                 '*x horizontal (column) address values (real numbers)\n',
                 ' '.join(xvol_extended.astype(str))+'\n',
                 '*y vertical (row) address values (real numbers)\n',
                 ' '.join(ytime.astype(str))+'\n',
                 '*z block address values (real numbers)\n',
                 ' '.join(zvol_extended.astype(str))+'\n',
                 '*table with output data (horizontally addressed by x, vertically by y)\n']
    else:
        newline=['*table source\n',
                 '*number of columns (x)\n',
                 str(len(xvol_extended))+'\n',
                 '*number of rows (y)\n',
                 str(len(ytime))+'\n',
                 '*x horizontal (column) address values (real numbers)\n',
                 ' '.join(xvol_extended.astype(str))+'\n',
                 '*y vertical (row) address values (real numbers)\n',
                 ' '.join(ytime.astype(str))+'\n',
                 '*table with output data (horizontally addressed by x, vertically by y)\n']
    lines=newline+lines
    with open(LVtablefile,'w') as f:
        f.writelines(lines)
    if len(zvol)>0:
        new_zvol=np.array(zvol_extended)
    else:
        new_zvol=np.array([-1,0,1])
    if timetopeak is not None:
        timetopeak=timetopeak
        LVtabletrfile = casename + "_"+cavity.lower()+"cirtabletr.txt"
        trdata=[]
        for Nvol in range(len(datatable)):
            spl = interpolate.splrep(ytime,datatable[Nvol])
            temp_maxpress = interpolate.splev(np.array([timetopeak]), spl)[0]
            temp_datatable=datatable[Nvol][ytime>timetopeak]
            if temp_maxpress<0:
                temp_maxpress=temp_maxpress*-1.
                temp_datatable=temp_datatable*-1.
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
        if len(zvol)>0:
            trdata=trdata.reshape((len(zvol),-1)).T
            if extendvoldata>0:
                trdata=np.pad(trdata,extend_num_of_points,mode='edge')
        else:
            trdata=np.concatenate((trdata.reshape((-1,1)),trdata.reshape((-1,1)),trdata.reshape((-1,1))),axis=1)
            if extendvoldata>0:
                trdata=np.pad(trdata,((extend_num_of_points,extend_num_of_points),(0,0)),mode='edge')
        np.savetxt(LVtabletrfile,trdata,fmt='%.9e')
        with open(LVtabletrfile,'r') as f:
            lines=f.readlines()
        newline=['*table source\n',
             '*number of columns (x)\n',
             str(trdata.shape[1])+'\n',
             '*number of rows (y)\n',
             str(len(xvol_extended))+'\n',
             '*x horizontal (column) address values (real numbers)\n',
             ' '.join(new_zvol.astype(str))+'\n',
             '*y vertical (row) address values (real numbers)\n',
             ' '.join(xvol_extended.astype(str))+'\n',
             '*table with output data (horizontally addressed by x, vertically by y)\n']
        lines=newline+lines
        with open(LVtabletrfile,'w') as f:
            f.writelines(lines)
            
    LVtablebasefile = casename + "_"+cavity.lower()+"cirtablebase.txt"
    
    datatable=[]
    for addstr in stackaddstr:
        if len(zvol)>0:
            datatable.append(np.loadtxt(loading_casename+"_Press_VolTime_"+cavity.upper()+"_base"+addstr+".txt").reshape((len(zvol),len(xvol))))
        else:
            datatable.append(np.loadtxt(loading_casename+"_Press_VolTime_"+cavity.upper()+"_base"+addstr+".txt"))
    if len(datatable)>1:
        datatable=np.concatenate(datatable,axis=-1)
    else:
        datatable=np.array(datatable[0])
    if len(zvol)>0:
        if extendvoldata>0:
            temp_smallx=np.zeros((extend_num_of_points,datatable.shape[0]))
            temp_largex=np.zeros((extend_num_of_points,datatable.shape[0]))
        for n in range(len(zvol)):
            datatable[n]=datatable[n][sortVolume]
            if extendvoldata>0:
                spl = interpolate.splrep(xvol,datatable[n])
                temp_smallx[0:extend_num_of_points,n] = interpolate.splev(xvol_extended[:extend_num_of_points], spl)[0:extend_num_of_points]
                temp_largex[0:extend_num_of_points,n] = interpolate.splev(xvol_extended[-extend_num_of_points:], spl)[0:extend_num_of_points]
                #spl=np.poly1d(np.polyfit(xvol[:2],datatable[n][:2],1))
                #temp_smallx[0:extend_num_of_points,n] = spl(xvol_extended[:extend_num_of_points])[0:extend_num_of_points]
                #spl=np.poly1d(np.polyfit(xvol[-2:],datatable[n][-2:],1))
                #temp_largex[0:extend_num_of_points,n] = spl(xvol_extended[-extend_num_of_points:])[0:extend_num_of_points]
        datatable=datatable[zsortVolume].T
        if extendvoldata>0:
            datatable=np.concatenate((temp_smallx,datatable,temp_largex),axis=0)
            temp_smallz=np.zeros((datatable.shape[0],extend_num_of_points))
            temp_largez=np.zeros((datatable.shape[0],extend_num_of_points))
            for n in range(datatable.shape[0]):
                spl = interpolate.splrep(zvol,datatable[n])
                temp_smallz[n,0:extend_num_of_points] = interpolate.splev(zvol_extended[:extend_num_of_points], spl)[0:extend_num_of_points]
                temp_largez[n,0:extend_num_of_points] = interpolate.splev(zvol_extended[-extend_num_of_points:], spl)[0:extend_num_of_points]
                #spl=np.poly1d(np.polyfit(zvol[:2],datatable[n][:2],1))
                #temp_smallz[n,0:extend_num_of_points] = spl(zvol_extended[:extend_num_of_points])[0:extend_num_of_points]
                #spl=np.poly1d(np.polyfit(zvol[-2:],datatable[n][-2:],1))
                #temp_largez[n,0:extend_num_of_points] = spl(zvol_extended[-extend_num_of_points:])[0:extend_num_of_points]
            datatable=np.concatenate((temp_smallz,datatable,temp_largez),axis=1)
            #popt, pcov = opt.curve_fit(exponential1D_log, zvol, np.log(datatable[0]-datatable[0].min()+1.), p0 = (1., 1., 1.))
            #amplitude_z=popt[0]
            #sigma_z=popt[1]
            #offset=popt[2]
            #popt, pcov = opt.curve_fit(exponential1D_log, xvol, np.log(datatable[:,0]-datatable[:,0].min()+1.), p0 = (1., 1., 1.))
            #amplitude_x=popt[0]
            #sigma_x=popt[1]
            #offset=offset+popt[2]
            #temp_zvol,temp_xvol=np.meshgrid(zvol,xvol)
            #initial_guess= (amplitude_z,amplitude_x, sigma_z, sigma_x, offset)
            #popt, pcov = opt.curve_fit(exponential2D_log, (temp_zvol.reshape(-1),temp_xvol.reshape(-1)), np.log(datatable-datatable.min()+1.).reshape(-1), p0 = initial_guess)
            #temp_zvol,temp_xvol=np.meshgrid(zvol_extended,xvol_extended)
            #temp_datatable = (np.exp(exponential2D_log((temp_zvol, temp_xvol), *popt))+datatable.min()-1.)
            #temp_datatable[extend_num_of_points:-extend_num_of_points,extend_num_of_points:-extend_num_of_points]=datatable
            #datatable=np.array(temp_datatable)
    else:
        datatable=datatable[sortVolume]
        if extendvoldata>0:
            spl = interpolate.splrep(xvol,datatable)
            temp_smallx = interpolate.splev(xvol_extended[:extend_num_of_points], spl)[0:extend_num_of_points]
            temp_largex = interpolate.splev(xvol_extended[-extend_num_of_points:], spl)[0:extend_num_of_points]
            #spl=np.poly1d(np.polyfit(xvol[:2],datatable[:2],1))
            #temp_smallx = spl(xvol_extended[:extend_num_of_points])[0:extend_num_of_points]
            #spl=np.poly1d(np.polyfit(xvol[-2:],datatable[-2:],1))
            #temp_largex = spl(xvol_extended[-extend_num_of_points:])[0:extend_num_of_points]
            datatable=np.concatenate((temp_smallx,datatable,temp_largex),axis=0)
            #popt, pcov = opt.curve_fit(exponential1D_log, xvol, np.log(datatable-datatable.min()+1.), p0 = (1., 1., 1.))
            #temp_datatable = np.exp(exponential1D_log(xvol_extended, *popt))+datatable.min()-1.
            #temp_datatable[extend_num_of_points:-extend_num_of_points]=datatable
            #datatable=np.array(temp_datatable)
        datatable=np.concatenate((datatable.reshape((-1,1)),datatable.reshape((-1,1)),datatable.reshape((-1,1))),axis=1)
    np.savetxt(LVtablebasefile,datatable,fmt='%.9e')
    with open(LVtablebasefile,'r') as f:
        lines=f.readlines()
    newline=['*table source\n',
             '*number of columns (x)\n',
             str(datatable.shape[1])+'\n',
             '*number of rows (y)\n',
             str(len(xvol_extended))+'\n',
             '*x horizontal (column) address values (real numbers)\n',
             ' '.join(new_zvol.astype(str))+'\n',
             '*y vertical (row) address values (real numbers)\n',
             ' '.join(xvol_extended.astype(str))+'\n',
             '*table with output data (horizontally addressed by x, vertically by y)\n']
    lines=newline+lines
    with open(LVtablebasefile,'w') as f:
        f.writelines(lines)
        
def generateLAtable(casename,period,timetopeak=None,verbose=True,stackaddstr=None,loading_casename=None,scale_maximum_fiber_tension=1.,adjust_relaxedPressure=True):
    generatetable("LA",casename,period,timetopeak=timetopeak,verbose=verbose,stackaddstr=stackaddstr,loading_casename=loading_casename,scale_maximum_fiber_tension=scale_maximum_fiber_tension,adjust_relaxedPressure=adjust_relaxedPressure)
def generateRAtable(casename,period,timetopeak=None,verbose=True,stackaddstr=None,loading_casename=None,scale_maximum_fiber_tension=1.,adjust_relaxedPressure=True):
    generatetable("RA",casename,period,timetopeak=timetopeak,verbose=verbose,stackaddstr=stackaddstr,loading_casename=loading_casename,scale_maximum_fiber_tension=scale_maximum_fiber_tension,adjust_relaxedPressure=adjust_relaxedPressure)
def generateLVtable(casename,period,timetopeak=None,verbose=True,stackaddstr=None,loading_casename=None,scale_maximum_fiber_tension=1.,adjust_relaxedPressure=True):
    generatetable("LV",casename,period,timetopeak=timetopeak,verbose=verbose,stackaddstr=stackaddstr,loading_casename=loading_casename,scale_maximum_fiber_tension=scale_maximum_fiber_tension,adjust_relaxedPressure=adjust_relaxedPressure)
def generateRVtable(casename,period,timetopeak=None,verbose=True,stackaddstr=None,loading_casename=None,scale_maximum_fiber_tension=1.,adjust_relaxedPressure=True):
    generatetable("RV",casename,period,timetopeak=timetopeak,verbose=verbose,stackaddstr=stackaddstr,loading_casename=loading_casename,scale_maximum_fiber_tension=scale_maximum_fiber_tension,adjust_relaxedPressure=adjust_relaxedPressure)