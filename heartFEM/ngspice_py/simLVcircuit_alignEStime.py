'''
File: simLVcircuit_alignEStime.py
Description: simulate circuit with ngspice
History:
    Date    Programmer SAR# - Description
    ---------- ---------- ----------------------------
  Author: w.x.chan@gmail.com         08MAR2021           - Created
  Author: w.x.chan@gmail.com         13APR2021           - v2.0.0
  Author: w.x.chan@gmail.com         21APR2021           - v2.1.0
  Author: w.x.chan@gmail.com         09JUN2021           - v3.0.0
                                                              -adjust EStime to within 0.1% od minimum volume
                                                              -adjust try_initLVvol=initLVvol*1.05
  Author: w.x.chan@gmail.com         28Sep2021           - v3.5.1
                                                              -add 2nd guess for adjust EDVol
'''
########################################################################
_version='3.5.1'
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

suffixDict={4:'T  ',3:'g  ',2:'meg',1:'k  ',0:' ',-1:'m  ',-2:'u  ',-3:'n  ',-4:'p  ',-5:'f  '}

def simLVcircuit_alignEStime(casename,stopTime,lvufile,period,targetEStime,init_timetopeaktension,try_timetopeaktension=None,lvinputvar='V',initLAvol=0,initRAvol=0,initLVvol=0,initRVvol=0,vla0=None,vra0=None,init_file=None,init_time=None,iterationNumber=100,verbose=True):

    if try_timetopeaktension is None:
        try_timetopeaktension=init_timetopeaktension
        adj_timetopeaktension=0.1*init_timetopeaktension
    elif isinstance(try_timetopeaktension,(int,float)):
        adj_timetopeaktension=0.1*init_timetopeaktension
    else:
        adj_timetopeaktension=try_timetopeaktension[1]
        try_timetopeaktension=try_timetopeaktension[0]
    tune_timetopeaktension=None
    last_EStime=0
    for n in range(iterationNumber):
        logger.info("    Trying timetopeak="+repr(try_timetopeaktension))
        case_dir,lvufilename = os.path.split(lvufile)
        timetopeak_from_to=[init_timetopeaktension,try_timetopeaktension]
        ngspice_py.simLVcircuit(casename,stopTime,lvufile,lvinputvar=lvinputvar,initLAvol=initLAvol,initRAvol=initRAvol,initLVvol=initLVvol,initRVvol=initRVvol,vla0=vla0,vra0=vra0,init_file=init_file,init_time=init_time,timetopeak_from_to=timetopeak_from_to,verbose=verbose)
        
        cir_results=np.loadtxt(case_dir+'/'+'circuit_results.txt',skiprows=1)[:,2:4]
        cir_results[:,0]*=1000. #set to ms
        while cir_results[-1,0]>=(2.*period):
            cir_results[:,0]-=period
        cir_results=cir_results[cir_results[:,0]>=0]
        cir_results=cir_results[cir_results[:,0]<period]
        minvolindex=np.argmin(cir_results[:,1])
        currentEStime=cir_results[np.nonzero(cir_results[:(minvolindex+1),1]<=cir_results[minvolindex,1]*1.001)[0][0],0]
        logger.info("      Current EStime="+repr(currentEStime)+', target EStime='+repr(targetEStime))
        if abs(currentEStime-targetEStime)<10 or adj_timetopeaktension<10:
            break
        elif abs(last_EStime-currentEStime)<5:
            return (try_timetopeaktension-tune_timetopeaktension*adj_timetopeaktension,adj_timetopeaktension)
        elif currentEStime>targetEStime:
            if tune_timetopeaktension is None:
                tune_timetopeaktension=-1
            elif tune_timetopeaktension>0:
                tune_timetopeaktension=-1
                adj_timetopeaktension=adj_timetopeaktension*0.33
            try_timetopeaktension=try_timetopeaktension+tune_timetopeaktension*adj_timetopeaktension
        elif currentEStime<targetEStime:
            if tune_timetopeaktension is None:
                tune_timetopeaktension=1
            elif tune_timetopeaktension<0:
                tune_timetopeaktension=1
                adj_timetopeaktension=adj_timetopeaktension*0.33
            try_timetopeaktension=try_timetopeaktension+tune_timetopeaktension*adj_timetopeaktension
        last_EStime=currentEStime
    
    return (try_timetopeaktension,adj_timetopeaktension)

def simLVcircuit_align_EDvol_and_EStime(casename,stopTime,lvufile,period,targetEStime,init_timetopeaktension,try_initLVvol=None,lvinputvar='V',initLAvol=0,initRAvol=0,initLVvol=0,initRVvol=0,vla0=None,vra0=None,init_file=None,init_time=None,iterationNumber=100,verbose=True):
    targetinitLVvol=initLVvol
    if try_initLVvol is None:
        try_initLVvol=initLVvol*1.05
        adj_try_initLVvol=0.05*targetinitLVvol
    elif isinstance(try_initLVvol,(int,float)):
        adj_try_initLVvol=0.05*targetinitLVvol
    else:
        adj_try_initLVvol=try_initLVvol[1]
        try_initLVvol=try_initLVvol[0]

    tune_initLVvol=None
    last_currentinitLVvol=0
    for n in range(iterationNumber):
        logger.info("Trying initLVvol="+repr(try_initLVvol))
        if targetEStime is None:
            ngspice_py.simLVcircuit(casename,stopTime,lvufile,lvinputvar=lvinputvar,initLAvol=initLAvol,initRAvol=initRAvol,initLVvol=try_initLVvol,initRVvol=initRVvol,vla0=vla0,vra0=vra0,init_file=init_file,init_time=init_time,verbose=verbose)
        else:
            if n==0:
                try_timetopeaktension=None
            else:
                try_timetopeaktension=list(try_timetopeaktension)
                try_timetopeaktension[1]=try_timetopeaktension[1]*3.
            try_timetopeaktension=simLVcircuit_alignEStime(casename,stopTime,lvufile,period,targetEStime,init_timetopeaktension,try_timetopeaktension=try_timetopeaktension,lvinputvar=lvinputvar,initLAvol=initLAvol,initRAvol=initRAvol,initLVvol=try_initLVvol,initRVvol=initRVvol,vla0=vla0,vra0=vra0,init_file=init_file,init_time=init_time,iterationNumber=iterationNumber,verbose=verbose)
        case_dir,outfilename = os.path.split(casename)
        cir_results=np.loadtxt(case_dir+'/'+'circuit_results.txt',skiprows=1)[:,2:4]
        cir_results[:,0]*=1000. #set to ms
        while cir_results[-1,0]>=(2.*period):
            cir_results[:,0]-=period
        cir_results=cir_results[cir_results[:,0]>=0]
        cir_results=cir_results[cir_results[:,0]<period]
        currentinitLVvol=cir_results[0,1]
        logger.info("  Current EDvol="+repr(currentinitLVvol)+', target EDvol='+repr(targetinitLVvol))
        if abs(currentinitLVvol-targetinitLVvol)<(targetinitLVvol*10.**-4.) or adj_try_initLVvol<(targetinitLVvol*10.**-4.):
            break
        elif abs(last_currentinitLVvol-currentinitLVvol)<(targetinitLVvol*10.**-4.):
            return (try_initLVvol-tune_initLVvol*adj_try_initLVvol,adj_try_initLVvol)
        elif currentinitLVvol>targetinitLVvol:
            if tune_initLVvol is None:
                tune_initLVvol=-1
                try_initLVvol=try_initLVvol*targetinitLVvol/currentinitLVvol
            elif tune_initLVvol>0:
                tune_initLVvol=-1
                adj_try_initLVvol=adj_try_initLVvol*0.33
                try_initLVvol=try_initLVvol+tune_initLVvol*adj_try_initLVvol
        elif currentinitLVvol<targetinitLVvol:
            if tune_initLVvol is None:
                tune_initLVvol=1
                try_initLVvol=try_initLVvol*targetinitLVvol/currentinitLVvol
            elif tune_initLVvol<0:
                tune_initLVvol=1
                adj_try_initLVvol=adj_try_initLVvol*0.33
                try_initLVvol=try_initLVvol+tune_initLVvol*adj_try_initLVvol
        last_currentinitLVvol=currentinitLVvol
    if abs(cir_results[0,-1]/cir_results[0,1]-1)>0.05:
        logger.warning("Circuit might not converged, last cycle: start EDvol="+repr(cir_results[0,1])+', end EDvol='+repr(cir_results[0,-1]))
    return (try_initLVvol,adj_try_initLVvol)
