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
  Author: w.x.chan@gmail.com         28Sep2021           - v3.6.0
                                                              -debug remove exit when not converged
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
from heartFEM import ngspice_py
import numpy as np
from scipy import interpolate
########################################################################

def simLVcircuit_alignEStime(casename,stopTime,sourcefileDict,period,targetEStime,init_timetopeaktension,try_timetopeaktension=None,initLAvol=0,initRAvol=0,initLVvol=0,initRVvol=0,vla0=None,vra0=None,init_file=None,init_time=None,iterationNumber=100,verbose=True):

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
        case_dir,lvufilename = os.path.split(sourcefileDict['Windkessel LV source file'])
        timetopeak_from_to=[init_timetopeaktension,try_timetopeaktension]
        ngspice_py.simLVcircuit(casename,stopTime,sourcefileDict,initLAvol=initLAvol,initRAvol=initRAvol,initLVvol=initLVvol,initRVvol=initRVvol,vla0=vla0,vra0=vra0,init_file=init_file,init_time=init_time,timetopeak_from_to=timetopeak_from_to,verbose=verbose)
        
        cir_results=np.loadtxt(case_dir+'/'+'circuit_results.txt',skiprows=1)[:,2:4]
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

def simLVcircuit_align_EDvol_and_EStime(casename,stopTime,sourcefileDict,period,targetEStime,init_timetopeaktension,try_initLVvol=None,initLAvol=0,initRAvol=0,initLVvol=0,initRVvol=0,vla0=None,vra0=None,init_file=None,init_time=None,iterationNumber=100,verbose=True):
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
            ngspice_py.simLVcircuit(casename,stopTime,sourcefileDict,initLAvol=initLAvol,initRAvol=initRAvol,initLVvol=try_initLVvol,initRVvol=initRVvol,vla0=vla0,vra0=vra0,init_file=init_file,init_time=init_time,verbose=verbose)
        else:
            if n==0:
                try_timetopeaktension=None
            else:
                try_timetopeaktension=list(try_timetopeaktension)
                try_timetopeaktension[1]=try_timetopeaktension[1]*3.
            try_timetopeaktension=simLVcircuit_alignEStime(casename,stopTime,sourcefileDict,period,targetEStime,init_timetopeaktension,try_timetopeaktension=try_timetopeaktension,initLAvol=initLAvol,initRAvol=initRAvol,initLVvol=try_initLVvol,initRVvol=initRVvol,vla0=vla0,vra0=vra0,init_file=init_file,init_time=init_time,iterationNumber=iterationNumber,verbose=verbose)
        case_dir,outfilename = os.path.split(casename)
        cir_results=np.loadtxt(case_dir+'/'+'circuit_results.txt',skiprows=1)[:,2:4]
        removePeriodCount=0
        while cir_results[-1,0]>=(2.*period):
            cir_results[:,0]-=period
            removePeriodCount+=1
        cir_results=cir_results[cir_results[:,0]>=0]
        cir_results=cir_results[cir_results[:,0]<period]
        currentinitLVvol=cir_results[0,1]
        #np.savetxt(casename+'extracted_circuit_results_'+str(n)+'.txt',cir_results)
        logger.info("  Current EDvol="+repr(currentinitLVvol)+', target EDvol='+repr(targetinitLVvol)+" , time="+repr(cir_results[0,0]+removePeriodCount*period))
        if abs(currentinitLVvol-targetinitLVvol)<(targetinitLVvol*10.**-4.) or adj_try_initLVvol<(targetinitLVvol*10.**-5.):
            break
        elif currentinitLVvol>targetinitLVvol:
            if tune_initLVvol is None:
                tune_initLVvol=-1
                try_initLVvol=try_initLVvol*targetinitLVvol/currentinitLVvol
            else:
                if tune_initLVvol>0:
                    tune_initLVvol=-1
                    adj_try_initLVvol=adj_try_initLVvol*0.33
                try_initLVvol=try_initLVvol+tune_initLVvol*adj_try_initLVvol
        elif currentinitLVvol<targetinitLVvol:
            if tune_initLVvol is None:
                tune_initLVvol=1
                try_initLVvol=try_initLVvol*targetinitLVvol/currentinitLVvol
            else:
                if tune_initLVvol<0:
                    tune_initLVvol=1
                    adj_try_initLVvol=adj_try_initLVvol*0.33
                try_initLVvol=try_initLVvol+tune_initLVvol*adj_try_initLVvol
        last_currentinitLVvol=currentinitLVvol
    if abs(cir_results[0,-1]/cir_results[0,1]-1)>0.05:
        logger.warning("Circuit might not converged, last cycle: start EDvol="+repr(cir_results[0,1])+', end EDvol='+repr(cir_results[0,-1]))
    return (try_initLVvol,adj_try_initLVvol)
