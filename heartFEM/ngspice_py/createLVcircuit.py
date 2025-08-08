'''
File: createLVcircuit.py
Description: creates a ngspice circuit from templete
History:
    Date    Programmer SAR# - Description
    ---------- ---------- ----------------------------
  Author: w.x.chan@gmail.com         08MAR2021           - Created
  Author: w.x.chan@gmail.com         08MAR2021           - v1.0.0
  Author: w.x.chan@gmail.com         13APR2021           - v2.0.0
  Author: w.x.chan@gmail.com         21APR2021           - v2.1.0
  Author: w.x.chan@gmail.com         10Jun2021           - v3.0.0
                                                              -added lvregurger/k/b and rvregurger/k/b for regurgitation on lv and rv to la and ra respectively
  Author: w.x.chan@gmail.com         31Aug2021           - v3.5.0
                                                              -added Aortic_stenosis to multiple resistances from LV to AO
  Author: w.x.chan@gmail.com         03Sep2021           - v3.5.1
                                                              -added stepTime, in micro seconds
                                                              - added switchvalve
  Author: w.x.chan@gmail.com         29Sep2021           - v3.6.0
                                                              -added aortic and pulmonary valve regurgitation
  Author: w.x.chan@gmail.com         09Nov2021           - v4.0.0
                                                              -remove aortic stenosis parameter
                                                              -change time unit to ms
  Author: w.x.chan@gmail.com         16FEB2024           - v4.3.0
                                                              -added 'Windkessel heart cycle control'
'''
########################################################################
_version='4.3.0'
import logging
logger = logging.getLogger(__name__)

import sys
import vtk
import os
import inspect
from . import ngspice_util
import numpy as np
########################################################################

def createLVcircuit(casename,paramDict,stepTime=0.01,skipVariableList=None,verbose=True):

    logger.info('*** createLVcircuit ***')
    if 'reference ngspice circuit filename' not in paramDict:
        refcirfilename="LV"
    else:
        refcirfilename=paramDict["reference ngspice circuit filename"]
    if skipVariableList is None:
        skipVariableList=[]
    cur_dir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))) 
    savePath=os.path.dirname(os.path.abspath(casename))
    LVcirfile = cur_dir + "/"+refcirfilename+".cir"
    cirfilename = casename + ".cir"
    
    logger.info('cur_dir ' +repr(cur_dir))
    cmd = "cp " + LVcirfile + " " + cirfilename
    os.system(cmd)
    for key in ngspice_util.keyToFill:
        if key not in skipVariableList:
            addstr='{:10.6f}   '.format(0)
            if key in paramDict:
                if isinstance(paramDict[key],(float,int)):
                    if paramDict[key]!=0:
                        suffixDict_ind=max(-5,min(4,int(np.floor(np.log10(abs(paramDict[key]))/3.))))
                        addstr='{:10.6f}'.format(paramDict[key]/(1000.**suffixDict_ind))+ngspice_util.suffixDict[suffixDict_ind]
                else:
                    addstr=str(paramDict[key])
            cmd = "sed -i.bak s/'<<"+key+">>'/'" + addstr + "'/g " + cirfilename
            os.system(cmd)
    for key in paramDict:
        if key[:4]=="cir_":
            if isinstance(paramDict[key],(float,int)):
                if paramDict[key]!=0:
                    suffixDict_ind=max(-5,min(4,int(np.floor(np.log10(abs(paramDict[key]))/3.))))
                    addstr='{:10.6f}'.format(paramDict[key]/(1000.**suffixDict_ind))+ngspice_util.suffixDict[suffixDict_ind]
            else:
                addstr=str(paramDict[key])
            cmd = "sed -i.bak s/'<<"+key+">>'/'" + addstr + "'/g " + cirfilename
            os.system(cmd)
    if "cycle" not in skipVariableList:
        cycle=paramDict['duration of one cardiac cycle in ms']
        cmd = "sed -i.bak s/'<<originalcycle>>'/'" + str(cycle) + "'/g " + cirfilename
        os.system(cmd)
        if paramDict['Windkessel heart cycle control'] is not None:
            cycle=paramDict['Windkessel heart cycle control']
        cmd = "sed -i.bak s/'<<cycle>>'/'" + str(cycle) + "'/g " + cirfilename
        os.system(cmd)
    lvufilename={'la':'laufile.txt','ra':'raufile.txt','lv':'lvufile.txt','rv':'rvufile.txt'}
    for side in ['la','ra','lv','rv']:
        if paramDict['Windkessel '+side.upper()+' source function']=='pulse':
            cmd = "sed -i.bak s/'<<"+side+"sourcemode>>'/'0'/g " + cirfilename
            os.system(cmd)
        elif paramDict['Windkessel '+side.upper()+' source function']=='fourier current':
            
            cmd = "sed -i.bak s/'<<"+side+"sourcemode>>'/'1'/g " + cirfilename
            os.system(cmd)
        if 'current' in paramDict['Windkessel '+side.upper()+' source function']:
            cmd = "sed -i.bak s/'<<"+side+"inputvar>>'/'I'/g " + cirfilename
            os.system(cmd)
        else:
            cmd = "sed -i.bak s/'<<"+side+"inputvar>>'/'V'/g " + cirfilename
            os.system(cmd)
        for tempstr,temp_paramStr in zip(['amp','peaktime','width'],['Windkessel {0:s} source pulse function peak pressure in mmHg','Windkessel {0:s} source pulse function time at peak pressure in ms','Windkessel {0:s} source pulse function pressure pulse width in ms']):
            if side+tempstr not in skipVariableList:
                cmd = "sed -i.bak s/'<<"+side+tempstr+">>'/'" + '{:10.6f}'.format(paramDict[temp_paramStr.format(side.upper())]) + "'/g " + cirfilename
                os.system(cmd)
        for n in range(4):
            if paramDict['Windkessel '+side.upper()+' source fourier function sine amplitude term '+str(n)]!=0:
                suffixDict_ind=max(-5,min(4,int(np.floor(np.log10(abs(paramDict['Windkessel '+side.upper()+' source fourier function sine amplitude term '+str(n)]))/3.))))
            else:
                suffixDict_ind=0
            rvfuncarg='{:10.6f}'.format(paramDict['Windkessel '+side.upper()+' source fourier function sine amplitude term '+str(n)]/(10.**(3*suffixDict_ind)))+ngspice_util.suffixDict[suffixDict_ind]
            
            cmd = "sed -i.bak s/'<<"+side+"uamp"+str(n+1)+">>'/'" + rvfuncarg + "'/g " + cirfilename
            os.system(cmd)
        
            rvfuncarg='{:10.6f}'.format(paramDict['Windkessel '+side.upper()+' source fourier function sine degree phase term '+str(n)])
            cmd = "sed -i.bak s/'<<"+side+"uphase"+str(n+1)+">>'/'" + rvfuncarg + "'/g " + cirfilename
            os.system(cmd)
        if 'time based' in paramDict['Windkessel '+side.upper()+' source function']:
            cmd = "sed -i.bak s/'<<"+side+"sourcemode>>'/'2'/g " + cirfilename
            os.system(cmd)
            cmd = "sed -i.bak s/'<<"+side+"ufile>>'/'temp_" +lvufilename[side]+ "'/g " + cirfilename
            os.system(cmd)
        elif '2D table' in paramDict['Windkessel '+side.upper()+' source function']:
            cmd = "sed -i.bak s/'<<"+side+"utablefile>>'/'temp_" +lvufilename[side] + "'/g " + cirfilename
            os.system(cmd)
            cmd = "sed -i.bak s/'<<"+side+"utablebasefile>>'/'temp_" +lvufilename[side][:-4]+'base.txt' + "'/g " + cirfilename
            os.system(cmd)
            if 'tracking of relaxation phase' in paramDict['Windkessel '+side.upper()+' source function']:
                cmd = "sed -i.bak s/'<<"+side+"sourcemode>>'/'4'/g " + cirfilename
                os.system(cmd)
                cmd = "sed -i.bak s/'<<"+side+"trtablefile>>'/'temp_"+lvufilename[side][:-4]+"tr.txt'/g " + cirfilename
                os.system(cmd)
            else:
                cmd = "sed -i.bak s/'<<"+side+"sourcemode>>'/'3'/g " + cirfilename
                os.system(cmd)
        elif '3D table' in paramDict['Windkessel '+side.upper()+' source function']:
            cmd = "sed -i.bak s/'<<"+side+"utablefile>>'/'temp_" +lvufilename[side] + "'/g " + cirfilename
            os.system(cmd)
            cmd = "sed -i.bak s/'<<"+side+"utablebasefile>>'/'temp_" +lvufilename[side][:-4]+'base.txt' + "'/g " + cirfilename
            os.system(cmd)
            if 'tracking of relaxation phase' in paramDict['Windkessel '+side.upper()+' source function']:
                cmd = "sed -i.bak s/'<<"+side+"sourcemode>>'/'6'/g " + cirfilename
                os.system(cmd)
                cmd = "sed -i.bak s/'<<"+side+"trtablefile>>'/'temp_"+lvufilename[side][:-4]+"tr.txt'/g " + cirfilename
                os.system(cmd)
            else:
                cmd = "sed -i.bak s/'<<"+side+"sourcemode>>'/'5'/g " + cirfilename
                os.system(cmd)
    for valve in ['lv','rv','aa','pa1']:
        if (valve+"regurgevalveratio" not in skipVariableList) and (valve+"regurgevalveratio" in paramDict):
            cmd = "sed -i.bak s/'<<"+valve+"regurgevalveratio>>'/'" + str(paramDict[valve+"regurgevalveratio"]) + "'/g " + cirfilename
            os.system(cmd)
    for side in ['la','ra','lv','rv']:
        if (side+"timetopeaktension" not in skipVariableList) and ('time to maximum '+side.upper()+' fiber tension in ms' in paramDict):
            peaktime=paramDict['time to maximum '+side.upper()+' fiber tension in ms']
            cmd = "sed -i.bak s/'<<"+side+"originaltimetopeaktension>>'/'" + str(peaktime)+ "'/g " + cirfilename
            os.system(cmd)
            if paramDict['Windkessel time to maximum '+side.upper()+' fiber tension control'] is not None:
                peaktime=paramDict['Windkessel time to maximum '+side.upper()+' fiber tension control']
            cmd = "sed -i.bak s/'<<"+side+"timetopeaktension>>'/'" + str(peaktime)+ "'/g " + cirfilename
            os.system(cmd)
    
    if "stepTime" not in skipVariableList:
        cmd = "sed -i.bak s/'<<stepTime>>'/'" + str(stepTime)+ "'/g " + cirfilename
        os.system(cmd)
    if skipVariableList is None:
        rvdatabase=np.loadtxt(casename+'_rvflowrate.txt')
        rvdata=rvdatabase.copy()
        for n in range(1,11):
            rvdata=np.concatenate((rvdata,rvdatabase+np.array([[paramDict['duration of one cardiac cycle in ms']*n,0.]])),axis=0)
        lvufilecontrol=casename+'_lvflowratecontrol.txt'
        np.savetxt(lvufilecontrol,rvdata)
        ngspice_util.simLVcircuit(casename,paramDict['duration of one cardiac cycle in ms']*10.,lvufilecontrol,lvinputvar='i')
        lvdata=np.loadtxt(casename+'_circuit.txt',skiprows=1)[:,:2]
        np.savetxt(casename+'_lvucontrol.txt',lvdata)

