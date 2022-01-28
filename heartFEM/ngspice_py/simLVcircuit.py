'''
File: simLVcircuit.py
Description: simulate circuit with ngspice
History:
    Date    Programmer SAR# - Description
    ---------- ---------- ----------------------------
  Author: w.x.chan@gmail.com         08MAR2021           - Created
  Author: w.x.chan@gmail.com         08MAR2021           - v1.0.0
  Author: w.x.chan@gmail.com         13APR2021           - v2.0.0
  Author: w.x.chan@gmail.com         21APR2021           - v2.1.0
  Author: w.x.chan@gmail.com         09JUN2021           - v3.0.0
                                                              -added remove "temp_"+lvufilename[:-4] files
  Author: w.x.chan@gmail.com         10JUL2021           - v3.1.2
                                                              - debug lvflowrate input run by adding voltage lvgnd2 to gnd as v(lv)
  Author: w.x.chan@gmail.com         09Nov2021           - v4.0.0
                                                              -change time unit to ms
                                                              -added flexible cavity source
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


def simLVcircuit(casename,stopTime,sourcefileDict,initLAvol=0,initRAvol=0,initLVvol=0,initRVvol=0,vla0=None,vra0=None,init_file=None,init_time=None,additionalVarDict=None,timetopeak_from_to=None,verbose=True):

    logger.info('*** simulateLVcircuit ***')

    cirfilename = casename + ".cir"
    cirtempfilename = casename + "_temp.cir"
    cirlogfilename = casename + "_cir.log"
    cmd = "cp " + cirfilename + " " + cirtempfilename
    os.system(cmd)
    
    cmd = "sed -i.bak s/'<<stopTime>>'/'" + '{:.2f}'.format(stopTime) + "'/g " + cirtempfilename
    os.system(cmd)
    lvufilename={'la':'laufile.txt','ra':'raufile.txt','lv':'lvufile.txt','rv':'rvufile.txt'}
    for side in ['la','ra','lv','rv']:
        if 'time based' in sourcefileDict['Windkessel '+side.upper()+' source function']:
            if os.path.isfile(sourcefileDict['Windkessel '+side.upper()+' source file']):
                case_dir,lvufilename[side] = os.path.split(sourcefileDict['Windkessel '+side.upper()+' source file'])
            else:
                data=np.array(sourcefileDict['Windkessel '+side.upper()+' source file'])
                case_dir=os.getcwd()
                lvufilename[side]=side+'ufile.txt'
                lvufile=os.getcwd()+'/'+lvufilename[side]
                np.savetxt(lvufile,data)
                
            if init_file is not None and init_time is not None:
                data=np.loadtxt(sourcefileDict['Windkessel '+side.upper()+' source file'])
                init_data=np.loadtxt(init_file)
                init_data=init_data[init_data[:,0]<init_time]
                data=np.concatenate((init_data,data+np.array([[init_time,0.]])),axis=0)
                np.savetxt(os.getcwd()+'/temp_'+lvufilename[side],data)
            else:
                cmd = "cp "+sourcefileDict['Windkessel '+side.upper()+' source file']+" " +'./temp_'+lvufilename[side]
                os.system(cmd)
            cmd = "sed -i.bak s/'<<lvufile>>'/'temp_" +lvufilename[side] + "'/g " + cirtempfilename
            os.system(cmd)
            if side=='lv':
                if vla0 is None:
                    if 'pressure' in sourcefileDict['Windkessel '+side.upper()+' source function']:
                        vla0=np.loadtxt(os.getcwd()+'/temp_'+lvufilename[side])[0,1]*0.99
                    else:
                        vla0=0.
                if vra0 is None:
                    if 'pressure' in sourcefileDict['Windkessel '+side.upper()+' source function']:
                        vra0=np.loadtxt(os.getcwd()+'/temp_'+lvufilename[side])[0,1]*1.01
                    else:
                        vra0=0.
        elif '2D table' in sourcefileDict['Windkessel '+side.upper()+' source function']:
            case_dir,lvufilename[side] = os.path.split(sourcefileDict['Windkessel '+side.upper()+' source file'])
            cmd = "cp "+sourcefileDict['Windkessel '+side.upper()+' source file']+" " +'./temp_'+lvufilename[side]
            os.system(cmd)
            if timetopeak_from_to is not None:
                if timetopeak_from_to[1]!=timetopeak_from_to[0]:
                    with open(os.getcwd()+'/temp_'+lvufilename[side],'r') as f:
                        lines=f.readlines()
                    ytime=np.loadtxt(sourcefileDict['Windkessel '+side.upper()+' source file'],skiprows=8,max_rows=1).reshape(-1)
                    setAbove=(ytime>=(timetopeak_from_to[0]))
                    setBelow=(ytime<(timetopeak_from_to[0]))
                    ytime[setAbove]=ytime[setAbove]+(timetopeak_from_to[1]-timetopeak_from_to[0])
                    ytime[setBelow]=ytime[setBelow]*(timetopeak_from_to[1]/timetopeak_from_to[0])
                    ytime=np.round(ytime,decimals=4)
                    lines[8]=' '.join(ytime.astype(str))+'\n'
                    with open(os.getcwd()+'/temp_'+lvufilename[side],'w') as f:
                        f.writelines(lines)
                if additionalVarDict is not None:
                    if 'timetopeaktension' not in additionalVarDict:
                        additionalVarDict['timetopeaktension']=str(timetopeak_from_to[1])
                else:
                    additionalVarDict={'timetopeaktension':str(timetopeak_from_to[1])}
            cmd = "cp "+sourcefileDict['Windkessel '+side.upper()+' source file'][:-4]+'base.txt'+" " +'./temp_'+lvufilename[side][:-4]+'base.txt'
            os.system(cmd)
            cmd = "sed -i.bak s/'<<lvutablefile>>'/'temp_" +lvufilename[side] + "'/g " + cirtempfilename
            os.system(cmd)
            cmd = "sed -i.bak s/'<<lvutablebasefile>>'/'temp_" +lvufilename[side][:-4]+'base.txt' + "'/g " + cirtempfilename
            os.system(cmd)
            if side=='lv':
                if vla0 is None or vra0 is None:
                        data0=np.loadtxt(os.getcwd()+'/temp_'+lvufilename[side][:-4]+'base.txt',skiprows=10)
                        x0=np.loadtxt(os.getcwd()+'/temp_'+lvufilename[side][:-4]+'base.txt',skiprows=8,max_rows=1)
                        
                        datafunc = interpolate.splrep(x0,data0[:,0])
                        temp_vlv = interpolate.splev(np.array([initLVvol]), datafunc)[0]
                        logger.info('LV pressure start is '+repr(temp_vlv))
                        if vla0 is None:
                            vla0 = temp_vlv*0.95
                        if vra0 is None:
                            vra0 = temp_vlv*1.05
            if 'tracking of relaxation phase' in sourcefileDict['Windkessel '+side.upper()+' source function']:
                cmd = "cp "+sourcefileDict['Windkessel '+side.upper()+' source file'][:-4]+'tr.txt'+" " +'./temp_'+lvufilename[side][:-4]+'tr.txt'
                os.system(cmd)

    cmd = "sed -i.bak s/'<<lainitvol>>'/'"+str(initLAvol)+"'/g " + cirtempfilename
    os.system(cmd)
    cmd = "sed -i.bak s/'<<rainitvol>>'/'"+str(initRAvol)+"'/g " + cirtempfilename
    os.system(cmd)
    cmd = "sed -i.bak s/'<<lvinitvol>>'/'"+str(initLVvol)+"'/g " + cirtempfilename
    os.system(cmd)
    cmd = "sed -i.bak s/'<<rvinitvol>>'/'"+str(initRVvol)+"'/g " + cirtempfilename
    os.system(cmd)
    cmd = "sed -i.bak s/'<<vla0>>'/'"+str(vla0)+"'/g " + cirtempfilename
    os.system(cmd)
    cmd = "sed -i.bak s/'<<vra0>>'/'"+str(vra0)+"'/g " + cirtempfilename
    os.system(cmd)
        
    case_dir,outfilename = os.path.split(casename)
    cmd = "sed -i.bak s/'<<outfile>>'/'" + outfilename+'circuit.txt' + "'/g " + cirtempfilename
    os.system(cmd)
    
    if additionalVarDict is not None:
        for key in additionalVarDict:
            cmd = "sed -i.bak s/'<<"+key+">>'/'"+str(additionalVarDict[key])+"'/g " + cirtempfilename
            os.system(cmd)
    cmd = "ngspice -o "+cirlogfilename+" -b " + cirtempfilename
    os.system(cmd)
    
    cmd = "mv " +'./'+outfilename+'circuit.txt '+case_dir+'/'+'circuit_results.txt'
    os.system(cmd)
    for side in ['la','ra','lv','rv']:
        cmd = "rm "+'./temp_'+lvufilename[side]
        os.system(cmd)
        if '2D table' in sourcefileDict['Windkessel '+side.upper()+' source function']:
            cmd = "rm "+'./temp_'+lvufilename[side][:-4]+'base.txt'
            os.system(cmd)
            if 'tracking of relaxation phase' in sourcefileDict['Windkessel '+side.upper()+' source function']:
                cmd = "rm "+'./temp_'+lvufilename[side][:-4]+'tr.txt'
                os.system(cmd)

def simcirFile(cirtempfilename,sourcefileDict,BCL=None,toFolder=None,init_file=None,init_time=None,timetopeak_from_to=None,additionalVarDict=None,editCircuit=None,identifyName=None):
    if identifyName is None:
        identifyName='newcircuit'
    if editCircuit is None:
        editCircuit={}
    newcirfilename=cirtempfilename[:-4]+'_'+identifyName+'.cir'
    cirlogfilename=newcirfilename[:-4]+'.log'
    cmd = "cp " + cirtempfilename + " " + newcirfilename
    os.system(cmd)
    casename,temp = os.path.split(cirtempfilename)
    with open(newcirfilename,'r') as f:
        lines=f.readlines()
    for n in range(len(lines)):
        if lines[n][:-1] in editCircuit.keys():
            lines[n]=editCircuit[lines[n][:-1]]+'\n'
        if lines[n][:6]=='wrdata':
            outfilename=lines[n].split(" ")[1]
    with open(newcirfilename,'w') as f:
        f.writelines(lines)
    
    lvufilename={'la':'laufile.txt','ra':'raufile.txt','lv':'lvufile.txt','rv':'rvufile.txt'}
    for side in ['la','ra','lv','rv']:
        if 'time based' in sourcefileDict['Windkessel '+side.upper()+' source function']:
            if os.path.isfile(sourcefileDict['Windkessel '+side.upper()+' source file']):
                case_dir,lvufilename[side] = os.path.split(sourcefileDict['Windkessel '+side.upper()+' source file'])
            else:
                data=np.array(sourcefileDict['Windkessel '+side.upper()+' source file'])
                case_dir=os.getcwd()
                lvufilename[side]=side+'ufile.txt'
                lvufile=os.getcwd()+'/'+lvufilename[side]
                np.savetxt(lvufile,data)
                
            if init_file is not None and init_time is not None:
                data=np.loadtxt(sourcefileDict['Windkessel '+side.upper()+' source file'])
                init_data=np.loadtxt(init_file)
                init_data=init_data[init_data[:,0]<init_time]
                data=np.concatenate((init_data,data+np.array([[init_time,0.]])),axis=0)
                np.savetxt(os.getcwd()+'/temp_'+lvufilename[side],data)
            else:
                cmd = "cp "+sourcefileDict['Windkessel '+side.upper()+' source file']+" " +'./temp_'+lvufilename[side]
                os.system(cmd)
            cmd = "sed -i.bak s/'<<lvufile>>'/'temp_" +lvufilename[side] + "'/g " + cirtempfilename
            os.system(cmd)
        elif '2D table' in sourcefileDict['Windkessel '+side.upper()+' source function']:
            case_dir,lvufilename[side] = os.path.split(sourcefileDict['Windkessel '+side.upper()+' source file'])
            cmd = "cp "+sourcefileDict['Windkessel '+side.upper()+' source file']+" " +'./temp_'+lvufilename[side]
            os.system(cmd)
            if timetopeak_from_to is not None:
                if timetopeak_from_to[1]!=timetopeak_from_to[0]:
                    with open(os.getcwd()+'/temp_'+lvufilename[side],'r') as f:
                        lines=f.readlines()
                    ytime=np.loadtxt(sourcefileDict['Windkessel '+side.upper()+' source file'],skiprows=8,max_rows=1).reshape(-1)
                    setAbove=(ytime>=(timetopeak_from_to[0]))
                    setBelow=(ytime<(timetopeak_from_to[0]))
                    ytime[setAbove]=ytime[setAbove]+(timetopeak_from_to[1]-timetopeak_from_to[0])
                    ytime[setBelow]=ytime[setBelow]*(timetopeak_from_to[1]/timetopeak_from_to[0])
                    ytime=np.round(ytime,decimals=4)
                    lines[8]=' '.join(ytime.astype(str))+'\n'
                    with open(os.getcwd()+'/temp_'+lvufilename[side],'w') as f:
                        f.writelines(lines)
            cmd = "cp "+sourcefileDict['Windkessel '+side.upper()+' source file'][:-4]+'base.txt'+" " +'./temp_'+lvufilename[side][:-4]+'base.txt'
            os.system(cmd)
            if 'tracking of relaxation phase' in sourcefileDict['Windkessel '+side.upper()+' source function']:
                cmd = "cp "+sourcefileDict['Windkessel '+side.upper()+' source file'][:-4]+'tr.txt'+" " +'./temp_'+lvufilename[side][:-4]+'tr.txt'
                os.system(cmd)
    
    if additionalVarDict is not None:
        for key in additionalVarDict:
            cmd = "sed -i.bak s/'<<"+key+">>'/'"+str(additionalVarDict[key])+"'/g " + newcirfilename
            os.system(cmd)
    cmd = "ngspice -o "+cirlogfilename+" -b " + newcirfilename
    os.system(cmd)
    folder,savecirresultfilename= os.path.split(newcirfilename)
    if toFolder is not None:
        folder=folder+'/'+toFolder
        os.makedirs(folder,exist_ok=True)
    savecirresultfile=folder+'/'+savecirresultfilename[:-4]+'_circuit_results.txt'
        
    cmd = "mv " +'./'+outfilename+' '+savecirresultfile
    os.system(cmd)
    if BCL is not None:
        ngspice_py.getLastcycleCircuitResults(BCL,folder,savecirresultfilename[:-4]+'_circuit_results')
    for side in ['la','ra','lv','rv']:
        cmd = "rm "+'./temp_'+lvufilename[side]
        os.system(cmd)
        if '2D table' in sourcefileDict['Windkessel '+side.upper()+' source function']:
            cmd = "rm "+'./temp_'+lvufilename[side][:-4]+'base.txt'
            os.system(cmd)
            if 'tracking of relaxation phase' in sourcefileDict['Windkessel '+side.upper()+' source function']:
                cmd = "rm "+'./temp_'+lvufilename[side][:-4]+'tr.txt'
                os.system(cmd)
        
