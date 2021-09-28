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

def createLVcircuit(casename,paramDict,stepTime=10,skipVariableList=None,verbose=True):

    logger.info('*** createLVcircuit ***')
    if skipVariableList is None:
        skipVariableList=[]
    cur_dir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))) 
    savePath=os.path.dirname(os.path.abspath(casename))
    LVcirfile = cur_dir + "/LV.cir"
    cirfilename = casename + ".cir"
    
    logger.info('cur_dir ' +repr(cur_dir))
    cmd = "cp " + LVcirfile + " " + cirfilename
    os.system(cmd)
    for key in keyToFill:
        if key not in skipVariableList:
            addstr='{:10.6f}   '.format(0)
            if key in paramDict:
                if isinstance(paramDict[key],(float,int)):
                    if paramDict[key]!=0:
                        if key in ['lvlvvalvk','lvlvvalvr']:
                            multiple_value=paramDict['Aortic_stenosis']
                        else:
                            multiple_value=1.
                        suffixDict_ind=max(-5,min(4,int(np.floor(np.log10(abs(paramDict[key]*multiple_value))/3.))))
                        addstr='{:10.6f}'.format(paramDict[key]*multiple_value/(1000.**suffixDict_ind))+suffixDict[suffixDict_ind]
            cmd = "sed -i.bak s/'<<"+key+">>'/'" + addstr + "'/g " + cirfilename
            os.system(cmd)
    if "cycle" not in skipVariableList:
        cmd = "sed -i.bak s/'<<cycle>>'/'" + '{:10.6f}'.format(paramDict['BCL']) + "m'/g " + cirfilename
        os.system(cmd)
    for side in ['l','r']:
        for tempstr in ['aamp','apeaktime','awidth']:
            if side+tempstr not in skipVariableList:
                cmd = "sed -i.bak s/'<<"+side+tempstr+">>'/'" + '{:10.6f}'.format(paramDict[side+tempstr]) + "'/g " + cirfilename
                os.system(cmd)
        if side+"vregurger" not in skipVariableList:
            cmd = "sed -i.bak s/'<<"+side+"vregurger>>'/'" + str(paramDict[side+"vregurger"]) + "'/g " + cirfilename
            os.system(cmd)
        if side+"vregurgevalveratio" not in skipVariableList:
            cmd = "sed -i.bak s/'<<"+side+"vregurgevalveratio>>'/'" + str(paramDict[side+"vregurgevalveratio"]) + "'/g " + cirfilename
            os.system(cmd)
    if "switchvalve" not in skipVariableList:
        cmd = "sed -i.bak s/'<<switchvalve>>'/'" + str(paramDict["switchvalve"])+ "'/g " + cirfilename
        os.system(cmd)
    if "rvfunc" not in skipVariableList:
        cmd = "sed -i.bak s/'<<rvfunc>>'/'" + paramDict['rvfunc'] + "'/g " + cirfilename
        os.system(cmd)
    if "timetopeaktension" not in skipVariableList:
        cmd = "sed -i.bak s/'<<timetopeaktension>>'/'" + str(paramDict['t0'])+ "m'/g " + cirfilename
        os.system(cmd)
    if "rvfuncarg" not in skipVariableList:
        for n in range(int(len(paramDict['rvfuncarg'])/2)):
            if n>0:
                cmd = "sed -i.bak s/'*Irvu"+str(n+1)+"'/'" + "Irvu"+str(n+1) + "'/g " + cirfilename
                os.system(cmd)
            suffixDict_ind=max(-5,min(4,int(np.floor(np.log10(abs(paramDict[paramDict['rvfuncarg'][n*2]]))/3.))))
            rvfuncarg='{:10.6f}'.format(paramDict[paramDict['rvfuncarg'][n*2]]/(1000.**suffixDict_ind))+suffixDict[suffixDict_ind]
            
            cmd = "sed -i.bak s/'<<rvuamp"+str(n+1)+">>'/'" + rvfuncarg + "'/g " + cirfilename
            os.system(cmd)
        
            rvfuncarg='{:10.6f}'.format(paramDict[paramDict['rvfuncarg'][n*2+1]])
            cmd = "sed -i.bak s/'<<rvuphase"+str(n+1)+">>'/'" + rvfuncarg + "'/g " + cirfilename
            os.system(cmd)
    
    if "stepTime" not in skipVariableList:
        cmd = "sed -i.bak s/'<<stepTime>>'/'" + str(stepTime)+'u'+ "'/g " + cirfilename
        os.system(cmd)
    if skipVariableList is None:
        rvdatabase=np.loadtxt(casename+'_rvflowrate.txt')
        rvdata=rvdatabase.copy()
        for n in range(1,11):
            rvdata=np.concatenate((rvdata,rvdatabase+np.array([[paramDict['BCL']*n/1000.,0.]])),axis=0)
        lvufilecontrol=casename+'_lvflowratecontrol.txt'
        np.savetxt(lvufilecontrol,rvdata)
        ngspice_py.simLVcircuit(casename,paramDict['BCL']*10.,lvufilecontrol,lvinputvar='i')
        lvdata=np.loadtxt(casename+'_circuit.txt',skiprows=1)[:,:2]
        np.savetxt(casename+'_lvucontrol.txt',lvdata)

