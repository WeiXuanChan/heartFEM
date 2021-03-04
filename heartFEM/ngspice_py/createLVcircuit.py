########################################################################

import sys
import vtk
import os
import inspect
import ngspice_py
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
    
suffixDict={4:'T  ',3:'g  ',2:'meg',1:'k  ',0:' ',-1:'m  ',-2:'u  ',-3:'m  ',-4:'p  ',-5:'f  '}

def createLVcircuit(casename,paramDict,verbose=True):

    if (verbose): print ('*** createLVcircuit ***')

    cur_dir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))) 
    savePath=os.path.dirname(os.path.abspath(casename))
    LVcirfile = cur_dir + "/LV.cir"
    cirfilename = casename + ".cir"
    
    print('cur_dir',cur_dir)
    cmd = "cp " + LVcirfile + " " + cirfilename
    os.system(cmd)
    for key in keyToFill:
        addstr='{:10.6f}   '.format(0)
        if key in paramDict:
            if isinstance(paramDict[key],(float,int)):
                if paramDict[key]!=0:
                    suffixDict_ind=max(-5,min(4,int(np.floor(np.log10(abs(paramDict[key]))/3.))))
                    addstr='{:10.6f}'.format(paramDict[key]/(1000.**suffixDict_ind))+suffixDict[suffixDict_ind]
        cmd = "sed -i.bak s/'<<"+key+">>'/'" + addstr + "'/g " + cirfilename
        os.system(cmd)
    cmd = "sed -i.bak s/'<<cycle>>'/'" + '{:10.6f}'.format(paramDict['BCL']) + "m'/g " + cirfilename
    os.system(cmd)
    for side in ['l','r']:
        for tempstr in ['aamp','apeaktime','awidth']:
            cmd = "sed -i.bak s/'<<"+side+tempstr+">>'/'" + '{:10.6f}'.format(paramDict[side+tempstr]) + "'/g " + cirfilename
            os.system(cmd)
    cmd = "sed -i.bak s/'<<rvfunc>>'/'" + paramDict['rvfunc'] + "'/g " + cirfilename
    os.system(cmd)
    
    cmd = "sed -i.bak s/'<<timetopeaktension>>'/'" + str(paramDict['t0'])+ "m'/g " + cirfilename
    os.system(cmd)
    
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
    

    cmd = "sed -i.bak s/'<<stepTime>>'/'" + '10u'+ "'/g " + cirfilename
    os.system(cmd)
    
    rvdatabase=np.loadtxt(casename+'_rvflowrate.txt')
    rvdata=rvdatabase.copy()
    for n in range(1,11):
        rvdata=np.concatenate((rvdata,rvdatabase+np.array([[paramDict['BCL']*n/1000.,0.]])),axis=0)
    lvufilecontrol=casename+'_lvflowratecontrol.txt'
    np.savetxt(lvufilecontrol,rvdata)
    ngspice_py.simLVcircuit(casename,paramDict['BCL']*10.,lvufilecontrol,lvinputvar='i')
    lvdata=np.loadtxt(casename+'_circuit.txt',skiprows=1)[:,:2]
    np.savetxt(casename+'_lvucontrol.txt',lvdata)
