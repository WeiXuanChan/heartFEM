'''
File: ngspice_unit.py
Description: load all class for ngspice_py
             Contains externally usable class
History:
    Date    Programmer SAR# - Description
    ---------- ---------- ----------------------------
  Author: w.x.chan@gmail.com         15DEC2022           - Created
  Author: w.x.chan@gmail.com         15DEC2022           - v4.0.0
'''
_version='4.0.0'

import numpy as np

Components=['lv','la','rv','ra','aa','ao1','ao2','ao3','ao4','br','ca','ub','he','inte','ivc','kid','leg','lung','pa1','pa2','plac','svc','uv','art','ven']
linkComponents=['aaao1','ao1ao2','ao2ao3','ao3ao4','pa1pa2','pa2lung','da','ao1ca','cabr','brsvc','ao1ub','ubsvc','ao3he','ao3inte','intehe','ao3kid','kidivc','ao4plac','placuv','ao4leg','legivc','uvhe','heivc','dv','svcra','ivcra','lungla','fo','raravalv','rvrvvalv','lvlvvalv','lalavalv','lvart','artven','venla']
keyToFill=[]
for comp in Components:
    keyToFill.append(comp+'c')

for comp in linkComponents:
    keyToFill.append(comp+'r')
    keyToFill.append(comp+'l')
    keyToFill.append(comp+'k')
    keyToFill.append(comp+'b')
    
suffixDict={4:'T  ',3:'g  ',2:'meg',1:'k  ',0:' ',-1:'m  ',-2:'u  ',-3:'n  ',-4:'p  ',-5:'f  '}

def getLastcycleCircuitResults(BCL,path,filename):
    #BCL in ms
    cir_results=np.loadtxt(path+'/'+filename+'.txt',skiprows=1)
    with open(path+'/'+filename+'.txt','r') as f:
        cir_results_header=f.readline()
    if cir_results_header[-1]=='\n':
        cir_results_header=cir_results_header[:-1]
    cir_results_header=cir_results_header[1:]
    while cir_results[-1,0]>=(2.*BCL):
        cir_results[:,0]-=BCL
    cir_results=cir_results[cir_results[:,0]>=0]
    cir_results=cir_results[cir_results[:,0]<BCL]
    np.savetxt(path+'/'+filename+'_lastcycle.txt',cir_results,header=cir_results_header,comments='#')
    return