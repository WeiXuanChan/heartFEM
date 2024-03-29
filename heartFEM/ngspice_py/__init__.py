'''
File: __init__.py
Description: load all class for ngspice_py
             Contains externally usable class
History:
    Date    Programmer SAR# - Description
    ---------- ---------- ----------------------------
  Author: w.x.chan@gmail.com         08MAR2021           - Created
  Author: w.x.chan@gmail.com         08MAR2021           - v1.0.0
                                                            -createLVcircuit v1.0.0
                                                            -simLVcircuit v1.0.0
                                                            -generateLVtable v1.0.0
  Author: w.x.chan@gmail.com         08MAR2021           - v2.0.0
                                                            -createLVcircuit v2.0.0
                                                            -simLVcircuit v2.0.0
                                                            -generateLVtable v2.0.0
                                                            -simLVcircuit_alignEStime v2.0.0
  Author: w.x.chan@gmail.com         21APR2021           - v2.1.0
                                                            -createLVcircuit v2.1.0
                                                            -simLVcircuit v2.1.0
                                                            -generateLVtable v2.1.0
                                                            -simLVcircuit_alignEStime v2.1.0
  Author: w.x.chan@gmail.com         12MAY2021           - v2.3.4
                                                            -createLVcircuit v2.1.0
                                                            -simLVcircuit v2.1.0
                                                            -generateLVtable v2.3.4
                                                            -simLVcircuit_alignEStime v2.1.0
  Author: w.x.chan@gmail.com         25MAY2021           - v2.4.0
                                                            -createLVcircuit v2.1.0
                                                            -simLVcircuit v2.1.0
                                                            -generateLVtable v2.3.4
                                                            -simLVcircuit_alignEStime v2.1.0
                                                            -LV.cir updated to include output of phaseTime
  Author: w.x.chan@gmail.com         31MAY2021           - v2.5.0
                                                            -createLVcircuit v2.1.0
                                                            -simLVcircuit v2.1.0
                                                            -generateLVtable v2.3.4
                                                            -simLVcircuit_alignEStime v2.1.0
                                                            -LV.cir debuged for track phase mode
  Author: w.x.chan@gmail.com         10JUN2021           - v3.0.0
                                                            -createLVcircuit v3.0.0
                                                            -simLVcircuit v2.1.0
                                                            -generateLVtable v2.3.4
                                                            -simLVcircuit_alignEStime v2.1.0
                                                            -LV.cir v3.0.0
  Author: w.x.chan@gmail.com         10JUL2021           - v3.1.2
                                                            -createLVcircuit v3.0.0
                                                            -simLVcircuit v3.1.2
                                                            -generateLVtable v2.3.4
                                                            -simLVcircuit_alignEStime v2.1.0
                                                            -LV.cir v3.1.2
  Author: w.x.chan@gmail.com         31Aug2021           - v3.5.0
                                                            -createLVcircuit v3.5.0
                                                            -simLVcircuit v3.1.2
                                                            -generateLVtable v3.5.0
                                                            -simLVcircuit_alignEStime v2.1.0
                                                            -LV.cir v3.5.0
  Author: w.x.chan@gmail.com         09Nov2021           - v4.0.0
                                                            -createLVcircuit v4.0.0
                                                            -simLVcircuit v4.0.0
                                                            -generateLVtable v4.0.0
                                                            -simLVcircuit_alignEStime v4.0.0
                                                            -LV.cir v4.0.0
'''
_version='4.0.0'

from .createLVcircuit            import *
from .simLVcircuit               import *
from .generateLVtable            import *
from .simLVcircuit_alignEStime   import *

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