case='/media/yaplab/DATA/fetalheartFEM/preSurgery/Case_088_all'

import logging
logging.basicConfig(level=logging.INFO,filename=case+'/optimiseHeart.txt')
import heartFEM
import numpy as np
import medImgProc.processFunc as pf
import os
from scipy import optimize


age=np.loadtxt(case+'/age.txt')
hf=heartFEM.LVclosed(defaultAge='fetal'+str(age))
normal_Avalve_k=hf.defaultParameters['lvlvvalvk']
hf.casename=case
hf.meshname='t0'
volumeFile='volume2.txt'
rvflowrateFile='t0_rvflowrate.txt'

hf.defaultParameters['BCL']=np.loadtxt(case+'/heartcycle.txt')
Laxis=np.loadtxt(case+'/Laxis.txt')
volumes=np.loadtxt(case+'/'+volumeFile)[:,1]
hf.defaultParameters['EDV_LV']=volumes.max()
hf.defaultParameters['ESV_LV']=volumes.min()

hf.defaultParameters['EDP_LV'] = 5.
hf.defaultParameters['endo_angle']=60.
hf.defaultParameters['epi_angle']=-60.
hf.defaultParameters['clip_ratio']=0.98
hf.defaultParameters['t0'] = 140.5



hf.defaultRunMode='LVbehaviorRun'
hf.defaultParameters['Laxis_X']=Laxis[0]
hf.defaultParameters['Laxis_Y']=Laxis[1]
hf.defaultParameters['Laxis_Z']=Laxis[2]

hf.setAtrialpulse(0,1,1)
hf.setRVcurrentInflowFourier(case+'/'+rvflowrateFile,fourierTerms=4,fineAdjustPhase=True)
hf.setManualVolumeFile(case+'/'+volumeFile)

pressureDiff_LVAO_LVLA=np.loadtxt(case+'/doppler_LVAO_LVLA_pressure.txt')

default_run={'EDP_LV':hf.defaultParameters['EDP_LV'],'aaao1r':hf.defaultParameters['aaao1r'],'T0_LV':hf.defaultParameters['T0_LV']}

#heartFEM.optimiseFiberAngle(hf,tryPressureRatio=3.188525895/5.)
#results=hf.readRunParameters(case+'/optimiseFiberAngle/best_fit')
#hf.defaultParameters['endo_angle']=results['endo_angle']
#hf.defaultParameters['epi_angle']=results['epi_angle']

#heartFEM.optimiseEDP(hf,pressureDiff_LVAO_LVLA[2,1],folder=None,try_stlPressure=5.)
#os.makedirs(case+'/LVBehavior',exist_ok=True)
#cmd = "cp -r " + case+'/optimiseEDP/best_fit'+'/* '+ case+'/LVBehavior'
#os.system(cmd)

hf.defaultParameters['EDP_LV']=hf.readRunParameters(case+'/LVBehavior')['EDP_LV']

#hf.runCount="LVBehavior"
#hf.LVbehaviorRun(unloadGeo=True,folderToLVbehavior=hf.runCount)

hf.defaultParameters['ESV_LV']=volumes.min()
pressureDiff_LVAO=pressureDiff_LVAO_LVLA[0,1]
heartFEM.optimiseAllWinkesselParameters(hf,hf.defaultParameters['EDV_LV']-hf.defaultParameters['ESV_LV'],pressureDiff_LVAO_LVLA[1,1],pressureDiff_LVAO,Q_regurge_flowratio=0.6,folderToLVbehavior='LVBehavior',runCycles=20)

for lastiter in range(0,10):
    if not(os.path.isfile(case+'/optimiseWinkesselParameters/'+str(lastiter+1)+'/init_WindScale/best_fit/circuit_results_lastcycle.txt')):
        break
results=hf.readRunParameters(case+'/optimiseWinkesselParameters/'+str(lastiter)+'/init_WindScale/best_fit')
hf.defaultParameters['lvlvvalvk']=results['lvlvvalvk']
hf.defaultParameters['lvregurgevalveratio']=results['lvregurgevalveratio']
hf.defaultParameters['Windkessel_scale_T0_LV']=results['Windkessel_scale_T0_LV']
diesease_Avalve_k=hf.defaultParameters['lvlvvalvk']
hf.defaultParameters['aaregurgevalveratio']=1
os.makedirs(case+'/simpostsurgery',exist_ok=True)
for n in np.linspace(0,1,num=10):
    hf.runCount="simpostsurgery/lvlvvalvk_{0:.2f}".format(n)
    hf.defaultParameters['lvlvvalvk']=diesease_Avalve_k*(1.-n)+n*normal_Avalve_k
    hf.fullWindkesselRun(unloadGeo=True,folderToLVbehavior="LVBehavior",runCycles=20)