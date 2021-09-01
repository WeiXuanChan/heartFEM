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
'''
########################################################################
_version='3.1.2'
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

def simLVcircuit(casename,stopTime,lvufile,lvinputvar='V',initLAvol=0,initRAvol=0,initLVvol=0,initRVvol=0,vla0=None,vra0=None,init_file=None,init_time=None,additionalVarDict=None,timetopeak_from_to=None,verbose=True):

    logger.info('*** simulateLVcircuit ***')

    cirfilename = casename + ".cir"
    cirtempfilename = casename + "_temp.cir"
    cirlogfilename = casename + "_cir.log"
    cmd = "cp " + cirfilename + " " + cirtempfilename
    os.system(cmd)
    
    cmd = "sed -i.bak s/'<<stopTime>>'/'" + '{:.2f}m'.format(stopTime) + "'/g " + cirtempfilename
    os.system(cmd)
    if lvinputvar=='i' or lvinputvar=='I':
        cmd = "sed -i.bak s/'<<lvinputvar>>'/'I'/g " + cirtempfilename
        os.system(cmd)
        cmd = "sed -i.bak s/'Vlvu'/'*Vlvu'/g " + cirtempfilename
        os.system(cmd)
        cmd = "sed -i.bak s/'*Elvu'/'Elvu'/g " + cirtempfilename
        os.system(cmd)
    else:
        cmd = "sed -i.bak s/'<<lvinputvar>>'/'V'/g " + cirtempfilename
        os.system(cmd)
    if lvinputvar in ['i','v','I','V']:
        if os.path.isfile(lvufile):
            case_dir,lvufilename = os.path.split(lvufile)
        else:
            data=np.array(lvufile)
            case_dir=os.getcwd()
            lvufilename='lvufile.txt'
            lvufile=os.getcwd()+'/'+lvufilename
            np.savetxt(lvufile,data)
            
        if init_file is not None and init_time is not None:
            data=np.loadtxt(lvufile)
            init_data=np.loadtxt(init_file)
            init_data=init_data[init_data[:,0]<init_time]
            data=np.concatenate((init_data,data+np.array([[init_time,0.]])),axis=0)
            np.savetxt(os.getcwd()+'/temp_'+lvufilename,data)
        else:
            cmd = "cp "+lvufile+" " +'./temp_'+lvufilename
            os.system(cmd)
        cmd = "sed -i.bak s/'<<lvufile>>'/'temp_" +lvufilename + "'/g " + cirtempfilename
        os.system(cmd)
        if vla0 is None:
            if lvinputvar in ['v','V']:
                vla0=np.loadtxt(os.getcwd()+'/temp_'+lvufilename)[0,1]*0.99
            else:
                vla0=0.
        if vra0 is None:
            if lvinputvar in ['v','V']:
                vra0=np.loadtxt(os.getcwd()+'/temp_'+lvufilename)[0,1]*1.01
            else:
                vra0=0.
    elif lvinputvar[:7]=='table2d':
        case_dir,lvufilename = os.path.split(lvufile)
        cmd = "cp "+lvufile+" " +'./temp_'+lvufilename
        os.system(cmd)
        if timetopeak_from_to is not None:
            if timetopeak_from_to[1]!=timetopeak_from_to[0]:
                with open(os.getcwd()+'/temp_'+lvufilename,'r') as f:
                    lines=f.readlines()
                ytime=np.loadtxt(lvufile,skiprows=8,max_rows=1).reshape(-1)
                setAbove=(ytime>=(timetopeak_from_to[0]/1000.))
                setBelow=(ytime<(timetopeak_from_to[0]/1000.))
                ytime[setAbove]=ytime[setAbove]+(timetopeak_from_to[1]-timetopeak_from_to[0])/1000.
                ytime[setBelow]=ytime[setBelow]*(timetopeak_from_to[1]/timetopeak_from_to[0])
                ytime=np.round(ytime,decimals=4)
                lines[8]=' '.join(ytime.astype(str))+'\n'
                with open(os.getcwd()+'/temp_'+lvufilename,'w') as f:
                    f.writelines(lines)
            if additionalVarDict is not None:
                if 'timetopeaktension' not in additionalVarDict:
                    additionalVarDict['timetopeaktension']=str(timetopeak_from_to[1])+'m'
            else:
                additionalVarDict={'timetopeaktension':str(timetopeak_from_to[1])+'m'}
        cmd = "cp "+lvufile[:-4]+'base.txt'+" " +'./temp_'+lvufilename[:-4]+'base.txt'
        os.system(cmd)
        cmd = "sed -i.bak s/'<<lvutablefile>>'/'temp_" +lvufilename + "'/g " + cirtempfilename
        os.system(cmd)
        cmd = "sed -i.bak s/'<<lvutablebasefile>>'/'temp_" +lvufilename[:-4]+'base.txt' + "'/g " + cirtempfilename
        os.system(cmd)
        cmd = "sed -i.bak s/'alvu1'/'*alvu1'/g " + cirtempfilename
        os.system(cmd)
        cmd = "sed -i.bak s/'*alvu20'/'alvu20'/g " + cirtempfilename
        os.system(cmd)
        cmd = "sed -i.bak s/'*Blvu2'/'Blvu2'/g " + cirtempfilename
        os.system(cmd)
        
        cmd = "sed -i.bak s/'.model pwl_input'/'*.model pwl_input'/g " + cirtempfilename
        os.system(cmd)
        cmd = "sed -i.bak s/'*.model lvutabmod'/'.model lvutabmod'/g " + cirtempfilename
        os.system(cmd)
        if vla0 is None or vra0 is None:
                data0=np.loadtxt(os.getcwd()+'/temp_'+lvufilename[:-4]+'base.txt',skiprows=10)
                x0=np.loadtxt(os.getcwd()+'/temp_'+lvufilename[:-4]+'base.txt',skiprows=8,max_rows=1)
                
                datafunc = interpolate.splrep(x0,data0[:,0])
                temp_vlv = interpolate.splev(np.array([initLVvol]), datafunc)[0]
                logger.info('LV pressure start is '+repr(temp_vlv))
                if vla0 is None:
                    vla0 = temp_vlv*0.95
                if vra0 is None:
                    vra0 = temp_vlv*1.05
        if lvinputvar=='table2dtrackrelaxphase':
            cmd = "sed -i.bak s/'*alvu22'/'alvu22'/g " + cirtempfilename
            os.system(cmd)
            cmd = "sed -i.bak s/'PARAM trackrelaxphase = 0'/'PARAM trackrelaxphase = 1'/g " + cirtempfilename
            os.system(cmd)
            cmd = "cp "+lvufile[:-4]+'tr.txt'+" " +'./temp_'+lvufilename[:-4]+'tr.txt'
            os.system(cmd)
            cmd = "sed -i.bak s/'<<lvtrtablefile>>'/'temp_"+lvufilename[:-4]+"tr.txt'/g " + cirtempfilename
            os.system(cmd)
        else:
            cmd = "sed -i.bak s/'*alvu21'/'alvu21'/g " + cirtempfilename
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
            cmd = "sed -i.bak s/'<<"+key+">>'/'"+additionalVarDict[key]+"'/g " + cirtempfilename
            os.system(cmd)
    cmd = "ngspice -o "+cirlogfilename+" -b " + cirtempfilename
    os.system(cmd)
    
    cmd = "mv " +'./'+outfilename+'circuit.txt '+case_dir+'/'+'circuit_results.txt'
    os.system(cmd)
    cmd = "rm "+'./temp_'+lvufilename
    os.system(cmd)
    if lvinputvar[:7]=='table2d':
        cmd = "rm "+'./temp_'+lvufilename[:-4]+'base.txt'
        os.system(cmd)
        if lvinputvar=='table2dtrackrelaxphase':
            cmd = "rm "+'./temp_'+lvufilename[:-4]+'tr.txt'
            os.system(cmd)
    #if lvinputvar=='i' or lvinputvar=='v':
    #    cmd = "rm " +cur_dir+'/temp_'+lvufilename
    #    os.system(cmd)

 
