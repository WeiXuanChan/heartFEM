########################################################################

import sys
import vtk
import os
import inspect
import ngspice_py
import numpy as np
from scipy import interpolate
########################################################################

suffixDict={4:'T  ',3:'g  ',2:'meg',1:'k  ',0:' ',-1:'m  ',-2:'u  ',-3:'m  ',-4:'p  ',-5:'f  '}

def simLVcircuit(casename,stopTime,lvufile,lvinputvar='V',initLAvol=0,initRAvol=0,initLVvol=0,initRVvol=0,vla0=None,vra0=None,init_file=None,init_time=None,verbose=True):

    if (verbose): print ('*** createLVcircuit ***')

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
                print(x0)
                print(data0[:,0])
                datafunc = interpolate.splrep(x0,data0[:,0])
                temp_vlv = interpolate.splev(np.array([initLVvol]), datafunc)[0]
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
    
    cmd = "ngspice -o "+cirlogfilename+" -b " + cirtempfilename
    os.system(cmd)
    
    cmd = "mv " +'./'+outfilename+'circuit.txt '+case_dir+'/'+outfilename+'_circuit.txt'
    os.system(cmd)

    #if lvinputvar=='i' or lvinputvar=='v':
    #    cmd = "rm " +cur_dir+'/temp_'+lvufilename
    #    os.system(cmd)

 
