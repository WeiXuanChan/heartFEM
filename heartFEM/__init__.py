'''
################################################
MIT License
Copyright (c) 2019 W. X. Chan
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
################################################
File: __init__.py
Description: load all class for heartFEM
             Contains externally usable class
History:
    Date    Programmer SAR# - Description
    ---------- ---------- ----------------------------
  Author: w.x.chan@gmail.com         08MAR2021           - Created
  Author: w.x.chan@gmail.com         08MAR2021           - v2.0.6
                                                            -ngspice_py v1.0.0
  Author: w.x.chan@gmail.com         08MAR2021           - v2.1.0
                                                            -ngspice_py v1.0.0
                                                            -add raise error when pressure in not a value
                                                            -add inverse calculate init mesh given diastolic mesh and pressure: function adjustMeshToPressure
  Author: w.x.chan@gmail.com         13APR2021           - v2.2.0
                                                            -ngspice_py v2.0.0
                                                            -tidy up to stable version
  Author: w.x.chan@gmail.com         21APR2021           - v2.3.0
                                                            -ngspice_py v2.1.0
                                                            -debug getLVbehaviour non-convergence
  Author: w.x.chan@gmail.com         21APR2021           - v2.3.1
                                                            -ngspice_py v2.1.0
                                                            -debug LVbehaviorRun for folderToLVbehavior=None where ngspice_py.generateLVtable use folderToLVbehavior
                                                            -debug solveTa to include changing heart.dt.dt
  Author: w.x.chan@gmail.com         21APR2021           - v2.3.3
                                                            -ngspice_py v2.1.0
                                                            -heartParameters v1.2.0
                                                            -debug iterativeRun to save mesh into folder and solve without incremental t_a
  Author: w.x.chan@gmail.com         12MAY2021           - v2.3.4
                                                            -ngspice_py v2.3.4
                                                            -heartParameters v1.2.0
  Author: w.x.chan@gmail.com         24MAY2021           - v2.4.0
                                                            -ngspice_py v2.4.0
                                                            -heartParameters v1.2.0
                                                            - added func getphaseFromFile and setManualPhaseTimeFile
                                                            - added manualPhaseTimeInput and setHeart to function: iterativeRun
                                                            - added 'self.manualVol=self.getVolumeFromFile' to setManualvolFile
  Author: w.x.chan@gmail.com         27MAY2021           - v2.4.1
                                                            -ngspice_py v2.4.0
                                                            -heartParameters v1.2.0
                                                            - debug iterativeRun with manualPhaseTimeInput=True
  Author: w.x.chan@gmail.com         31MAY2021           - v2.5.0
                                                            -ngspice_py v2.5.0
                                                            -heartParameters v1.2.0
                                                            - added inverse Heart to unloaded geometry with deformation gradient
  Author: w.x.chan@gmail.com         09JUN2021           - v3.0.0
                                                            -ngspice_py v3.0.0
                                                            -heartParameters v3.0.0
                                                            -lcleeHeart v3.0.0
                                                            -added readRunParameters
                                                            -added writeRunParameters
                                                            -added option to continue getLVbehavior if isinstance(runParameters,str)
                                                            -debug generateLVtable in LVbehaviorRun when folderToLVbehavior is None
                                                            - added save last cycle of circuit results in LVbehaviorRun
                                                            -read vla0 and vra0 for LVbehaviorRun
                                                            - remove case when ['ES_time'] not in runParameters in LVbehaviorRun
                                                            - added fullWindkesselRun in modes to run full windkessel only
  Author: w.x.chan@gmail.com         08JUL2021           - v3.1.0
                                                            -ngspice_py v3.0.0
                                                            -heartParameters v3.1.0
                                                            -lcleeHeart v3.1.0
                                                            -added outOfplaneDeg for meshing
  Author: w.x.chan@gmail.com         10JUL2021           - v3.1.2
                                                            -ngspice_py v3.1.2
                                                            -heartParameters v3.1.0
                                                            -lcleeHeart v3.1.0
                                                            -added outOfplaneDeg for meshing
  Author: w.x.chan@gmail.com         10JUL2021           - v3.3.0
                                                            -debugged for adjust fibers after ALE.move
                                                            - remove "self.manualVol is not None" condition to not do unloading
                                                            -ngspice_py v3.1.2
                                                            -heartParameters v3.3.0
                                                            -lcleeHeart v3.3.0
                                                            -removed outOfplaneDeg and added 'fiberSheetletAngle','fiberSheetletWidth','radialFiberAngle', set fiberlength for meshing as sarcomere length ['lr']
  Author: w.x.chan@gmail.com         28JUL2021           - v3.4.2
                                                            - added fenicsResultWriter
                                                            -ngspice_py v3.4.0
                                                            -heartParameters v3.3.0
                                                            -lcleeHeart v3.4.2
  Author: w.x.chan@gmail.com         12Aug2021           - v3.5.0
                                                            - added functions: optimiseFiberAngle
                                                            -debug cost_function to read _PV.txt if fail for other files
                                                            -debug optimiser_linker to +1 when runSimulation.runCount is str
                                                            -debug optimiser_linker to print status and cost in runParameters.txt
                                                            -debug setRVcurrentInflowFourier phase to be in degrees to match ngspice
                                                            -added set_kwargs for optimiser_linker
                                                            - debug consolidated stress output to be in main folder
                                                            - added costs: 'StrokeVolume' and 'MAXLVLApressure'
                                                            - added 'scale Windkessel active pressure' to sclae LVbehavior during table generation for ngspice
                                                            -added 'Aortic_stenosis' to multiple resistances from LV to AO
                                                            - change setRVcurrentFourier to setRVcurrentInflowFourier to prevent ambiguity
                                                            -ngspice_py v3.5.0
                                                            -heartParameters v3.5.0
                                                            -lcleeHeart v3.4.2
  Author: w.x.chan@gmail.com         01Sep2021           - v3.5.1
                                                            - added run one time of fullWinkessel at the start of step 2 in optimiseWinkesselParameters
                                                            - added function optimiseAllWinkesselParameters
                                                            -added option 'fineAdjustPhase' to LVclosed.setRVcurrentInflowFourier
                                                            -ngspice_py v3.5.0
                                                            -heartParameters v3.5.0
                                                            -lcleeHeart v3.4.2
  Author: w.x.chan@gmail.com         05Oct2021           - v3.6.0
                                                            -added set_EVP_with_LA_systolic_volume
                                                            -added optimiseStiffness
                                                            -added trackphase for iterativeRun
                                                            -ngspice_py v3.5.0
                                                            -heartParameters v3.6.0
                                                            -lcleeHeart v3.6.0
  Author: w.x.chan@gmail.com         05Oct2021           - v3.6.1
                                                            -ngspice_py v3.5.0
                                                            -heartParameters v3.6.0
                                                            -lcleeHeart v3.6.1
  Author: w.x.chan@gmail.com         28Oct2021           - v3.6.3
                                                            -in iterativeRun, save PV_.txt again to make sure of output  
                                                            -ngspice_py v3.5.0
                                                            -heartParameters v3.6.0
                                                            -lcleeHeart v3.6.1
  Author: w.x.chan@gmail.com         02Dec2021           - v3.6.4
                                                            -debug iterative run with unloadGeo=True , save and load heart at different folder  
                                                            -ngspice_py v3.5.0
                                                            -heartParameters v3.6.0
                                                            -lcleeHeart v3.6.1
  Author: w.x.chan@gmail.com         05Nov2021           - v4.0.0
                                                            -heartModelTranslator v4.0.0 
                                                                -new
                                                            -ngspice_py v4.0.0
                                                            -heartParameters v4.0.0
                                                            -lcleeHeart v4.0.0
  Author: w.x.chan@gmail.com         28Jan2022           - v4.0.1
                                                              -added notAlignEDV to LVbehaviorRun and fullWindkesselRun as an option not to align EDV
                                                            -heartModelTranslator v4.0.0 
                                                            -ngspice_py v4.0.0
                                                            -heartParameters v4.0.0
                                                            -lcleeHeart v4.0.0
  Author: w.x.chan@gmail.com         06Dec2022           - v4.1.0
                                                            -heartModelTranslator v4.1.0 
                                                            -ngspice_py v4.0.0
                                                            -heartParameters v4.1.0
                                                            -lcleeHeart v4.1.0
  Author: w.x.chan@gmail.com         07FEB2023           - v4.2.0
                                                            -heartModelTranslator v4.2.0 
                                                            -ngspice_py v4.0.0
                                                            -heartParameters v4.1.0
                                                            -lcleeHeart v4.1.0
  Author: w.x.chan@gmail.com         16FEB2023           - v4.3.0
                                                            -heartModelTranslator v4.3.0 
                                                            -ngspice_py v4.3.0
                                                            -heartParameters v4.3.0
                                                            -lcleeHeart v4.1.0
  Author: w.x.chan@gmail.com         16FEB2023           - v4.3.1
                                                              - add velstep to be indicative of a factor scale from minESV to max EDV
                                                            -heartModelTranslator v4.3.0 
                                                            -ngspice_py v4.3.0
                                                            -heartParameters v4.3.0
                                                            -lcleeHeart v4.1.0
  Author: w.x.chan@gmail.com         08NOV2023           - v4.3.2
                                                            -heartModelTranslator v4.3.0 
                                                            -ngspice_py v4.3.0
                                                            -heartParameters v4.3.0
                                                            -lcleeHeart v4.3.2
  Author: w.x.chan@gmail.com         07Jun2024           - v4.4.0
                                                            -heartModelTranslator v4.4.0 
                                                            -ngspice_py v4.3.0
                                                            -heartParameters v4.4.0
                                                            -lcleeHeart v4.4.0
  Author: w.x.chan@gmail.com         02Jan2025           - v4.5.0
                                                              - added unloaded volume to fix any chamber to be a specific unloaded volume
                                                            -heartModelTranslator v4.4.0 
                                                            -ngspice_py v4.3.0
                                                            -heartParameters v4.4.0
                                                            -lcleeHeart v4.4.0
  Author: w.x.chan@gmail.com         26Feb2025           - v4.5.1
                                                              - remove square brackets for printing PV_.txt
                                                              - enable continue_from_break with iterative run
                                                            -heartModelTranslator v4.4.0 
                                                            -ngspice_py v4.3.0
                                                            -heartParameters v4.4.0
                                                            -lcleeHeart v4.4.0
  Author: w.x.chan@gmail.com         08Aug2025           - v4.5.2
                                                              - added output results for LVbehaviour
                                                              - added results tag 
                                                            -heartModelTranslator v4.4.0 
                                                            -ngspice_py v4.3.0
                                                            -heartParameters v4.4.0
                                                            -lcleeHeart v4.4.0
'''
_version='4.5.1'
import logging
logger = logging.getLogger('heartFEM v'+_version)
logger.info('heartFEM version '+_version)

import sys

import os as os
import shutil
try:
    import dolfin as fenics
except:
    import fenics as fenics
try:
    from multiprocessing import Pool as Pool
    #from pathos.multiprocessing import ProcessingPool as Pool
    #from multiprocessing.pool import ThreadPool
    import subprocess
except:
    pass
import numpy as np
from scipy import interpolate as scipyinterpolate
from scipy.optimize import curve_fit
from scipy import optimize
from matplotlib import pylab as plt
from petsc4py import PETSc

import math
import csv
import re
import vtk
from . import ngspice_py
from . import heartParameters
from . import heartModelTranslator
#from mpi4py import MPI as pyMPI

from fenicstools.Probe import Probes
from fenicstools.Probe import array as Probes_array

import motionSegmentation.BsplineFourier as bsf
 
WindkesselComponents=['lv','la','rv','ra','aa','ao1','ao2','ao3','ao4','br','ca','ub','he','inte','ivc','kid','leg','lung','pa1','pa2','plac','svc','uv']
WindkessellinkComponents=['aaao1','ao1ao2','ao2ao3','ao3ao4','pa1pa2','pa2lung','da','ao1ca','cabr','brsvc','ao1ub','ubsvc','ao3he','ao3inte','intehe','ao3kid','kidivc','ao4plac','placuv','ao4leg','legivc','uvhe','heivc','dv','svcra','ivcra','lungla','fo','raravalv','rvrvvalv','lvlvvalv','rvrvvalv']
defaultAgeScalePower={'defaultr':-1.,'pa2lungr':-1.2,'lunglar':-1.2,'cabrr':-1.1,'brsvcr':-1.1,'dvr':-0.55,
                      'defaultl':-0.33,
                      'defaultc':1.33,'brc':1.471,'lungc':1.6,'ra':0.5,'la':0.5,
                      'defaultk':0.,'fok':-0.6,'dak': -2.5,'dvk':-0.88,'raravalvk':-1.33,'rvrvvalvk':-1.33,'lalavalvk':-1.33,'lvlvvalvk':-1.33,
                      'defaultb':0.}

std_multicpu_getLVbehavior_file=["import sys\n",
                                 "import os\n",
                                 "import logging\n",
                                 "logger = logging.getLogger('multicpu')\n",
                                 "dir_path = os.path.split(sys.argv[0])[0]\n",
                                 "if os.path.isfile(dir_path+'/multipcu_LVBehavior/'+str(int(sys.argv[1]))+'.npy'):\n",
                                 "    sys.exit()\n",
                                 "logging.basicConfig(filename=dir_path+'/multipcu_LVBehavior/'+sys.argv[1]+'.log', level=logging.DEBUG)\n"
                                 "if len(sys.argv)>2:\n",
                                 "    for n in range(len(sys.argv)-1,1,-1):\n",
                                 "        if not(sys.argv[n] in sys.path):\n",
                                 "            logger.info('added path '+sys.argv[n])\n",
                                 "            sys.path.insert(0,sys.argv[n])\n",
                                 "else:\n",
                                 "    logger.info('no added path '+repr(sys.argv))\n",
                                 "import heartFEM\n",
                                 "heartFEM.multicpu_getLVbehavior(dir_path,int(sys.argv[1]))\n"]
std_multicpu_iterativerun_file=["import sys\n",
                                 "import os\n",
                                 "import logging\n",
                                 "logger = logging.getLogger('multicpu')\n",
                                 "dir_path = os.path.split(sys.argv[0])[0]\n",
                                 "if os.path.isfile(dir_path+'/multipcu_iterative/'+str(int(sys.argv[1]))+'/PV.png'):\n",
                                 "    sys.exit()\n",
                                 "logging.basicConfig(filename=dir_path+'/multipcu_iterative/'+sys.argv[1]+'.log', level=logging.DEBUG)\n"
                                 "if len(sys.argv)>2:\n",
                                 "    for n in range(len(sys.argv)-1,1,-1):\n",
                                 "        if not(sys.argv[n] in sys.path):\n",
                                 "            logger.info('added path '+sys.argv[n])\n",
                                 "            sys.path.insert(0,sys.argv[n])\n",
                                 "else:\n",
                                 "    logger.info('no added path '+repr(sys.argv))\n",
                                 "import heartFEM\n",
                                 "heartFEM.multicpu_iterativerun(dir_path,int(sys.argv[1]))\n"]
std_pvd_line='    <DataSet timestep="{1:d}" part="0" file="{0:s}{1:06d}.vtu" />\n'
def obb(pts):
    ca = np.cov(pts,y = None,rowvar = 0,bias = 1)
    v, vect = np.linalg.eig(ca)
    tvect = np.transpose(vect)
    pts_rotated = np.dot(pts,vect)
    min_axis = np.min(pts_rotated,axis=0)
    max_axis = np.max(pts_rotated,axis=0)
    boxLengthfromCenter=(max_axis-min_axis)*0.5
    center = (min_axis+max_axis)*0.5
    temp=np.dot(np.concatenate((center.reshape((1,-1)),np.diag(boxLengthfromCenter)),axis=0),tvect)
    return (temp[0],temp[1:],vect) #return center and lengthfrom center (vectors), rotation matrix to align box to x,y,z
def kabsch(from_a,to_b):
    '''
    returns rigid transformation of b=Ra+u, return R,u
    '''
    centroid_P=from_a.mean(axis=0)
    centroid_Q=to_b.mean(axis=0)
    P=from_a-centroid_P
    Q=to_b-centroid_Q
    u,s,vh=np.linalg.svd(P.T.dot(Q))
    d=np.sign(np.linalg.det(vh.T.dot(u.T)))
    R=vh.T.dot(np.array([[1,0,0,],[0,1,0],[0,0,d]]).dot(u.T))
    u=centroid_Q-R.dot(centroid_P.reshape((-1,1))).reshape(-1)
    return (R,u)
def pointsMeansqDiff_kabsch(from_a,to_b,averageFunc=None):
    R,u=kabsch(from_a,to_b)
    new_a=R.dot(from_a.T).T+u.reshape((1,-1))
    if averageFunc is None:
        return np.mean(np.sum((new_a-to_b)**2.,axis=1))
    else:
        return averageFunc(np.sum((new_a-to_b)**2.,axis=1))
class multicpu_cmd:
    def __init__(self,folder,insert_path_for_heartFEM=None):
        self.folder=folder
        self.insert_path_for_heartFEM=insert_path_for_heartFEM
    def __call__(self,sample):
        if isinstance(self.insert_path_for_heartFEM,str):
            cmd=[self.insert_path_for_heartFEM]
        elif self.insert_path_for_heartFEM is None:
            cmd=[]
        else:
            cmd=self.insert_path_for_heartFEM[::-1]
        cmd.insert(0,str(sample))
        cmd.insert(0,self.folder+"/multicpurun.py")
        cmd.insert(0,'python3')
        tool_subprocess = subprocess.run(cmd)
            
def multicpu_getLVbehavior(folderPath,step):
    import pickle
    with open(folderPath+"/multipcu_LVBehavior/"+'getLVBehavior_dictionary.pkl', 'rb') as f:
        kwargs = pickle.load(f)
    recreateKwargs={}
    behaviorkwargs={}
    for key in kwargs:
        if key[:9]=='LVclosed_':
            recreateKwargs[key[9:]]=kwargs[key]
        else:
            behaviorkwargs[key]=kwargs[key]
    runParameters=heartParameters.heartParameters()
    runParameters.readParameters(folderPath+"/runParameters.txt")

    LVclosed_object=LVclosed(defaultParameters=runParameters,defaultAge=recreateKwargs['LVage'],heartModel=recreateKwargs['heartModelStr'])
    for key in recreateKwargs:
        setattr(LVclosed_object,key, recreateKwargs[key])
    result=LVclosed_object.getLVbehavior(multicpu_step=step,**behaviorkwargs)
    np.save(folderPath+"/multipcu_LVBehavior/"+str(step)+".npy", result, allow_pickle=False)
    return;
def multicpu_iterativerun(folderPath,step):
    import pickle
    with open(folderPath+"/multipcu_iterative/"+'iterativerun_dictionary.pkl', 'rb') as f:
        kwargs = pickle.load(f)
    recreateKwargs={}
    behaviorkwargs={}
    for key in kwargs:
        if key[:9]=='LVclosed_':
            recreateKwargs[key[9:]]=kwargs[key]
        else:
            behaviorkwargs[key]=kwargs[key]

    LVclosed_object=LVclosed(defaultParameters=None,defaultAge=recreateKwargs['LVage'],heartModel=recreateKwargs['heartModelStr'])
    for key in recreateKwargs:
        if key=="casename":
            setattr(LVclosed_object,key, recreateKwargs[key]+"/"+str(step))
        else:
            setattr(LVclosed_object,key, recreateKwargs[key])
    if LVclosed_object.manualVolFile is not None:
        LVclosed_object.manualVol=LVclosed_object.getVolumeFromFile
    timelist=np.loadtxt(folderPath+"/multipcu_iterative/"+str(step)+"/runtimelist.txt")
    result=LVclosed_object.iterativeRun(editParameters={"read_from_file":folderPath+"/multipcu_iterative/"+str(step)},runTimeList=timelist.reshape(-1),**behaviorkwargs)
    return;
class LVclosed:
    def __new__(cls, *args, **kwargs):
        return super().__new__(cls)
    def __init__(self,defaultParameters=None,defaultAge='fetal28',heartModel='lcleeHeart'):
        if defaultParameters is None:
            defaultParameters={}
        self.defaultRunMode='iterativeRun'
        self.defaultRun_kwargs=None
        #define defaults
        self.meshname = "t0"
        self.casename = os.getcwd()
        self.numofCycle=12
        self.LVage=defaultAge
        self.heartMotionFile=None
        self.heartModelStr=heartModel
        if 'path of folder with case' in defaultParameters:
            self.casename=defaultParameters['path of folder with case']
        if 'name of mesh' in defaultParameters:
            self.meshname=defaultParameters['name of mesh']
        self.defaultParameters=heartParameters.heartParameters(defaultParameters=defaultParameters,defaultAge=defaultAge)
        if heartModel=='heArt.BiV':
            if 'Windkessel LV source function' not in defaultParameters:
                self.defaultParameters['Windkessel LV source function']='pressure 3D table'
            if 'Windkessel RV source function' not in defaultParameters:
                self.defaultParameters['Windkessel RV source function']='pressure 3D table'
        elif heartModel=='AVHeart':
            if 'Windkessel LV source function' not in defaultParameters:
                self.defaultParameters['Windkessel LV source function']='pressure 3D table'
            if 'Windkessel LA source function' not in defaultParameters:
                self.defaultParameters['Windkessel LA source function']='pressure 3D table'
        #self.heartModel=heartModelTranslator.HeartModel(self.heartModelStr,self.defaultParameters)
        self.manualVolFile=None
        self.manualVolFiletimecol=None
        self.manualVolFilevolcol=None
        self.manualVolFiletimescale=None
        self.manualVolFilevolscale=None
        self.manualVol=None
        self.manualPhaseTime=None
        
        self.constrains=['Kspring','translation','rotation']
        self.runCount=0
    def paraToRecreateObj(self):
        para ={'defaultRunMode':self.defaultRunMode,
               'defaultRun_kwargs':self.defaultRun_kwargs,
               'meshname':self.meshname,
               'casename':self.casename,
               'numofCycle':self.numofCycle,
               'LVage':self.LVage,
               'heartMotionFile':self.heartMotionFile,
               'heartModelStr':self.heartModelStr,
               'manualVolFile':self.manualVolFile,
               'manualVolFiletimecol':self.manualVolFiletimecol,
               'manualVolFilevolcol':self.manualVolFilevolcol,
               'manualVolFiletimescale':self.manualVolFiletimescale,
               'manualVolFilevolscale':self.manualVolFilevolscale,
               'manualPhaseTime':self.manualPhaseTime,
               'constrains':self.constrains,
               'runCount':self.runCount}
        return para
    def setHeartrate(self,beatsperminute):
        self.defaultParameters.setHeartrate(beatsperminute)
    def setHeartMotionFile(self,filename,timetotimestepScale=1.,timetotimestepShift=0.):
        self.timetotimestepScale=timetotimestepScale
        self.timetotimestepShift=timetotimestepShift
        self.heartMotionFile=filename
    def evaluateSTLCoords(self,coords,timestep):
        b=bsf.BsplineFourier(self.heartMotionFile)
        newcoords=np.zeros((coords.shape[0],4))
        newcoords[:,:3]=coords
        newcoords[:,3]=timestep
        return b.getCoordFromRef(newcoords)
    def getVolumeFromFile(self,time):
        if isinstance(self.manualVolFile,list):
            data=[]
            for n in range(len(self.manualVolFile)):
                try:
                    data_temp=np.loadtxt(self.manualVolFile[n])
                except:
                    data_temp=np.loadtxt(self.manualVolFile[n],skiprows=1)
                data.append(data_temp)
        elif isinstance(self.manualVolFile,dict):
            data=[]
            for chamber in ['LV','RV','LA','RA','Mitral','Aortic','Tricuspid','Pulmonic']:
                if chamber in self.manualVolFile:
                    try:
                        data_temp=np.loadtxt(self.manualVolFile[chamber])
                    except:
                        data_temp=np.loadtxt(self.manualVolFile[chamber],skiprows=1)
                    data.append(data_temp)
        else:
            try:
                data=[np.loadtxt(self.manualVolFile)]
            except:
                data=[np.loadtxt(self.manualVolFile,skiprows=1)]
        timedata=[]
        voldata=[]
        for n in range(len(data)):
            timedata.append(data[n][...,self.manualVolFiletimecol]*self.manualVolFiletimescale)
            voldata.append(data[n][...,self.manualVolFilevolcol]*self.manualVolFilevolscale)
        if isinstance(time,(int,float)):
            singleValued=True
            time=np.array([time])
        else:
            singleValued=False
            time=time.copy()
        return_result=[]
        
        for m in range(len(data)):
            for n in range(len(time)):
                while time[n]>timedata[m][-1]:
                    time[n]-=self.defaultParameters['duration of one cardiac cycle in ms']
                while time[n]<timedata[m][0]:
                    time[n]+=self.defaultParameters['duration of one cardiac cycle in ms']
            spl=scipyinterpolate.splrep(timedata[m], voldata[m])
            if singleValued:
                return_result.append( scipyinterpolate.splev(time, spl)[0])
            else:
                return_result.append( scipyinterpolate.splev(time, spl))
        if singleValued:
            return np.array(return_result)
        else:
            return np.array(return_result).T 
    def setManualVolumeFile(self,volfile,timecol=0,timescale=1.,volcol=1,volscale=1.):
        self.manualVolFile=volfile
        self.manualVolFiletimecol=timecol
        self.manualVolFilevolcol=volcol
        self.manualVolFiletimescale=timescale
        self.manualVolFilevolscale=volscale
        self.manualVol=self.getVolumeFromFile
        
    def getphaseFromFile(self,time):
        try:
            data=np.loadtxt(self.manualPhaseTimeFile)
        except:
            data=np.loadtxt(self.manualPhaseTimeFile,skiprows=1)
        timedata=data[:,self.manualPhaseTimeFiletimecol]*self.manualPhaseTimeFiletimescale
        phaseTimedata=data[:,self.manualPhaseTimeFilephaseTimecol]*self.manualPhaseTimeFilephaseTimescale
        while time>timedata[-1]:
            time-=self.defaultParameters['duration of one cardiac cycle in ms']
        if time<timedata[0]:
            raise Exception('time '+str(time)+' , is not inside file given:'+self.manualVolFile)
        spl=scipyinterpolate.splrep(timedata, phaseTimedata)
        return scipyinterpolate.splev(np.array([time]), spl)[0]
    def setManualPhaseTimeFile(self,phaseTimefile,timecol=0,timescale=1.,phaseTimecol=1,phaseTimescale=1.):
        self.manualPhaseTimeFile=phaseTimefile
        self.manualPhaseTimeFiletimecol=timecol
        self.manualPhaseTimeFilephaseTimecol=phaseTimecol
        self.manualPhaseTimeFiletimescale=timescale
        self.manualPhaseTimeFilephaseTimescale=phaseTimescale
        self.manualPhaseTime=self.getphaseFromFile
    def setDefaultWindkessel(self,modelString):
        self.defaultParameters.setDefaultWindkessel(modelString)
        if modelString[:5]=='fetal':
            try:
                self.setRVcurrentInflowFourier(self.casename+'/'+self.meshname+'_rvflowrate.txt')
            except Exception as e:
                logger.warning(e)
                logger.warning('RV fourier not loaded '+self.casename+'/'+self.meshname+'_rvflowrate.txt')
    def setPulse(self,cavity,amp,peaktime,width):
        #atrial_side = "r" or "l" or 'lr" or "rl"'A peak pressure','A time at peak pressure','A pressure pulse width'
        self.defaultParameters.setParameter('Windkessel '+cavity+' source function','pulse')
        self.defaultParameters.setParameter('Windkessel '+cavity+' source pulse function peak pressure in mmHg',amp)
        self.defaultParameters.setParameter('Windkessel '+cavity+' source pulse function peak pressure in mmHg',peaktime)
        self.defaultParameters.setParameter('Windkessel '+cavity+' source pulse function peak pressure in mmHg',width)
    def setLApulse(self,amp,peaktime,width):
        self.setPulse('LA',amp,peaktime,width)
    def setRApulse(self,amp,peaktime,width):
        self.setPulse('RA',amp,peaktime,width)
    def setLVpulse(self,amp,peaktime,width):
        self.setPulse('LV',amp,peaktime,width)
    def setRVpulse(self,amp,peaktime,width):
        self.setPulse('RV',amp,peaktime,width)
    def setCurrentInflowFourier(self,cavity,filename,fourierTerms=1,fineAdjustPhase=False):
        #outflow from RV
        fourierTerms=int(fourierTerms)
        data=np.loadtxt(filename)
        period=2.*data[-1,0]-data[-2,0]
        def fourierFit(x, *a):
            ret = 0.
            for deg in range(fourierTerms):
                ret += a[deg*2] * np.sin((deg+1) *2.* np.pi / period * x+a[deg*2+1])
            return ret
        popt, pcov = curve_fit(fourierFit, data[:,0], data[:,1],p0=np.ones(2*fourierTerms))
        if fineAdjustPhase:
            try_x=0.
            adj_x=period/30.
            y=fourierFit(try_x, *popt)
            if y>0:
                tune=1
            else:
                tune=-1
            while y!=0 and adj_x>(period*10**-6.):
                if y*tune<0:
                    tune*=-1
                    adj_x*=0.3
                try_x+=tune*adj_x
                y=fourierFit(try_x, *popt)
            adjust_x=-try_x
            logger.info("Adjusted "+cavity+"currentInflowFourier by "+str(adjust_x)+" radian")
        else:
            adjust_x=0.
        self.defaultParameters['Windkessel '+cavity+' source function']='fourier current'
        for n in range(4):
            self.defaultParameters['Windkessel '+cavity+' source fourier function sine amplitude term '+str(n)]=0.
            self.defaultParameters['Windkessel '+cavity+' source fourier function sine degree phase term '+str(n)]=0.
        for n in range(len(popt)):
            popt_temp=popt[n]
            if n%2==1:
                if adjust_x!=0:
                    popt_temp=popt_temp+(int(n/2.)+1) *2.* np.pi / period * adjust_x
                while popt_temp<-np.pi:
                    popt_temp+=2.*np.pi
                while popt_temp>np.pi:
                    popt_temp-=2.*np.pi
                popt_temp*=180./np.pi
                self.defaultParameters['Windkessel '+cavity+' source fourier function sine degree phase term '+str(int(n/2))]=popt_temp
            else:
                self.defaultParameters['Windkessel '+cavity+' source fourier function sine amplitude term '+str(int(n/2))]=popt_temp
                
    def setLAcurrentInflowFourier(self,filename,fourierTerms=1,fineAdjustPhase=False):
        self.setCurrentInflowFourier('LA',filename,fourierTerms=fourierTerms,fineAdjustPhase=fineAdjustPhase)
    def setLVcurrentInflowFourier(self,filename,fourierTerms=1,fineAdjustPhase=False):
        self.setCurrentInflowFourier('LV',filename,fourierTerms=fourierTerms,fineAdjustPhase=fineAdjustPhase)
    def setRAcurrentInflowFourier(self,filename,fourierTerms=1,fineAdjustPhase=False):
        self.setCurrentInflowFourier('RA',filename,fourierTerms=fourierTerms,fineAdjustPhase=fineAdjustPhase)
    def setRVcurrentInflowFourier(self,filename,fourierTerms=1,fineAdjustPhase=False):
        self.setCurrentInflowFourier('RV',filename,fourierTerms=fourierTerms,fineAdjustPhase=fineAdjustPhase)
    def scaleWinkessel(self,scaleDict,compstr=''):
        self.defaultParameters.scaleWinkessel(scaleDict,compstr=compstr)
    def scaleWinkesselwithAge(self,ageInWeeks,poweradjustDict=None,compstr=''):
        self.defaultParameters.scaleWinkesselwithAge(ageInWeeks,poweradjustDict=poweradjustDict,compstr=compstr)
    def resetCount(self,count=0):
        self.runCount=count
    def addConstrain(self,constrain_strings):
        if isinstance(constrain_strings,str):
            constrain_strings=[constrain_strings]
        for string in constrain_strings:
            if string not in self.constrains:
                self.constrains+=[string]
    def removeConstrain(self,constrain_strings):
        if isinstance(constrain_strings,str):
            constrain_strings=[constrain_strings]
        for string in constrain_strings:
            if string in self.constrains:
                self.constrains.remove(string)
    def changeDefaultParameters(self,editParameters=None):
        if editParameters is None:
            return
        else:
            for key in self.defaultParameters:
                if key in editParameters:
                    self.defaultParameters[key]=editParameters[key]
    def set_EVP_with_LA_min_volume(self,LAminVol):
        ageInWeeks=float(self.LVage[5:])
        LA_unloadedVolume=self.defaultParameters.getFetalPopulation_LA_unloadedVolume(ageInWeeks)
        self.defaultParameters['LV end diastolic pressure in mmHg']=(LAminVol-LA_unloadedVolume)/self.defaultParameters['lac']
        logger.info('with age'+repr(ageInWeeks)+' ,set EDP as '+repr(self.defaultParameters['LV end diastolic pressure in mmHg'])+' mmHg')
    def generateMesh(self,params,toRunCountFolder=False):
        params=params+{'path of folder with case':self.casename,'filename of stl':self.meshname}
        if toRunCountFolder:
            if isinstance(toRunCountFolder,str):
                if toRunCountFolder[0]!='/':
                    toRunCountFolder='/'+toRunCountFolder
            else:
                toRunCountFolder='/'+str(self.runCount)
            params=params+{'subfolder to save files':toRunCountFolder}
            self.defaultParameters+=heartModelTranslator.generateMesh(self.heartModelStr,params)
        else:
            self.defaultParameters+=heartModelTranslator.generateMesh(self.heartModelStr,params)
    def runWindkessel(self,comm,tstep,dt_dt,*args):
        raise Exception("Running Windkessel while doing FEA is depreciated")
        lvufile=self.casename+"/"+str(self.runCount)+'/LVgenerator.txt'
        lvudata=np.loadtxt(self.casename+"/"+str(self.runCount)+"/PV_.txt")
        lvucontroldata=np.loadtxt(self.casename+'/'+self.meshname+'_lvucontrol.txt')
        lvudata2=np.concatenate((lvucontroldata[:-1],lvudata[:,2]+np.array([[lvucontroldata[-1,0],0]])),axis=0)
        np.savetxt(lvufile,lvudata2)
        ngspice_py.simLVcircuit(self.casename+'/'+self.meshname,lvucontroldata[-1,0]+tstep+dt_dt,lvufile)
        self.manualVolFile=self.casename+'/'+self.meshname+'_circuit.txt'
        V_cav=self.getVolumeFromFile(lvucontroldata[-1,0]+tstep+dt_dt)
        return V_cav
    def solveLoadActive(self,heart,num=50):
        Tmax=heart.get_ReferenceFiberStrength()
        heart.set_CavityVolume(heart.get_CavityVolume())
        for n in np.linspace(0.,1.,num=num):
            heart.set_ReferenceFiberStrength(Tmax*n)
            heart.solve()
        return 1
    def solveVolume(self,heart,targetVolume,voladj=0.05):
        if isinstance(targetVolume,(int,float)):
            targetVolume=[targetVolume]
        targetVolume=np.array(targetVolume)
        tuneVol=[]
        for adjInd in range(len(targetVolume)):
            if heart.get_CavityVolume()[adjInd]>(targetVolume[adjInd]):
            	tuneVol.append(-1)
            elif heart.get_CavityVolume()[adjInd]<(targetVolume[adjInd]):
                tuneVol.append(1)
            else:
                tuneVol.append(0)
        tuneVol=np.array(tuneVol)
        prev_vol=targetVolume*0.
        if logger.isEnabledFor(logging.DEBUG):
            logger.debug("Starting PLV = "+repr( heart.get_CavityPressure())+" VLV = "+repr(heart.get_CavityVolume()))
            logger.debug("Target vol = "+repr(targetVolume))
        while np.any((heart.get_CavityVolume()*tuneVol)<(targetVolume*tuneVol)):
            pressadjust_voladj=np.maximum(0.01,0.75**(np.log(np.maximum(1.,np.abs(heart.get_CavityPressure())/10.))/np.log(2.)))
            vol_all=heart.get_CavityVolume()
            for adjInd in range(len(targetVolume)):
                if (((1+tuneVol[adjInd]*voladj*pressadjust_voladj[adjInd])*heart.get_CavityVolume()[adjInd])*tuneVol[adjInd])<(targetVolume[adjInd]*tuneVol[adjInd]):
                    vol_all[adjInd]=vol_all[adjInd]* (1+tuneVol[adjInd]*voladj*pressadjust_voladj[adjInd])
                else:
                    vol_all[adjInd]=targetVolume[adjInd]
            heart.set_CavityVolume(vol_all)
            if logger.isEnabledFor(logging.DEBUG):
                logger.debug("Start solve with PLV = "+repr( heart.get_CavityPressure())+" VLV = "+repr(heart.get_CavityVolume()))
                logger.debug("Solve for vol = "+repr(vol_all))
            heart.solve()
            p_cav = heart.get_CavityPressure()
            V_cav = heart.get_CavityVolume()
            if isinstance(p_cav,(float,int)):
                logger.info("PLV = "+repr(p_cav)+" VLV = "+repr(V_cav))
            elif isinstance(p_cav,(list,np.ndarray)):
        	    for n in range(len(p_cav)):
        	        if not(isinstance(p_cav[n],(float,int))):
        	            raise Exception('Unable to converge to volume '+repr(targetVolume))
        	        else:
        	            logger.info("Chamber "+str(n)+", PLV = "+repr(p_cav[n])+" VLV = "+repr(V_cav[n]))
            else:
                raise Exception('Unable to converge to volume '+repr(targetVolume))
            if np.all(np.abs(V_cav/targetVolume-1.)<10**-4):
                break
            if np.all(np.abs((V_cav-prev_vol)/targetVolume)<10**-4) :
                logger.warning("unable to reach target "+repr(targetVolume))
                break
            prev_vol=V_cav
        return 1
    def solveTa(self,heart,targetTa,peaktime=None):
        if heart.get_activeTime()>(targetTa):
        	tuneTa=-1
        elif heart.get_activeTime()<(targetTa):
        	tuneTa=1
        else:
        	tuneTa=0
        while (heart.get_activeTime()*tuneTa)<(targetTa*tuneTa):
        	t=heart.get_activeTime()
        	if (t >= 0.0 and t < 4.0):
        		dt = 0.50
        	elif peaktime is not None:
        		if (t >= (peaktime-0.5) and t < (peaktime+0.5)):
        		    dt = 3.
        		else:
        		    dt = 1.0
        	else :
        		dt = 1.0
        	if ((t+dt*tuneTa)*tuneTa)<(targetTa*tuneTa):
        		heart.set_activeTime(t+dt*tuneTa)
        	else:
        		heart.set_activeTime(targetTa)
            
        	heart.set_dTime(abs(heart.get_activeTime()-t))
        	logger.info("heart.activeTime = "+repr(heart.get_activeTime())+" heart.dTime = "+repr(heart.get_dTime()))
        	heart.solve()
        	p_cav = heart.get_CavityPressure()
        	V_cav = heart.get_CavityVolume()
        	if isinstance(p_cav,(float,int)):
        		logger.info("PLV = "+repr(p_cav)+" VLV = "+repr(V_cav))
        	elif isinstance(p_cav,(list,np.ndarray)):
        	    for n in range(len(p_cav)):
        	        if not(isinstance(p_cav[n],(float,int))):
        	            raise Exception('Unable to converge to active time '+repr(targetTa))
        	        else:
        	            logger.info("Chamber "+str(n)+", PLV = "+repr(p_cav[n])+" VLV = "+repr(V_cav[n]))
        	else:
        		raise Exception('Unable to converge to active time '+repr(targetTa))  
        return 1
    def solvePressure(self,heart,targetPressure,voladj=0.05,minVolumeBreak=0.,maxVolumeBreak=float('inf'),fixed_volume=None):
        
        
        if isinstance(targetPressure,(int,float)):
            targetPressure=[targetPressure]
        if isinstance(minVolumeBreak,(int,float)):
            minVolumeBreak=[minVolumeBreak]*len(targetPressure)
        if isinstance(maxVolumeBreak,(int,float)):
            maxVolumeBreak=[maxVolumeBreak]*len(targetPressure)
        if fixed_volume is None:
            fixed_volume=np.array([0.]*len(targetPressure))
        #target pressure in mmHg
        tuneVol=[]
        for n in range(len(targetPressure)):
            if heart.get_CavityPressure()[n]>(targetPressure[n]):
            	tuneVol.append(-1)
            elif heart.get_CavityPressure()[n]<(targetPressure[n]):
                tuneVol.append(1)
            else:
                tuneVol.append(0)
        tuneVol=np.array(tuneVol)
        for n in range(len(targetPressure)):
            if fixed_volume[n]!=0:
                if heart.get_CavityVolume()[n]>fixed_volume[n]:
                    tuneVol[n]=-1.
                elif heart.get_CavityVolume()[n]<fixed_volume[n]:
                    tuneVol[n]=1.
                else:
                    tuneVol[n]=0.
        if isinstance(voladj,(int,float)):
            voladj=[voladj]*len(targetPressure)
        voladj=np.array(voladj)
        heart.set_CavityVolume(heart.get_CavityVolume())
        count=0
        if np.any(targetPressure<0):
            addcount=-1
        else:
            addcount=1
        if logger.isEnabledFor(logging.DEBUG):
            debugout=0
            while os.path.isdir(self.casename+"/debug"+str(debugout)):
                debugout+=1
            debugWriter=fenicsResultWriter(self.casename+"/debug"+str(debugout),["displacement","strain_energy_density"])
        while True:#np.any((heart.get_CavityPressure()*tuneVol)<(targetPressure*tuneVol)):
            if np.any(heart.get_CavityVolume()<minVolumeBreak) or np.any(heart.get_CavityVolume()>maxVolumeBreak):
                return 0
            setVol=heart.get_CavityVolume()*(1+tuneVol*voladj)
            for n in range(len(targetPressure)):
                if fixed_volume[n]!=0 and abs(setVol[n]/fixed_volume[n]-1)<voladj[n]:
                    setVol[n]=fixed_volume[n]
            heart.set_CavityVolume(setVol)
            if logger.isEnabledFor(logging.DEBUG):
                logger.debug("Start solve with PLV = "+repr( heart.get_CavityPressure())+" VLV = "+repr(heart.get_CavityVolume()))
                logger.debug("Solve for vol = "+repr(setVol))
                if count==0:
                    logger.debug(repr(heart.set_debug_call(debugWriter,heart)))
                    logger.debug(repr(heart.set_debug_call(heart.get_CavityVolume,True)))
                    logger.debug(repr(heart.set_debug_call(heart.get_CavityPressure,True)))
            heart.solve()
            #debugWriter(heart)##!!!!debug
            p_cav = heart.get_CavityPressure()
            V_cav = heart.get_CavityVolume()
            if isinstance(p_cav,(float,int)):
                logger.info("PLV = "+repr(p_cav)+" VLV = "+repr(V_cav))	
            elif isinstance(p_cav,(list,np.ndarray)):
        	    for n in range(len(p_cav)):
        	        if not(isinstance(p_cav[n],(float,int))):
        	            raise Exception('Unable to converge to pressure '+repr(targetPressure)+' mmHg')
        	        else:
        	            logger.info("Chamber "+str(n)+", PLV = "+repr(p_cav[n])+" VLV = "+repr(V_cav[n]))
            else:
                raise Exception('Unable to converge to pressure '+repr(targetPressure)+' mmHg')
            for n in range(len(targetPressure)):
                if fixed_volume[n]!=0:
                    if heart.get_CavityVolume()[n]>fixed_volume[n]:
                        tuneVol[n]=-1.
                    elif heart.get_CavityVolume()[n]<fixed_volume[n]:
                        tuneVol[n]=1.
                    else:
                        tuneVol[n]=0.
                else:
                    if (heart.get_CavityPressure()[n]*tuneVol[n])>(targetPressure[n]*tuneVol[n]):
                        voladj[n]=voladj[n]*0.33
                        tuneVol[n]=tuneVol[n]*-1
                        logging.debug('Change tuneVol '+str(n)+' to '+repr(tuneVol[n])+'  voladj='+repr(voladj[n]))
                    if tuneVol[n]==0 and heart.get_CavityPressure()[n]!=(targetPressure[n]):
                        if heart.get_CavityPressure()[n]>(targetPressure[n]):
                            tuneVol[n]=-1
                        else:
                            tuneVol[n]=1
            if np.abs(targetPressure[fixed_volume==0]).min()/np.abs(targetPressure[fixed_volume==0]).max() < 10**-2:
                temp_pressure_base=np.abs(targetPressure[fixed_volume==0]).max()
                temp_condition=np.all(np.abs(heart.get_CavityPressure()[fixed_volume==0]-targetPressure[fixed_volume==0])/temp_pressure_base<10**-4) 
            else:
                temp_condition=np.all(np.abs(heart.get_CavityPressure()[fixed_volume==0]/targetPressure[fixed_volume==0]-1.)<10**-4)    
            if temp_condition and np.all(np.abs(heart.get_CavityVolume()[fixed_volume>0]/fixed_volume[fixed_volume>0]-1.)<10**-4):
                break
            elif voladj.max()<10**-4:
                break
            if count>10000:
                raise Exception('Due to solving for negative pressure, Maximum (10000) iterations is set and reached. Unable to converge.')
            else:
                count+=addcount
        return 1
    def getLVbehavior(self,runParameters=None,meshname=None,meshFolder=None,volstep=50,extendVol=0.1,minESV=None,maxEDV=None,saveaddstr='',voladj=0.1,starttime=0,endtime=None,outputResultList=None,toRunCountFolder=False,multicpu_step=None,returnVolspaceOnly=False,saveFile=False):
        if toRunCountFolder:
            if isinstance(toRunCountFolder,str):
                if toRunCountFolder[0]!="/":
                    toRunCountFolder="/"+toRunCountFolder
            else:
                toRunCountFolder="/"+str(self.runCount)
        else:
            toRunCountFolder=""
        if meshFolder is None:
            meshFolder=''
        elif meshFolder[0]!='/':
            meshFolder='/'+meshFolder
        if runParameters is None:
            runParameters=self.defaultParameters
            trycontinueprevious=False
        elif isinstance(runParameters,str):
            if runParameters[0]!='/':
                runParameters='/'+runParameters
            runParameters=self.readRunParameters(self.casename+runParameters)
            trycontinueprevious=True
        else:
            trycontinueprevious=False
        if meshname is None:
            meshname=self.meshname
        if self.heartModelStr=='heArt.BiV':
            project_form_compiler_parameters={"representation":"uflacs","quadrature_degree":4}
        else:
            project_form_compiler_parameters={"representation":"uflacs"}
        if outputResultList is None:
            outputResultList=[]
        elif isinstance(outputResultList,str):
            outputResultList=[outputResultList]
        
        heart=heartModelTranslator.HeartModel(self.heartModelStr,runParameters+{'path of folder with case':self.casename+meshFolder,'name of mesh':meshname})
        
        
        ########################### Fenics's Newton  #########################################################
        heart.set_CavityVolume(heart.get_CavityVolume())

        cavityNames=heart.get_CavityNames()
        resultWriter=fenicsResultWriter(self.casename+toRunCountFolder,outputResultList,project_form_compiler_parameters=project_form_compiler_parameters,tag=["time"]+["vol:"+x for x in cavityNames])
        # Closed loop cycle
        BCL = runParameters['duration of one cardiac cycle in ms']
        
        
        #######set volume solve space############# 
        EDV={}
        for cavity in cavityNames:
            if cavity+' end diastolic volume in mL' in runParameters:
                EDV[cavity]=runParameters[cavity+' end diastolic volume in mL']
        if maxEDV is None:
            maxEDV={}
            for cavity in cavityNames:
                maxEDV[cavity]=EDV[cavity]*(1+extendVol)
        ESV={}
        for cavity in cavityNames:
            if cavity+' end systolic volume in mL' in runParameters:
                ESV[cavity]=runParameters[cavity+' end systolic volume in mL']
            else:
                ESV[cavity]=getattr(heart,'get_'+cavity+'Volume')()
        if minESV is None:
            minESV={}
            for cavity in cavityNames:
                minESV[cavity]=ESV[cavity]*(1-extendVol)
        
        if isinstance(volstep,(list,np.ndarray)):
            volstep=np.sort(volstep)
            volumeSpace=np.zeros((0,len(volstep)))
            for cavity in cavityNames:
                volumeSpace=np.concatenate((volumeSpace,[minESV[cavity]+(maxEDV[cavity]-minESV[cavity])*volstep[::-1]]),axis=0)
        elif isinstance(volstep,dict):
            volstep_temp=volstep
            volstep={}
            volstep_len=None
            for key in volstep_temp:
                volstep[key]=np.sort(volstep_temp[key])
                if volstep_len is None:
                    volstep_len=len(volstep[key])
                elif volstep_len!=len(volstep[key]):
                    raise Exception("Length of volstep for different cavity have to be the same.")
            volumeSpace=np.zeros((0,volstep_len))
            for cavity in cavityNames:
                volumeSpace=np.concatenate((volumeSpace,[minESV[cavity]+(maxEDV[cavity]-minESV[cavity])*volstep[cavity][::-1]]),axis=0)
        else:
            volumeSpace=np.zeros((0,volstep))
            for cavity in cavityNames:
                volumeSpace=np.concatenate((volumeSpace,[minESV[cavity]*(maxEDV[cavity]/minESV[cavity])**(np.linspace(0,1,num=volstep))[::-1]]),axis=0)

        
        ####### Loading phase to MAX volume ####################################################
        #self.solveVolume(heart,volumeSpace[0])

        timeSpace=[]
        tstep=0
        while(tstep < 2*BCL):
            timeSpace.append(tstep)
            if(tstep >= 0.0 and tstep < 4.0):
           		tstep += 0.50
            elif (tstep >= (runParameters['time to maximum LV fiber tension in ms']-0.5) and tstep < (runParameters['time to maximum LV fiber tension in ms']+0.5)):
           		tstep += 3
            else :
           		tstep += 1.0
        timeSpace=np.array(timeSpace)
        timeSpace=timeSpace[timeSpace>=starttime]
        if endtime is not None:
            if np.count_nonzero(timeSpace<=endtime)<1:
                timeSpace=timeSpace[:1]
            else:
                timeSpace=timeSpace[timeSpace<=endtime]
        if returnVolspaceOnly:
            if saveFile:
                np.savetxt(self.casename+toRunCountFolder+"/"+meshname+"_Press_volumeSpace"+saveaddstr+".txt",np.array(volumeSpace))
                np.savetxt(self.casename+toRunCountFolder+"/"+meshname+"_Press_timeSpace"+saveaddstr+".txt",np.array(timeSpace))
            return (volumeSpace,timeSpace)
        if(fenics.MPI.rank(heart.get_Comm()) == 0):
            if len(timeSpace)>1 and (multicpu_step is None or multicpu_step==0):
                np.savetxt(self.casename+toRunCountFolder+"/"+meshname+"_Press_volumeSpace"+saveaddstr+".txt",np.array(volumeSpace))
                np.savetxt(self.casename+toRunCountFolder+"/"+meshname+"_Press_timeSpace"+saveaddstr+".txt",np.array(timeSpace))
        if len(timeSpace)==1:
            volPress=np.zeros((volumeSpace.shape[1]**volumeSpace.shape[0],2))
            volPress[:,0]=volumeSpace.copy()
            if abs(timeSpace[0])<=0.5:
                volPressStr='EDPVR'
            elif abs(timeSpace[0]-runParameters['time to maximum LV fiber tension in ms'])<=0.5:
                volPressStr='ESPVR'
            else:
                volPressStr='t'+str(int(timeSpace[0]))+'PVR'
        logger.info('========START============')
        
        press_volTime_base=[]
        press_volTime=[]
        for n in range(volumeSpace.shape[0]):
            press_volTime_base.append(np.zeros(volumeSpace.shape[1]**volumeSpace.shape[0]))
            press_volTime.append(np.zeros((volumeSpace.shape[1]**volumeSpace.shape[0],len(timeSpace))))
        startVolumeIndex=0
        
        if trycontinueprevious and (multicpu_step is None):
            try:
                volumeSpace2=np.loadtxt(self.casename+toRunCountFolder+"/"+meshname+"_Press_volumeSpace"+saveaddstr+".txt")
                timeSpace2=np.loadtxt(self.casename+toRunCountFolder+"/"+meshname+"_Press_timeSpace"+saveaddstr+".txt")
                if np.all(np.abs(volumeSpace2-volumeSpace)<(runParameters['LV end diastolic volume in mL']/(10.**3.))) and np.all(np.abs(timeSpace2-timeSpace)<0.01):
                    for n in range(volumeSpace.shape[0]):
                        press_volTime[n]=np.loadtxt(self.casename+toRunCountFolder+"/"+meshname+"_Press_VolTime_"+cavityNames[n]+saveaddstr+".txt")
                        press_volTime_base[n]=np.loadtxt(self.casename+toRunCountFolder+"/"+meshname+"_Press_VolTime_"+cavityNames[n]+"_base"+saveaddstr+".txt")
                    startVolumeIndex=np.nonzero(press_volTime_base==0.)[0][0]
                    volumeSpace=volumeSpace2
                    timeSpace=timeSpace2
            except:
                pass
        if multicpu_step is not None:
            startVolumeIndex=multicpu_step
            total_volume_to_calculate=1+startVolumeIndex
        else:
            total_volume_to_calculate=volumeSpace.shape[1]**volumeSpace.shape[0]
        
        for curr_volN in range(startVolumeIndex,total_volume_to_calculate):
            if curr_volN==0 or len(timeSpace)>1:
                heart=heartModelTranslator.HeartModel(self.heartModelStr,runParameters+{'path of folder with case':self.casename+meshFolder,'name of mesh':meshname})
                self.solveLoadActive(heart)
            press_volTime_temp=[]
            samePressure=0
            urInd=np.unravel_index(curr_volN, [volumeSpace.shape[1]]*volumeSpace.shape[0],order='F')
            currentVolValue=[]
            for m in range(volumeSpace.shape[0]):
                currentVolValue.append(volumeSpace[m,urInd[m]])
            self.solveVolume(heart,currentVolValue)
            for curr_timeN in range(len(timeSpace)):
                t = timeSpace[curr_timeN]
                logger.info('vol:'+repr(curr_volN)+'/'+repr(total_volume_to_calculate)+'  time:'+repr(curr_timeN)+'/'+repr(len(timeSpace))+" t_a="+repr(t))
                self.solveTa(heart,t,peaktime=runParameters['time to maximum LV fiber tension in ms'])
                logger.info("PLV = "+repr(heart.get_CavityPressure())+" VLV = "+repr(heart.get_CavityVolume()))
                press_volTime_temp.append(heart.get_CavityPressure())
                resultWriter(heart,tag=dict(zip(["time"]+["vol:"+x for x in cavityNames], [t]+currentVolValue)))
                if len(press_volTime_temp)>1:
                    if np.abs(press_volTime_temp[-1]-press_volTime_temp[-2]).max()<10**-6:
                        samePressure+=1
                    else:
                        samePressure=0
                    if samePressure>1:
                        break
            press_volTime_temp=np.array(press_volTime_temp).reshape((-1,volumeSpace.shape[0]))
            for n in range(volumeSpace.shape[0]):
                press_volTime_base[n][curr_volN]=press_volTime_temp[0,n]
                press_volTime[n][curr_volN,:len(press_volTime_temp[:,n])]=press_volTime_temp[:,n]-press_volTime_base[n][curr_volN]
            if len(timeSpace)==1:
                volPress[curr_volN,0]=heart.get_CavityVolume()
                volPress[curr_volN,1]=press_volTime_base[0][curr_volN]+press_volTime[0][curr_volN,0]
            if(fenics.MPI.rank(heart.get_Comm()) == 0):
                if (multicpu_step is None):
                    if len(timeSpace)>1:
                        for cavityN in range(len(cavityNames)):
                            np.savetxt(self.casename+toRunCountFolder+"/"+meshname+"_Press_VolTime"+"_"+cavityNames[cavityN]+saveaddstr+".txt",press_volTime[cavityN],header=str(curr_volN)+'solved rowVolm: '+str(total_volume_to_calculate)+'\ncolTime: '+str(len(timeSpace)))
                            np.savetxt(self.casename+toRunCountFolder+"/"+meshname+"_Press_VolTime"+"_"+cavityNames[cavityN]+"_base"+saveaddstr+".txt",press_volTime_base[cavityN])
                    else:
                        np.savetxt(self.casename+toRunCountFolder+"/"+meshname+"_"+volPressStr+saveaddstr+".txt",volPress,header="Volume, Pressure")

        if(fenics.MPI.rank(heart.get_Comm()) == 0):
            if (multicpu_step is None):
                if len(timeSpace)>1:
                    for cavityN in range(len(cavityNames)):
                        np.savetxt(self.casename+toRunCountFolder+"/"+meshname+"_Press_VolTime"+"_"+cavityNames[cavityN]+saveaddstr+".txt",press_volTime[cavityN],header=str(curr_volN)+'solved rowVolm: '+str(total_volume_to_calculate)+'\ncolTime: '+str(len(timeSpace)))
                        np.savetxt(self.casename+toRunCountFolder+"/"+meshname+"_Press_VolTime"+"_"+cavityNames[cavityN]+"_base"+saveaddstr+".txt",press_volTime_base[cavityN])
                else:
                    np.savetxt(self.casename+toRunCountFolder+"/"+meshname+"_"+volPressStr+saveaddstr+".txt",volPress,header="Volume, Pressure")
        if multicpu_step is not None:
            if len(timeSpace)>1:
                volPress=[]
                for cavityN in range(len(cavityNames)):
                    volPress.append(press_volTime_base[cavityN][multicpu_step]+press_volTime[cavityN][multicpu_step])
            else:
                volPress=volPress[multicpu_step]
            return np.array(volPress)
    def savehdf5(self,heart,casename,meshname,savename,saveFolder=''):
        #mesh = fenics.Mesh()
        f = fenics.HDF5File(fenics.MPI.comm_world, casename+'/'+meshname+".hdf5", 'r') 
        #f.read(mesh, meshname, False)
        facetboundaries = fenics.MeshFunction("size_t", heart.get_Mesh(), 2)
        f.read(facetboundaries, meshname+"/"+"facetboundaries")
        
        fiberFS = fenics.VectorFunctionSpace(heart.get_Mesh(), 'DG', 0)
        VQuadelem = fenics.VectorElement("Quadrature", heart.get_Mesh().ufl_cell(), degree=4, quad_scheme="default") #segmentation fault happened due to degree = 2. degree = 4 should be implied. I changed it. (Joy) 
        VQuadelem._quad_scheme = 'default'
        for e in VQuadelem.sub_elements():
        	e._quad_scheme = 'default'
        fiberFS = fenics.FunctionSpace(heart.get_Mesh(), VQuadelem)
        f0 = fenics.Function(fiberFS)
        s0 = fenics.Function(fiberFS)
        n0 = fenics.Function(fiberFS)
        c0 = fenics.Function(fiberFS)
        l0 = fenics.Function(fiberFS)
        r0 = fenics.Function(fiberFS)
        matid = fenics.MeshFunction('size_t',heart.get_Mesh(), 3, heart.get_Mesh().domains())
        facetboundaries = fenics.MeshFunction("size_t", heart.get_Mesh(), 2)
        edgeboundaries = fenics.MeshFunction("size_t", heart.get_Mesh(), 1)
        
        f.read(f0, meshname+"/"+"eF")
        f.read(s0, meshname+"/"+"eS")
        f.read(n0, meshname+"/"+"eN")
        f.read(c0, meshname+"/"+"eC")
        f.read(l0, meshname+"/"+"eL")
        f.read(r0, meshname+"/"+"eR")
        f.read(matid, meshname+"/"+"matid")
        f.read(facetboundaries, meshname+"/"+"facetboundaries")
        f.read(edgeboundaries, meshname+"/"+"edgeboundaries")
        
        cstrainFactor=None
        activefiberstrength=None
        activefiberdelay=None
        if f.has_dataset(meshname+"/"+"cstrainfactor"):
            cstrainFactor = fenics.Function(fiberFS)
            f.read(cstrainFactor, meshname+"/"+"cstrainfactor")
        if f.has_dataset(meshname+"/"+"activefiberstrength"):
            activefiberstrength = fenics.Function(fiberFS)
            f.read(activefiberstrength, meshname+"/"+"activefiberstrength")
        if f.has_dataset(meshname+"/"+"activefiberdelay"):
            activefiberdelay = fenics.Function(fiberFS)
            f.read(activefiberdelay, meshname+"/"+"activefiberdelay")
        
        
        matid_filename = casename+'/'+savename + "_matid.pvd"
        fenics.File(matid_filename) << matid
        f = fenics.HDF5File(heart.get_Mesh().mpi_comm(), casename+saveFolder+'/'+savename+".hdf5", 'w')
        f.write(heart.get_Mesh(), savename)
        f.close()
       	
        f = fenics.HDF5File(heart.get_Mesh().mpi_comm(), casename+saveFolder+'/'+savename+".hdf5", 'a') 
        f.write(facetboundaries, savename+"/"+"facetboundaries") 
        f.write(edgeboundaries, savename+"/"+"edgeboundaries") 
        f.write(matid, savename+"/"+"matid") 
        f.write(f0, savename+"/"+"eF") 
        f.write(s0, savename+"/"+"eS") 
        f.write(n0, savename+"/"+"eN")
        f.write(c0, savename+"/"+"eC") 
        f.write(l0, savename+"/"+"eL") 
        f.write(r0, savename+"/"+"eR")
       
        if cstrainFactor is not None:
            f.write(cstrainFactor, savename+"/"+"cstrainfactor")
        if activefiberstrength is not None:
            f.write(activefiberstrength, savename+"/"+"activefiberstrength")
        if cstrainFactor is not None:
            f.write(cstrainFactor, savename+"/"+"cstrainFactor")
        f.close()
       
       
        fenics.File(casename+saveFolder+'/'+savename + "_facetboundaries"+".pvd") << facetboundaries
        fenics.File(casename+saveFolder+'/'+savename + "_edgeboundaries"+".pvd") << edgeboundaries
        fenics.File(casename+saveFolder+'/'+savename + "_mesh" + ".pvd") << heart.get_Mesh()
        fenics.File(casename+saveFolder+'/'+savename + "_matid" +".pvd") << matid
        
    def adjustMeshToPressure(self,casename,meshname,targetPressure,savename=None,softAdjust=float('inf'),iterationNumber=float('inf'),initHeart=None,prevDeform=None,saveFolder='',tolerance=10.**-4.):
        #set usePrevious to False or 0 to not use previous mesh for deformation, set True to use previous mesh only and set a number to use previous mesh a number of times is it doesnt converge
        runParameters=self.defaultParameters.duplicate()
        runParameters['maximum LV fiber tension in Pa']=0.0
        refHeart=heartModelTranslator.HeartModel(self.heartModelStr,runParameters+{'path of folder with case':casename,'name of mesh':meshname})
        
        if initHeart is None:
            heart=heartModelTranslator.HeartModel(self.heartModelStr,runParameters+{'path of folder with case':casename,'name of mesh':meshname})
            prevDeformVectors=heart.get_DisplacementResult().vector()[:].copy()*0.
        else:
            heart=initHeart
            if prevDeform is None:
                prevDeformVectors=heart.get_DisplacementResult().vector()[:].copy()*0.
            else:
                prevDeformVectors=prevDeform
        targetVolume=refHeart.get_CavityVolume()
        
        prevMeshCoords=heart.get_MeshCoordinates().copy()
        self.solvePressure(heart,targetPressure,voladj=runParameters['solver volume adjustment'],maxVolumeBreak=heart.get_CavityVolume()*softAdjust)
        deformVectors=heart.get_DisplacementResult().vector()[:].copy()
        maxadjustmentsq=None
        newheartVol=heart.get_CavityVolume()
        prevHeartVol=heart.get_CavityVolume()
        
        temp_fiber_vec_file=runParameters['LV geometrical fiber direction variation file']
        count=0
        unloadresultWriter=fenicsResultWriter(self.casename+"/unloadsave_neg",["displacement"])
        while np.any(np.abs(heart.get_CavityPressure() /targetPressure-1)>tolerance) or np.any(np.abs(heart.get_CavityVolume()/targetVolume-1)>tolerance):##alternate convergence criteria##(np.mean((heart.get_MeshCoordinates()-du.vector()[:].reshape((-1,3))-refHeart.get_MeshCoordinates())**2.).max()**1.5/refHeart.get_CavityVolume())<10.**-6.:
            V = fenics.VectorFunctionSpace(heart.get_Mesh(), "CG", 1)
            du = fenics.project(heart.get_DisplacementResult(),V)
            if maxadjustmentsq is None:
                maxadjustmentsq=np.sum(du.vector()[:].reshape((-1,3))**2.,axis=1).max()**1.5/newheartVol
                ratio=1.
            else:
                ratio=min(1.,maxadjustmentsq/np.sum((heart.get_MeshCoordinates()+du.vector()[:].reshape((-1,3))-refHeart.get_MeshCoordinates())**2.,axis=1).max()**1.5/newheartVol)
            
            
            newheartVol=heart.get_CavityVolume()
            logging.info('Adjusting mesh... Max adjustment ='+repr(np.sum(du.vector()[:].reshape((-1,3)),axis=1).max()))
            logging.info('ratio='+repr(ratio))
            logging.info('target volume='+repr(targetVolume))
            logging.info('prevheart volume='+repr(prevHeartVol))
            logging.info('newheart volume ='+repr(heart.get_CavityVolume()))
            prevHeartVol=heart.get_CavityVolume()
            unloadresultWriter(heart)
            
            refHeart=heartModelTranslator.HeartModel(self.heartModelStr,runParameters+{'path of folder with case':casename,'name of mesh':meshname})         
            newDeform=refHeart.get_DisplacementResult().copy(deepcopy=True)
            newDeform.vector()[:]*=0.
            newDeform.vector()[:]-=heart.get_DisplacementResult().vector()[:]
            if isinstance(runParameters['LV geometrical fiber direction variation file'],str):
                temp_fiber_vec_file=runParameters['LV geometrical fiber direction variation file'][:-4]+'_tempadjmesh'+runParameters['LV geometrical fiber direction variation file'][-4:]
                fiber_temp=heart.get_StretchTensorFunction()
                temp_u=heart.get_DisplacementFunction()
                old_fiber=np.loadtxt(runParameters['LV geometrical fiber direction variation file'])
                temp_probes_f=np.empty((old_fiber.shape[0],9),dtype=float)
                temp_probes_u=np.empty((old_fiber.shape[0],3),dtype=float)
                fiber_temp=fenics.project(fiber_temp,fenics.TensorFunctionSpace(heart.get_Mesh(),'CG',2))
                fiber_temp.set_allow_extrapolation(True)
                temp_u=fenics.project(temp_u,fenics.VectorFunctionSpace(heart.get_Mesh(),'CG',2))
                temp_u.set_allow_extrapolation(True)
                for n in range(old_fiber.shape[0]):
                    fiber_temp.eval(temp_probes_f[n],old_fiber[n,:3])
                    temp_u.eval(temp_probes_u[n],old_fiber[n,:3])
                temp_probes_f=-temp_probes_f.reshape((-1,3,3))
                temp_probes_f=np.matmul(temp_probes_f,old_fiber[:,:3].reshape((-1,3,1))).reshape((-1,3))
                temp_probes_f=temp_probes_f/np.linalg.norm(temp_probes_f,axis=-1,keepdims=True)
                old_fiber[:,:3]=old_fiber[:,:3]-temp_probes_u
                old_fiber[:,3:]=temp_probes_f[:]
                np.savetxt(temp_fiber_vec_file,old_fiber)
            else:
                temp_fiber_vec_file=runParameters['LV geometrical fiber direction variation file']
            if isinstance(runParameters['LV geometrical fiber delay variation file'],str):
                temp_u=heart.get_DisplacementFunction()
                old_fiber=np.loadtxt(runParameters['LV geometrical fiber delay variation file'])
                temp_probes_u=np.empty((old_fiber.shape[0],3),dtype=float)
                temp_u=fenics.project(temp_u,fenics.VectorFunctionSpace(heart.get_Mesh(),'CG',2))
                temp_u.set_allow_extrapolation(True)
                for n in range(old_fiber.shape[0]):
                    temp_u.eval(temp_probes_u[n],old_fiber[n,:3])
                old_fiber[:,:3]=old_fiber[:,:3]-temp_probes_u
                temp_fiber_delay_file=runParameters['LV geometrical fiber delay variation file'][:-4]+'_tempadjmesh'+runParameters['LV geometrical fiber delay variation file'][-4:]
                np.savetxt(temp_fiber_delay_file,old_fiber)
            else:
                temp_fiber_delay_file=runParameters['LV geometrical fiber delay variation file']
            if isinstance(runParameters['LV geometrical fiber strength variation file'],str):
                temp_u=heart.get_DisplacementFunction()
                old_fiber=np.loadtxt(runParameters['LV geometrical fiber strength variation file'])
                temp_probes_u=np.empty((old_fiber.shape[0],3),dtype=float)
                temp_u=fenics.project(temp_u,fenics.VectorFunctionSpace(heart.get_Mesh(),'CG',2))
                temp_u.set_allow_extrapolation(True)
                for n in range(old_fiber.shape[0]):
                    temp_u.eval(temp_probes_u[n],old_fiber[n,:3])
                old_fiber[:,:3]=old_fiber[:,:3]-temp_probes_u
                temp_fiber_strength_file=runParameters['LV geometrical fiber strength variation file'][:-4]+'_tempadjmesh'+runParameters['LV geometrical fiber strength variation file'][-4:]
                np.savetxt(temp_fiber_strength_file,old_fiber)
            else:
                temp_fiber_strength_file=runParameters['LV geometrical fiber strength variation file']
            if isinstance(runParameters['LV geometrical overall stiffness constant variation file'],str):
                temp_u=heart.get_DisplacementFunction()
                old_fiber=np.loadtxt(runParameters['LV geometrical overall stiffness constant variation file'])
                temp_probes_u=np.empty((old_fiber.shape[0],3),dtype=float)
                temp_u=fenics.project(temp_u,fenics.VectorFunctionSpace(heart.get_Mesh(),'CG',2))
                temp_u.set_allow_extrapolation(True)
                for n in range(old_fiber.shape[0]):
                    temp_u.eval(temp_probes_u[n],old_fiber[n,:3])
                old_fiber[:,:3]=old_fiber[:,:3]-temp_probes_u
                temp_fiber_stffvar_file=runParameters['LV geometrical overall stiffness constant variation file'][:-4]+'_tempadjmesh'+runParameters['LV geometrical overall stiffness constant variation file'][-4:]
                np.savetxt(temp_fiber_stffvar_file,old_fiber)
            else:
                temp_fiber_stffvar_file=runParameters['LV geometrical overall stiffness constant variation file']
            fenics.ALE.move(refHeart.get_Mesh(),newDeform)
            self.savehdf5(refHeart,casename,meshname,'tempadjmesh_'+meshname,saveFolder=saveFolder)
            heartModelTranslator.changeFiberAngles(self.heartModelStr,runParameters+{'path of folder with case':casename,'name of mesh':'tempadjmesh_'+meshname,'subfolder to save files':saveFolder,'subfolder to load files':saveFolder,'LV geometrical fiber direction variation file':temp_fiber_vec_file,'LV geometrical fiber delay variation file':temp_fiber_delay_file,'LV geometrical fiber strength variation file':temp_fiber_strength_file,'LV geometrical overall stiffness constant variation file':temp_fiber_stffvar_file})
            prevDeformVectors=deformVectors.copy()
            heart=heartModelTranslator.HeartModel(self.heartModelStr,runParameters+{'path of folder with case':casename+saveFolder,'name of mesh':'tempadjmesh_'+meshname})
            logging.info('newheart starting volume ='+repr(heart.get_CavityVolume()))
            self.solvePressure(heart,targetPressure,voladj=runParameters['solver volume adjustment'],maxVolumeBreak=newheartVol*softAdjust)
            deformVectors=heart.get_DisplacementResult().vector()[:].copy()
            maxcoordchange=np.sum((heart.get_MeshCoordinates()-prevMeshCoords)**2.,axis=1).max()**1.5
            prevMeshCoords=heart.get_MeshCoordinates().copy()
            count+=1
            if (maxcoordchange/min(targetVolume))<10**-8:
                break
            else:
                logging.info('max coord change ratio= '+repr(maxcoordchange/targetVolume))
            if count>iterationNumber:
                logging.warning('Maximum Iteration reached for adjustMeshToPressure')
                break
        logging.info('Final heart volume ='+repr(heart.get_CavityVolume()))
        if savename is not None:
            self.savehdf5(heart,casename,meshname,savename,saveFolder=saveFolder)
            heartModelTranslator.changeFiberAngles(self.heartModelStr,runParameters+{'path of folder with case':casename,'name of mesh':savename,'subfolder to save files':saveFolder,'subfolder to load files':saveFolder,'LV geometrical fiber direction variation file':temp_fiber_vec_file,'LV geometrical fiber delay variation file':temp_fiber_delay_file,'LV geometrical fiber strength variation file':temp_fiber_strength_file,'LV geometrical overall stiffness constant variation file':temp_fiber_stffvar_file})
        return heart
    def getHeartCoords(self,casename,meshname,time):
        timestep=self.timetotimestepScale*time+self.timetotimestepShift
        heart=heartModelTranslator.HeartModel(self.heartModelStr,self.defaultParameters+{'path of folder with case':casename,'name of mesh':meshname})
        coords=np.copy(heart.get_MeshCoordinates())
        newcoords=np.zeros((coords.shape[0],4))
        newcoords[:,:3]=coords
        newcoords[:,3]=1.
        transform_mat=np.loadtxt(casename+'/'+meshname+'_Tmatrix_fromFEniCStoSTL.txt')
        coords=transform_mat.dot(newcoords.T)[:3,:].T
        coords=self.evaluateSTLCoords(coords,timestep)
        transform_mat=np.loadtxt(casename+'/'+meshname+'_Tmatrix_fromSTLtoFEniCS.txt')
        newcoords[:,:3]=coords
        newcoords[:,3]=1.
        coords=transform_mat.dot(newcoords.T)[:3,:].T
        return coords
    def adjustMeshToPressureWithInverse(self,casename,meshname,targetPressure,savename=None,softAdjust=float('inf'),iterationNumber=float('inf'),initHeart=None,prevDeform=None,saveFolder='',tolerance=10.**-4.,unloadedVolume=None):
        #set usePrevious to False or 0 to not use previous mesh for deformation, set True to use previous mesh only and set a number to use previous mesh a number of times is it doesnt converge
        runParameters=self.defaultParameters.duplicate()
        runParameters['maximum LV fiber tension in Pa']=0.0
        if isinstance(unloadedVolume,dict):
            if self.heartModelStr=='heArt.BiV':
                chambers=['LV','RV']
            elif self.heartModelStr=='AVHeart':
                if runParameters['Mitral valve channel control'] == False:
                    chambers=['LV','LA']
                else:
                    chambers=['LV','LA','Mitral']
            else:
                chambers=['LV']
            fixed_volume=[]
            for chamber in chambers:
                try:
                    fixed_volume.append(unloadedVolume[chamber])
                except:
                    fixed_volume.append(0.)
            fixed_volume=np.array(fixed_volume)
        elif isinstance(unloadedVolume,np.ndarray):#in this case targetPressure is the stiffness
            fixed_volume=unloadedVolume.copy()
            runParameters["strain energy density function coefficient in Pa"]=targetPressure
        else:
            fixed_volume=np.array([0.]*len(targetPressure))
        refHeart=heartModelTranslator.HeartModel(self.heartModelStr,runParameters+{'path of folder with case':casename,'name of mesh':meshname})
        unloadresultWriter=fenicsResultWriter(self.casename+"/unloadsave",["displacement"])
        targetVolume=refHeart.get_CavityVolume()
        logging.info('Fixed Volume = '+repr(fixed_volume))
        if not(isinstance(targetPressure,float)):
            new_targetPressure=targetPressure.copy()
        else:
            new_targetPressure=np.array([targetPressure]*len(fixed_volume))
        new_targetPressure[fixed_volume!=0]=0.
        logging.info('Target Volume = '+repr(targetVolume)+', target pressure = '+repr(new_targetPressure))
        if np.all(fixed_volume!=0.):
            targetMeshCoords=refHeart.get_MeshCoordinates().copy()
            heart=heartModelTranslator.HeartModel(self.heartModelStr,runParameters+{'path of folder with case':casename,'name of mesh':meshname})
            self.solveVolume(heart,unloadedVolume,voladj=runParameters['solver volume adjustment'])
        else:
            targetMeshCoords=refHeart.get_MeshCoordinates().copy()
            self.solvePressure(refHeart,new_targetPressure,voladj=runParameters['solver volume adjustment'],fixed_volume=fixed_volume)
            finv=fenics.grad(refHeart.get_DisplacementResult())+fenics.Identity(refHeart.get_DisplacementFunction().geometric_dimension())
            heart=heartModelTranslator.HeartModel(self.heartModelStr,runParameters+{'path of folder with case':casename,'name of mesh':meshname,'heart to solve for inverse':True})
            heart.get_InverseStretchTensorFunction().vector()[:]=fenics.project(finv,refHeart.get_InverseStretchTensorSpace()).vector()[:].copy()
    
            heart.solve()
        unloadresultWriter(heart)
        deform=fenics.project(heart.get_DisplacementResult(),heart.get_DisplacementSpace())
        deform_vectors=deform.vector()[:].copy()
        if isinstance(runParameters['LV geometrical fiber direction variation file'],str):
            temp_fiber_vec_file=runParameters['LV geometrical fiber direction variation file'][:-4]+'_tempadjmesh'+runParameters['LV geometrical fiber direction variation file'][-4:]
            fiber_temp=heart.get_StretchTensorFunction()
            temp_u=heart.get_DisplacementFunction()
            old_fiber=np.loadtxt(runParameters['LV geometrical fiber direction variation file'])
            temp_probes_f=np.empty((old_fiber.shape[0],9),dtype=float)
            temp_probes_u=np.empty((old_fiber.shape[0],3),dtype=float)
            fiber_temp=fenics.project(fiber_temp,fenics.TensorFunctionSpace(heart.get_Mesh(),'CG',2))
            fiber_temp.set_allow_extrapolation(True)
            temp_u=fenics.project(temp_u,fenics.VectorFunctionSpace(heart.get_Mesh(),'CG',2))
            temp_u.set_allow_extrapolation(True)
            for n in range(old_fiber.shape[0]):
                fiber_temp.eval(temp_probes_f[n],old_fiber[n,:3])
                temp_u.eval(temp_probes_u[n],old_fiber[n,:3])
            temp_probes_f=temp_probes_f.reshape((-1,3,3))
            temp_probes_f=np.matmul(temp_probes_f,old_fiber[:,:3].reshape((-1,3,1))).reshape((-1,3))
            temp_probes_f=temp_probes_f/np.linalg.norm(temp_probes_f,axis=-1,keepdims=True)
            old_fiber[:,:3]=old_fiber[:,:3]+temp_probes_u
            old_fiber[:,3:]=temp_probes_f[:]
            np.savetxt(temp_fiber_vec_file,old_fiber)
            store_temp_fiber_vec_file=old_fiber.copy()
        else:
            temp_fiber_vec_file=runParameters['LV geometrical fiber direction variation file']
        if isinstance(runParameters['LV geometrical fiber delay variation file'],str):
            temp_u=heart.get_DisplacementFunction()
            old_fiber=np.loadtxt(runParameters['LV geometrical fiber delay variation file'])
            temp_probes_u=np.empty((old_fiber.shape[0],3),dtype=float)
            temp_u=fenics.project(temp_u,fenics.VectorFunctionSpace(heart.get_Mesh(),'CG',2))
            temp_u.set_allow_extrapolation(True)
            for n in range(old_fiber.shape[0]):
                temp_u.eval(temp_probes_u[n],old_fiber[n,:3])
            old_fiber[:,:3]=old_fiber[:,:3]+temp_probes_u
            temp_fiber_delay_file=runParameters['LV geometrical fiber delay variation file'][:-4]+'_tempadjmesh'+runParameters['LV geometrical fiber delay variation file'][-4:]
            np.savetxt(temp_fiber_delay_file,old_fiber)
            store_temp_fiber_delay_file=old_fiber.copy()
        else:
            temp_fiber_delay_file=runParameters['LV geometrical fiber delay variation file']
        if isinstance(runParameters['LV geometrical fiber strength variation file'],str):
            temp_u=heart.get_DisplacementFunction()
            old_fiber=np.loadtxt(runParameters['LV geometrical fiber strength variation file'])
            temp_probes_u=np.empty((old_fiber.shape[0],3),dtype=float)
            temp_u=fenics.project(temp_u,fenics.VectorFunctionSpace(heart.get_Mesh(),'CG',2))
            temp_u.set_allow_extrapolation(True)
            for n in range(old_fiber.shape[0]):
                temp_u.eval(temp_probes_u[n],old_fiber[n,:3])
            old_fiber[:,:3]=old_fiber[:,:3]+temp_probes_u
            temp_fiber_strength_file=runParameters['LV geometrical fiber strength variation file'][:-4]+'_tempadjmesh'+runParameters['LV geometrical fiber strength variation file'][-4:]
            np.savetxt(temp_fiber_strength_file,old_fiber)
            store_temp_fiber_strength_file=old_fiber.copy()
        else:
            temp_fiber_strength_file=runParameters['LV geometrical fiber strength variation file']
        if isinstance(runParameters['LV geometrical overall stiffness constant variation file'],str):
            temp_u=heart.get_DisplacementFunction()
            old_fiber=np.loadtxt(runParameters['LV geometrical overall stiffness constant variation file'])
            temp_probes_u=np.empty((old_fiber.shape[0],3),dtype=float)
            temp_u=fenics.project(temp_u,fenics.VectorFunctionSpace(heart.get_Mesh(),'CG',2))
            temp_u.set_allow_extrapolation(True)
            for n in range(old_fiber.shape[0]):
                temp_u.eval(temp_probes_u[n],old_fiber[n,:3])
            old_fiber[:,:3]=old_fiber[:,:3]+temp_probes_u
            temp_fiber_stffvar_file=runParameters['LV geometrical overall stiffness constant variation file'][:-4]+'_tempadjmesh'+runParameters['LV geometrical overall stiffness constant variation file'][-4:]
            np.savetxt(temp_fiber_stffvar_file,old_fiber)
            store_temp_fiber_stffvar_file=old_fiber.copy()
        else:
            temp_fiber_stffvar_file=runParameters['LV geometrical overall stiffness constant variation file']
        fenics.ALE.move(heart.get_Mesh(),deform)
        self.savehdf5(heart,casename,meshname,'tempadjmesh_'+meshname,saveFolder=saveFolder)
        heartModelTranslator.changeFiberAngles(self.heartModelStr,runParameters+{'path of folder with case':casename,'name of mesh':'tempadjmesh_'+meshname,'subfolder to save files':saveFolder,'subfolder to load files':saveFolder,'LV geometrical fiber direction variation file':temp_fiber_vec_file,'LV geometrical fiber delay variation file':temp_fiber_delay_file,'LV geometrical fiber strength variation file':temp_fiber_strength_file,'LV geometrical overall stiffness constant variation file':temp_fiber_stffvar_file})
        if not(np.all(fixed_volume!=0.)) and np.any(fixed_volume==0.):
            heart=heartModelTranslator.HeartModel(self.heartModelStr,runParameters+{'path of folder with case':casename+saveFolder,'name of mesh':'tempadjmesh_'+meshname})
            temp_target_volume=heart.get_CavityVolume().copy()
            temp_target_volume[fixed_volume!=0.]=fixed_volume[fixed_volume!=0.]
            self.solveVolume(heart,temp_target_volume,voladj=runParameters['solver volume adjustment'])

            deform=fenics.project(heart.get_DisplacementResult(),heart.get_DisplacementSpace())
            deform_vectors=deform.vector()[:].copy()
            if isinstance(runParameters['LV geometrical fiber direction variation file'],str):
                temp_fiber_vec_file=runParameters['LV geometrical fiber direction variation file'][:-4]+'_tempadjmesh'+runParameters['LV geometrical fiber direction variation file'][-4:]
                fiber_temp=heart.get_StretchTensorFunction()
                temp_u=heart.get_DisplacementFunction()
                old_fiber=np.loadtxt(temp_fiber_vec_file)
                temp_probes_f=np.empty((old_fiber.shape[0],9),dtype=float)
                temp_probes_u=np.empty((old_fiber.shape[0],3),dtype=float)
                fiber_temp=fenics.project(fiber_temp,fenics.TensorFunctionSpace(heart.get_Mesh(),'CG',2))
                fiber_temp.set_allow_extrapolation(True)
                temp_u=fenics.project(temp_u,fenics.VectorFunctionSpace(heart.get_Mesh(),'CG',2))
                temp_u.set_allow_extrapolation(True)
                for n in range(old_fiber.shape[0]):
                    fiber_temp.eval(temp_probes_f[n],old_fiber[n,:3])
                    temp_u.eval(temp_probes_u[n],old_fiber[n,:3])
                temp_probes_f=temp_probes_f.reshape((-1,3,3))
                temp_probes_f=np.matmul(temp_probes_f,old_fiber[:,:3].reshape((-1,3,1))).reshape((-1,3))
                temp_probes_f=temp_probes_f/np.linalg.norm(temp_probes_f,axis=-1,keepdims=True)
                old_fiber[:,:3]=old_fiber[:,:3]+temp_probes_u
                old_fiber[:,3:]=temp_probes_f[:]
                np.savetxt(temp_fiber_vec_file,old_fiber)
            else:
                temp_fiber_vec_file=runParameters['LV geometrical fiber direction variation file']
            if isinstance(runParameters['LV geometrical fiber delay variation file'],str):
                temp_u=heart.get_DisplacementFunction()
                temp_fiber_delay_file=runParameters['LV geometrical fiber delay variation file'][:-4]+'_tempadjmesh'+runParameters['LV geometrical fiber delay variation file'][-4:]
                old_fiber=np.loadtxt(temp_fiber_delay_file)
                temp_probes_u=np.empty((old_fiber.shape[0],3),dtype=float)
                temp_u=fenics.project(temp_u,fenics.VectorFunctionSpace(heart.get_Mesh(),'CG',2))
                temp_u.set_allow_extrapolation(True)
                for n in range(old_fiber.shape[0]):
                    temp_u.eval(temp_probes_u[n],old_fiber[n,:3])
                old_fiber[:,:3]=old_fiber[:,:3]+temp_probes_u
                np.savetxt(temp_fiber_delay_file,old_fiber)
                store_temp_fiber_delay_file=old_fiber.copy()
            else:
                temp_fiber_delay_file=runParameters['LV geometrical fiber delay variation file']
            if isinstance(runParameters['LV geometrical fiber strength variation file'],str):
                temp_u=heart.get_DisplacementFunction()
                temp_fiber_strength_file=runParameters['LV geometrical fiber strength variation file'][:-4]+'_tempadjmesh'+runParameters['LV geometrical fiber strength variation file'][-4:]
                old_fiber=np.loadtxt(temp_fiber_strength_file)
                temp_probes_u=np.empty((old_fiber.shape[0],3),dtype=float)
                temp_u=fenics.project(temp_u,fenics.VectorFunctionSpace(heart.get_Mesh(),'CG',2))
                temp_u.set_allow_extrapolation(True)
                for n in range(old_fiber.shape[0]):
                    temp_u.eval(temp_probes_u[n],old_fiber[n,:3])
                old_fiber[:,:3]=old_fiber[:,:3]+temp_probes_u
                np.savetxt(temp_fiber_strength_file,old_fiber)
                store_temp_fiber_strength_file=old_fiber.copy()
            else:
                temp_fiber_strength_file=runParameters['LV geometrical fiber strength variation file']
            if isinstance(runParameters['LV geometrical overall stiffness constant variation file'],str):
                temp_u=heart.get_DisplacementFunction()
                temp_fiber_stffvar_file=runParameters['LV geometrical overall stiffness constant variation file'][:-4]+'_tempadjmesh'+runParameters['LV geometrical overall stiffness constant variation file'][-4:]
                old_fiber=np.loadtxt(temp_fiber_stffvar_file)
                temp_probes_u=np.empty((old_fiber.shape[0],3),dtype=float)
                temp_u=fenics.project(temp_u,fenics.VectorFunctionSpace(heart.get_Mesh(),'CG',2))
                temp_u.set_allow_extrapolation(True)
                for n in range(old_fiber.shape[0]):
                    temp_u.eval(temp_probes_u[n],old_fiber[n,:3])
                old_fiber[:,:3]=old_fiber[:,:3]+temp_probes_u
                np.savetxt(temp_fiber_stffvar_file,old_fiber)
                store_temp_fiber_stffvar_file=old_fiber.copy()
            else:
                temp_fiber_stffvar_file=runParameters['LV geometrical overall stiffness constant variation file']
            fenics.ALE.move(heart.get_Mesh(),deform)
            self.savehdf5(heart,casename,meshname,'tempadjmesh_'+meshname,saveFolder=saveFolder)
            heartModelTranslator.changeFiberAngles(self.heartModelStr,runParameters+{'path of folder with case':casename,'name of mesh':'tempadjmesh_'+meshname,'subfolder to save files':saveFolder,'subfolder to load files':saveFolder,'LV geometrical fiber direction variation file':temp_fiber_vec_file,'LV geometrical fiber delay variation file':temp_fiber_delay_file,'LV geometrical fiber strength variation file':temp_fiber_strength_file,'LV geometrical overall stiffness constant variation file':temp_fiber_stffvar_file})
            
        maxadjustmentsq=None
        newheartVol=heart.get_CavityVolume()
        prevHeartVol=heart.get_CavityVolume()
        
        refHeart=heartModelTranslator.HeartModel(self.heartModelStr,runParameters+{'path of folder with case':casename+saveFolder,'name of mesh':'tempadjmesh_'+meshname})
        logging.info('New start Volume = '+repr(refHeart.get_CavityVolume()))
        prevMeshCoords=refHeart.get_MeshCoordinates().copy()
        if np.all(fixed_volume!=0.):
            self.solveVolume(refHeart,targetVolume,voladj=runParameters['solver volume adjustment'])
        else:
            self.solvePressure(refHeart,new_targetPressure,voladj=runParameters['solver volume adjustment'],fixed_volume=fixed_volume)
        meanSqError=pointsMeansqDiff_kabsch(refHeart.get_MeshCoordinates(),targetMeshCoords)
        logging.info('meanSqError= '+repr(meanSqError))
        count=0
        while np.any(np.abs(refHeart.get_CavityPressure()[fixed_volume==0] /new_targetPressure[fixed_volume==0]-1)>tolerance) or np.any(np.abs(refHeart.get_CavityVolume()/targetVolume-1)>tolerance):##alternate convergence criteria##(np.mean((heart.get_MeshCoordinates()-du.vector()[:].reshape((-1,3))-refHeart.get_MeshCoordinates())**2.).max()**1.5/refHeart.get_CavityVolume())<10.**-6.:
            if count>=iterationNumber:
                logging.warning('Maximum Iteration reached for adjustMeshToPressure')
                break
            if np.all(fixed_volume!=0.):
                break
            new_deform_vectors=fenics.project(refHeart.get_DisplacementResult(),refHeart.displacementSpace).vector()[:]
            originalHeart=heartModelTranslator.HeartModel(self.heartModelStr,runParameters+{'path of folder with case':casename,'name of mesh':meshname})
            inv_deform=fenics.project(originalHeart.get_DisplacementResult(),originalHeart.get_DisplacementSpace())
            inv_deform.vector()[:]=deform_vectors+new_deform_vectors
            finv=fenics.grad(inv_deform)+fenics.Identity(originalHeart.get_DisplacementFunction().geometric_dimension())
            
            
            heart=heartModelTranslator.HeartModel(self.heartModelStr,runParameters+{'path of folder with case':casename+saveFolder,'name of mesh':'tempadjmesh_'+meshname,'heart to solve for inverse':True})
            heart.get_InverseStretchTensorFunction().vector()[:]=fenics.project(finv,originalHeart.get_InverseStretchTensorSpace()).vector()[:].copy()
            heart.solve()
            unloadresultWriter(heart)
            deform=fenics.project(heart.get_DisplacementResult(),heart.get_DisplacementSpace())
            deform_vectors=deform_vectors+deform.vector()[:].copy()
            
            deform=fenics.project(originalHeart.get_DisplacementResult(),originalHeart.get_DisplacementSpace())
            deform.vector()[:]=deform_vectors.copy()
            if isinstance(runParameters['LV geometrical fiber direction variation file'],str):
                temp_fiber_vec_file=runParameters['LV geometrical fiber direction variation file'][:-4]+'_tempadjmesh'+runParameters['LV geometrical fiber direction variation file'][-4:]
                fiber_temp=heart.get_StretchTensorFunction()
                temp_u=heart.get_DisplacementFunction()
                old_fiber=np.loadtxt(runParameters['LV geometrical fiber direction variation file'])
                temp_probes_f=np.empty((old_fiber.shape[0],9),dtype=float)
                temp_probes_u=np.empty((old_fiber.shape[0],3),dtype=float)
                fiber_temp=fenics.project(fiber_temp,fenics.TensorFunctionSpace(heart.get_Mesh(),'CG',2))
                fiber_temp.set_allow_extrapolation(True)
                temp_u=fenics.project(temp_u,fenics.VectorFunctionSpace(heart.get_Mesh(),'CG',2))
                temp_u.set_allow_extrapolation(True)
                for n in range(old_fiber.shape[0]):
                    fiber_temp.eval(temp_probes_f[n],old_fiber[n,:3])
                    temp_u.eval(temp_probes_u[n],old_fiber[n,:3])
                temp_probes_f=temp_probes_f.reshape((-1,3,3))
                temp_probes_f=np.matmul(temp_probes_f,old_fiber[:,:3].reshape((-1,3,1))).reshape((-1,3))
                temp_probes_f=temp_probes_f/np.linalg.norm(temp_probes_f,axis=-1,keepdims=True)
                old_fiber[:,:3]=old_fiber[:,:3]+temp_probes_u
                old_fiber[:,3:]=temp_probes_f[:]
                np.savetxt(temp_fiber_vec_file,old_fiber)
            else:
                temp_fiber_vec_file=runParameters['LV geometrical fiber direction variation file']
            if isinstance(runParameters['LV geometrical fiber delay variation file'],str):
                temp_u=heart.get_DisplacementFunction()
                old_fiber=np.loadtxt(runParameters['LV geometrical fiber delay variation file'])
                temp_probes_u=np.empty((old_fiber.shape[0],3),dtype=float)
                temp_u=fenics.project(temp_u,fenics.VectorFunctionSpace(heart.get_Mesh(),'CG',2))
                temp_u.set_allow_extrapolation(True)
                for n in range(old_fiber.shape[0]):
                    temp_u.eval(temp_probes_u[n],old_fiber[n,:3])
                old_fiber[:,:3]=old_fiber[:,:3]+temp_probes_u
                temp_fiber_delay_file=runParameters['LV geometrical fiber delay variation file'][:-4]+'_tempadjmesh'+runParameters['LV geometrical fiber delay variation file'][-4:]
                np.savetxt(temp_fiber_delay_file,old_fiber)
            else:
                temp_fiber_delay_file=runParameters['LV geometrical fiber delay variation file']
            if isinstance(runParameters['LV geometrical fiber strength variation file'],str):
                temp_u=heart.get_DisplacementFunction()
                old_fiber=np.loadtxt(runParameters['LV geometrical fiber strength variation file'])
                temp_probes_u=np.empty((old_fiber.shape[0],3),dtype=float)
                temp_u=fenics.project(temp_u,fenics.VectorFunctionSpace(heart.get_Mesh(),'CG',2))
                temp_u.set_allow_extrapolation(True)
                for n in range(old_fiber.shape[0]):
                    temp_u.eval(temp_probes_u[n],old_fiber[n,:3])
                old_fiber[:,:3]=old_fiber[:,:3]+temp_probes_u
                temp_fiber_strength_file=runParameters['LV geometrical fiber strength variation file'][:-4]+'_tempadjmesh'+runParameters['LV geometrical fiber strength variation file'][-4:]
                np.savetxt(temp_fiber_strength_file,old_fiber)
            else:
                temp_fiber_strength_file=runParameters['LV geometrical fiber strength variation file']
            if isinstance(runParameters['LV geometrical overall stiffness constant variation file'],str):
                temp_u=heart.get_DisplacementFunction()
                old_fiber=np.loadtxt(runParameters['LV geometrical overall stiffness constant variation file'])
                temp_probes_u=np.empty((old_fiber.shape[0],3),dtype=float)
                temp_u=fenics.project(temp_u,fenics.VectorFunctionSpace(heart.get_Mesh(),'CG',2))
                temp_u.set_allow_extrapolation(True)
                for n in range(old_fiber.shape[0]):
                    temp_u.eval(temp_probes_u[n],old_fiber[n,:3])
                old_fiber[:,:3]=old_fiber[:,:3]+temp_probes_u
                temp_fiber_stffvar_file=runParameters['LV geometrical overall stiffness constant variation file'][:-4]+'_tempadjmesh'+runParameters['LV geometrical overall stiffness constant variation file'][-4:]
                np.savetxt(temp_fiber_stffvar_file,old_fiber)
            else:
                temp_fiber_stffvar_file=runParameters['LV geometrical overall stiffness constant variation file']
            fenics.ALE.move(originalHeart.get_Mesh(),deform)
            self.savehdf5(originalHeart,casename,meshname,'tempadjmesh_'+meshname,saveFolder=saveFolder)
            heartModelTranslator.changeFiberAngles(self.heartModelStr,runParameters+{'path of folder with case':casename,'name of mesh':'tempadjmesh_'+meshname,'subfolder to save files':saveFolder,'subfolder to load files':saveFolder,'LV geometrical fiber direction variation file':temp_fiber_vec_file,'LV geometrical fiber delay variation file':temp_fiber_delay_file,'LV geometrical fiber strength variation file':temp_fiber_strength_file,'LV geometrical overall stiffness constant variation file':temp_fiber_stffvar_file})
            if not(np.all(fixed_volume!=0.)) and np.any(fixed_volume==0.):
                heart=heartModelTranslator.HeartModel(self.heartModelStr,runParameters+{'path of folder with case':casename+saveFolder,'name of mesh':'tempadjmesh_'+meshname})
                temp_target_volume=heart.get_CavityVolume().copy()
                temp_target_volume[fixed_volume!=0.]=fixed_volume[fixed_volume!=0.]
                self.solveVolume(heart,temp_target_volume,voladj=runParameters['solver volume adjustment'])

                deform=fenics.project(heart.get_DisplacementResult(),heart.get_DisplacementSpace())
                deform_vectors=deform.vector()[:].copy()
                if isinstance(runParameters['LV geometrical fiber direction variation file'],str):
                    temp_fiber_vec_file=runParameters['LV geometrical fiber direction variation file'][:-4]+'_tempadjmesh'+runParameters['LV geometrical fiber direction variation file'][-4:]
                    fiber_temp=heart.get_StretchTensorFunction()
                    temp_u=heart.get_DisplacementFunction()
                    old_fiber=np.loadtxt(temp_fiber_vec_file)
                    temp_probes_f=np.empty((old_fiber.shape[0],9),dtype=float)
                    temp_probes_u=np.empty((old_fiber.shape[0],3),dtype=float)
                    fiber_temp=fenics.project(fiber_temp,fenics.TensorFunctionSpace(heart.get_Mesh(),'CG',2))
                    fiber_temp.set_allow_extrapolation(True)
                    temp_u=fenics.project(temp_u,fenics.VectorFunctionSpace(heart.get_Mesh(),'CG',2))
                    temp_u.set_allow_extrapolation(True)
                    for n in range(old_fiber.shape[0]):
                        fiber_temp.eval(temp_probes_f[n],old_fiber[n,:3])
                        temp_u.eval(temp_probes_u[n],old_fiber[n,:3])
                    temp_probes_f=temp_probes_f.reshape((-1,3,3))
                    temp_probes_f=np.matmul(temp_probes_f,old_fiber[:,:3].reshape((-1,3,1))).reshape((-1,3))
                    temp_probes_f=temp_probes_f/np.linalg.norm(temp_probes_f,axis=-1,keepdims=True)
                    old_fiber[:,:3]=old_fiber[:,:3]+temp_probes_u
                    old_fiber[:,3:]=temp_probes_f[:]
                    np.savetxt(temp_fiber_vec_file,old_fiber)
                else:
                    temp_fiber_vec_file=runParameters['LV geometrical fiber direction variation file']
                if isinstance(runParameters['LV geometrical fiber delay variation file'],str):
                    temp_u=heart.get_DisplacementFunction()
                    temp_fiber_delay_file=runParameters['LV geometrical fiber delay variation file'][:-4]+'_tempadjmesh'+runParameters['LV geometrical fiber delay variation file'][-4:]
                    old_fiber=np.loadtxt(temp_fiber_delay_file)
                    temp_probes_u=np.empty((old_fiber.shape[0],3),dtype=float)
                    temp_u=fenics.project(temp_u,fenics.VectorFunctionSpace(heart.get_Mesh(),'CG',2))
                    temp_u.set_allow_extrapolation(True)
                    for n in range(old_fiber.shape[0]):
                        temp_u.eval(temp_probes_u[n],old_fiber[n,:3])
                    old_fiber[:,:3]=old_fiber[:,:3]+temp_probes_u
                    np.savetxt(temp_fiber_delay_file,old_fiber)
                    store_temp_fiber_delay_file=old_fiber.copy()
                else:
                    temp_fiber_delay_file=runParameters['LV geometrical fiber delay variation file']
                if isinstance(runParameters['LV geometrical fiber strength variation file'],str):
                    temp_u=heart.get_DisplacementFunction()
                    temp_fiber_strength_file=runParameters['LV geometrical fiber strength variation file'][:-4]+'_tempadjmesh'+runParameters['LV geometrical fiber strength variation file'][-4:]
                    old_fiber=np.loadtxt(temp_fiber_strength_file)
                    temp_probes_u=np.empty((old_fiber.shape[0],3),dtype=float)
                    temp_u=fenics.project(temp_u,fenics.VectorFunctionSpace(heart.get_Mesh(),'CG',2))
                    temp_u.set_allow_extrapolation(True)
                    for n in range(old_fiber.shape[0]):
                        temp_u.eval(temp_probes_u[n],old_fiber[n,:3])
                    old_fiber[:,:3]=old_fiber[:,:3]+temp_probes_u
                    np.savetxt(temp_fiber_strength_file,old_fiber)
                    store_temp_fiber_strength_file=old_fiber.copy()
                else:
                    temp_fiber_strength_file=runParameters['LV geometrical fiber strength variation file']
                if isinstance(runParameters['LV geometrical overall stiffness constant variation file'],str):
                    temp_u=heart.get_DisplacementFunction()
                    temp_fiber_stffvar_file=runParameters['LV geometrical overall stiffness constant variation file'][:-4]+'_tempadjmesh'+runParameters['LV geometrical overall stiffness constant variation file'][-4:]
                    old_fiber=np.loadtxt(temp_fiber_stffvar_file)
                    temp_probes_u=np.empty((old_fiber.shape[0],3),dtype=float)
                    temp_u=fenics.project(temp_u,fenics.VectorFunctionSpace(heart.get_Mesh(),'CG',2))
                    temp_u.set_allow_extrapolation(True)
                    for n in range(old_fiber.shape[0]):
                        temp_u.eval(temp_probes_u[n],old_fiber[n,:3])
                    old_fiber[:,:3]=old_fiber[:,:3]+temp_probes_u
                    np.savetxt(temp_fiber_stffvar_file,old_fiber)
                    store_temp_fiber_stffvar_file=old_fiber.copy()
                else:
                    temp_fiber_stffvar_file=runParameters['LV geometrical overall stiffness constant variation file']
                fenics.ALE.move(heart.get_Mesh(),deform)
                self.savehdf5(heart,casename,meshname,'tempadjmesh_'+meshname,saveFolder=saveFolder)
                heartModelTranslator.changeFiberAngles(self.heartModelStr,runParameters+{'path of folder with case':casename,'name of mesh':'tempadjmesh_'+meshname,'subfolder to save files':saveFolder,'subfolder to load files':saveFolder,'LV geometrical fiber direction variation file':temp_fiber_vec_file,'LV geometrical fiber delay variation file':temp_fiber_delay_file,'LV geometrical fiber strength variation file':temp_fiber_strength_file,'LV geometrical overall stiffness constant variation file':temp_fiber_stffvar_file})
               
            refHeart=heartModelTranslator.HeartModel(self.heartModelStr,runParameters+{'path of folder with case':casename+saveFolder,'name of mesh':'tempadjmesh_'+meshname})
            logging.info('New start Volume = '+repr(refHeart.get_CavityVolume()))
            self.solvePressure(refHeart,new_targetPressure,voladj=runParameters['solver volume adjustment'],fixed_volume=fixed_volume)
            meanSqError=pointsMeansqDiff_kabsch(refHeart.get_MeshCoordinates(),targetMeshCoords)
            logging.info('meanSqError= '+repr(meanSqError))
            maxcoordchange=pointsMeansqDiff_kabsch(refHeart.get_MeshCoordinates(),prevMeshCoords,averageFunc=np.max)**1.5
            prevMeshCoords=refHeart.get_MeshCoordinates().copy()
            count+=1
            if (maxcoordchange/targetVolume)<10**-8:
                break
            else:
                logging.info('max coord change ratio= '+repr(maxcoordchange/targetVolume))
        logging.info('Final heart volume ='+repr(refHeart.get_CavityVolume()))
        logging.info('Final heart pressure ='+repr(refHeart.get_CavityPressure()))
        if savename is not None:
            self.savehdf5(refHeart,casename,meshname,savename,saveFolder=saveFolder)
            heartModelTranslator.changeFiberAngles(self.heartModelStr,runParameters+{'path of folder with case':casename,'name of mesh':savename,'subfolder to save files':saveFolder,'subfolder to load files':saveFolder,'LV geometrical fiber direction variation file':temp_fiber_vec_file,'LV geometrical fiber delay variation file':temp_fiber_delay_file,'LV geometrical fiber strength variation file':temp_fiber_strength_file,'LV geometrical overall stiffness constant variation file':temp_fiber_stffvar_file})
        return refHeart
    def getUnloadedGeometry(self,editParameters=None,casename=None,meshname=None,targetPressure=None,targetVolume=None,tryPressureRatio=0.2,targetMeshCoords=None,savename=None,iterationNumber=0,toRunCountFolder=False,inverseHeart=True,unloadedVolume=None,unload_loop_num=float('inf')):
        #when targetMeshCoords is not None, set target pressure as the maximum pressure to try
        #= mesh.coordinates()[:] at target volume
            
        runParameters=self.defaultParameters.duplicate()
        runParameters['maximum LV fiber tension in Pa']=0.0
        if editParameters is not None:
            for key in runParameters:
                if key in editParameters:
                    runParameters[key]=editParameters[key]
        else:
            editParameters={}
        if unloadedVolume is not None:
            logger.warning('Setting unloadedVolume switches process to optimising "strain energy density function coefficient in Pa".')
            inputStiffness=runParameters["strain energy density function coefficient in Pa"]
        if toRunCountFolder:
            if isinstance(toRunCountFolder,str):
                if toRunCountFolder[0]!='/':
                    toRunCountFolder="/"+toRunCountFolder
            else:
                toRunCountFolder="/"+str(self.runCount)
        else:
            toRunCountFolder=''
        if casename is None:
            casename=self.casename
        if meshname is None:
            meshname=self.meshname
        if isinstance(runParameters['LV geometrical fiber direction variation file'],str):
            temp_fiber_vec_file=runParameters['LV geometrical fiber direction variation file'][:-4]+'_tempadjmesh'+runParameters['LV geometrical fiber direction variation file'][-4:]
        else:
            temp_fiber_vec_file=runParameters['LV geometrical fiber direction variation file']
        if isinstance(runParameters['LV geometrical fiber delay variation file'],str):
            temp_fiber_delay_file=runParameters['LV geometrical fiber delay variation file'][:-4]+'_tempadjmesh'+runParameters['LV geometrical fiber delay variation file'][-4:]
        else:
            temp_fiber_delay_file=runParameters['LV geometrical fiber delay variation file']
        if isinstance(runParameters['LV geometrical fiber strength variation file'],str):
            temp_fiber_strength_file=runParameters['LV geometrical fiber strength variation file'][:-4]+'_tempadjmesh'+runParameters['LV geometrical fiber strength variation file'][-4:]
        else:
            temp_fiber_strength_file=runParameters['LV geometrical fiber strength variation file']
        if isinstance(runParameters['LV geometrical overall stiffness constant variation file'],str):
            temp_fiber_stffvar_file=runParameters['LV geometrical overall stiffness constant variation file'][:-4]+'_tempadjmesh'+runParameters['LV geometrical overall stiffness constant variation file'][-4:]
        else:
            temp_fiber_stffvar_file=runParameters['LV geometrical overall stiffness constant variation file']
        
        if targetPressure is None:
            if self.heartModelStr=='heArt.BiV':
                targetPressure=np.array([runParameters['LV end diastolic pressure in mmHg'],runParameters['RV end diastolic pressure in mmHg']])
            elif self.heartModelStr=='AVHeart':
                if runParameters['Mitral valve channel control'] == False:
                    targetPressure=np.array([runParameters['LV end diastolic pressure in mmHg'],runParameters['LA end diastolic pressure in mmHg']])
                else:
                    targetPressure=np.array([runParameters['LV end diastolic pressure in mmHg'],runParameters['LA end diastolic pressure in mmHg'],runParameters['Mitral end diastolic pressure in mmHg']])
            else:
                targetPressure=np.array([runParameters['LV end diastolic pressure in mmHg']])
        elif isinstance(targetPressure,(int,float)):
            if self.heartModelStr in ['heArt.BiV']:
                targetPressure=[targetPressure]*2
            elif self.heartModelStr in ['AVHeart']:
                if runParameters['Mitral valve channel control'] == False:
                    targetPressure=[targetPressure]*3
                else:
                    targetPressure=[targetPressure]*2
            else:
                targetPressure=[targetPressure]
        if targetPressure is not None:
            targetPressure=np.array(targetPressure)
        for n in range(targetPressure.shape[0]):
            if targetPressure[n] is None:
                targetPressure[n]=0.
        condition_boolean=targetPressure!=0
        if targetVolume is None:
            if self.heartModelStr=='heArt.BiV':
                targetVolume=np.array([runParameters['LV end diastolic volume in mL'],runParameters['RV end diastolic volume in mL']])
            elif self.heartModelStr=='AVHeart':
                if runParameters['Mitral valve channel control'] == False:
                    targetVolume=np.array([runParameters['LV end diastolic volume in mL'],runParameters['LA end diastolic volume in mL']])
                else:
                    targetVolume=np.array([runParameters['LV end diastolic volume in mL'],runParameters['LA end diastolic volume in mL'],runParameters['Mitral end diastolic volume in mL']])
            else:
                targetVolume=np.array([runParameters['LV end diastolic volume in mL']])
        elif isinstance(targetVolume,(int,float)):
            targetVolume=[targetVolume]*len(targetPressure)
        targetVolume=np.array(targetVolume)
        if unloadedVolume is not None:
            if isinstance(unloadedVolume,(int,float)):
                if self.heartModelStr=='heArt.BiV':
                    unloadedVolume=unloadedVolume*np.array([runParameters['LV end systolic volume in mL'],runParameters['RV end systolic volume in mL']])
                elif self.heartModelStr=='AVHeart':
                    if runParameters['Mitral valve channel control'] == False:
                        unloadedVolume=unloadedVolume*np.array([runParameters['LV end systolic volume in mL'],runParameters['LA end systolic volume in mL']])
                    else:
                        unloadedVolume=unloadedVolume*np.array([runParameters['LV end systolic volume in mL'],runParameters['LA end systolic volume in mL'],runParameters['Mitral end systolic volume in mL']])
                else:
                    unloadedVolume=unloadedVolume*np.array([runParameters['LV end diastolic volume in mL']])
            if isinstance(unloadedVolume,dict):
                try:
                    if self.heartModelStr=='heArt.BiV':
                        unloadedVolume=np.array([unloadedVolume['LV'],unloadedVolume['RV']])
                    elif self.heartModelStr=='AVHeart':
                        if runParameters['Mitral valve channel control'] == False:
                            unloadedVolume=np.array([unloadedVolume['LV'],unloadedVolume['LA']])
                        else:
                            unloadedVolume=np.array([unloadedVolume['LV'],unloadedVolume['LA'],unloadedVolume['Mitral']])
                    else:
                        unloadedVolume=np.array([unloadedVolume['LV']])
                except:
                    pass
        if unloadedVolume is not None and not(isinstance(unloadedVolume,dict)):
            tryPressureRatio=1.
        elif isinstance(tryPressureRatio,(int,float)):
            tryPressureRatio=np.array([tryPressureRatio]*len(targetPressure))
        heart=heartModelTranslator.HeartModel(self.heartModelStr,runParameters+{'path of folder with case':casename,'name of mesh':meshname})
        print('Start unload with volume =',heart.get_CavityVolume())
        if unloadedVolume is not None:
            if not(isinstance(unloadedVolume,dict)):
                tryPressure=tryPressureRatio*inputStiffness
            else:
                tryPressure=tryPressureRatio*targetPressure
            adjtryPressure=np.minimum(tryPressureRatio*0.33,np.minimum((1-tryPressureRatio)*0.33,0.1))
        elif targetMeshCoords is not None:
            self.solveVolume(heart,targetVolume,voladj=runParameters['solver volume adjustment'])
            fenics.ALE.move(heart.get_Mesh(),heart.get_DisplacementResult())
            meshRMS_error=np.array([pointsMeansqDiff_kabsch(heart.get_MeshCoordinates(),targetMeshCoords),float('inf')])
            #meshRMS_error=np.array([pointsMeansqDiff_kabsch(heart.get_MeshCoordinates(),targetMeshCoords)*0.5,float('inf')])
            meshRMS_pressurebound=np.array([0,targetPressure])
            tryPressure=targetPressure
            adjtryPressure=targetPressure
        else:
            tryPressure=tryPressureRatio*targetPressure
            adjtryPressure=np.minimum(tryPressureRatio*0.33,np.minimum((1-tryPressureRatio)*0.33,0.1))
        if inverseHeart or unloadedVolume is not None:
            heart=self.adjustMeshToPressureWithInverse(casename,meshname,tryPressure,savename=savename,iterationNumber=iterationNumber,saveFolder=toRunCountFolder,unloadedVolume=unloadedVolume)
        else:
            heart=self.adjustMeshToPressure(casename,meshname,tryPressure,savename=savename,iterationNumber=iterationNumber,saveFolder=toRunCountFolder)
        self.savehdf5(heart,casename,meshname,'tempadjmesh_'+meshname,saveFolder=toRunCountFolder)
        self.savehdf5(heart,casename,meshname,savename,saveFolder=toRunCountFolder)#backup save
        heartModelTranslator.changeFiberAngles(self.heartModelStr,runParameters+{'path of folder with case':casename,'name of mesh':'tempadjmesh_'+meshname,'subfolder to save files':toRunCountFolder,'subfolder to load files':toRunCountFolder,'LV geometrical fiber direction variation file':temp_fiber_vec_file,'LV geometrical fiber delay variation file':temp_fiber_delay_file,'LV geometrical fiber strength variation file':temp_fiber_strength_file,'LV geometrical overall stiffness constant variation file':temp_fiber_stffvar_file})
        
        heart=heartModelTranslator.HeartModel(self.heartModelStr,runParameters+{'path of folder with case':casename+toRunCountFolder,'name of mesh':'tempadjmesh_'+meshname})
        self.solveVolume(heart,targetVolume,voladj=runParameters['solver volume adjustment'])
        
        if unloadedVolume is not None and not(isinstance(unloadedVolume,dict)):
            if np.sum(heart.get_CavityPressure()[condition_boolean] /targetPressure[condition_boolean]-1)>0:
                tunePressure=-1
            elif np.sum(heart.get_CavityPressure()[condition_boolean] /targetPressure[condition_boolean]-1)<0:
                tunePressure=1
            else:
                tunePressure=0
        else:
            tunePressure=[]
            for n in range(len(targetPressure)):
                if condition_boolean[n] is False:
                    tunePressure.append(0)
                elif heart.get_CavityPressure()[n] <targetPressure[n]:
                    tunePressure.append(1)
                else:
                    tunePressure.append(-1)
            tunePressure=np.array(tunePressure)
        count=0
        while np.any(np.abs(heart.get_CavityPressure()[condition_boolean] /targetPressure[condition_boolean]-1)>10**-4.):
            if unloadedVolume is not None and not(isinstance(unloadedVolume,dict)):
                if abs(np.sum(heart.get_CavityPressure()[condition_boolean] /targetPressure[condition_boolean]-1))<10**-4.:
                    break
            if count>=unload_loop_num:
                break
            if unloadedVolume is not None and not(isinstance(unloadedVolume,dict)):
                tryPressureRatio=tryPressureRatio+tunePressure*adjtryPressure
                tryPressure=tryPressureRatio*inputStiffness
                logging.info('Trying Stiffness '+repr(tryPressure))
            elif targetMeshCoords is not None:
                fenics.ALE.move(heart.get_Mesh(),heart.get_DisplacementResult())
                temperror=pointsMeansqDiff_kabsch(heart.get_MeshCoordinates(),targetMeshCoords)
                replaceInd=np.argmax(meshRMS_error)
                meshRMS_error[replaceInd]=temperror
                meshRMS_pressurebound[replaceInd]=tryPressure
                replaceInd=np.argmax(meshRMS_error)
                tryPressure=np.array([meshRMS_pressurebound[replaceInd]*2./3.+meshRMS_pressurebound[1-replaceInd]/3.]*len(targetPressure))
                adjtryPressure=meshRMS_pressurebound[1]-meshRMS_pressurebound[0]
            else:
                tryPressureRatio=tryPressureRatio+tunePressure*adjtryPressure
                tryPressure=tryPressureRatio*targetPressure
                logging.info('Trying Pressure '+repr(tryPressure))
            if inverseHeart or unloadedVolume is not None:
                heart=self.adjustMeshToPressureWithInverse(casename,meshname,tryPressure,savename=savename,iterationNumber=iterationNumber,saveFolder=toRunCountFolder,unloadedVolume=unloadedVolume)
            else:
                heart=self.adjustMeshToPressure(casename,meshname,tryPressure,savename=savename,iterationNumber=iterationNumber,saveFolder=toRunCountFolder,tolerance=10.**-3)
            self.savehdf5(heart,casename,meshname,'tempadjmesh_'+meshname,saveFolder=toRunCountFolder)
            heartModelTranslator.changeFiberAngles(self.heartModelStr,runParameters+{'path of folder with case':casename,'name of mesh':'tempadjmesh_'+meshname,'subfolder to save files':toRunCountFolder,'subfolder to load files':toRunCountFolder,'LV geometrical fiber direction variation file':temp_fiber_vec_file,'LV geometrical fiber delay variation file':temp_fiber_delay_file,'LV geometrical fiber strength variation file':temp_fiber_strength_file,'LV geometrical overall stiffness constant variation file':temp_fiber_stffvar_file})
            heart=heartModelTranslator.HeartModel(self.heartModelStr,runParameters+{'path of folder with case':casename+toRunCountFolder,'name of mesh':'tempadjmesh_'+meshname})
            self.solveVolume(heart,targetVolume,voladj=runParameters['solver volume adjustment'])
            if unloadedVolume is not None and not(isinstance(unloadedVolume,dict)):
                if tunePressure*np.sum(heart.get_CavityPressure()[condition_boolean] /targetPressure[condition_boolean]-1)>0:
                    tunePressure=-1*tunePressure
                    adjtryPressure*=0.33
            elif targetMeshCoords is not None:
                logging.info('Trying Pressure '+repr(tryPressure)+' between '+repr(meshRMS_pressurebound))
                logging.info('Errors '+repr(meshRMS_error))
            else:
                for n in range(len(targetPressure)):
                    if condition_boolean[n] is False:
                        continue
                    elif (tunePressure[n]*heart.get_CavityPressure()[n] )>(tunePressure[n]*targetPressure[n]):
                        tunePressure[n]*=-1
                        adjtryPressure[n]*=0.33
                    #elif tunePressure[n]==-1 and (tryPressureRatio[n]<1.5*adjtryPressure[n]):
                    #    adjtryPressure[n]*=0.33
            if np.all(np.logical_and(tryPressure<10.**-4. , tunePressure==-1)):
                logging.warning("Trying Negative pressures.")#break
            if np.all(adjtryPressure<10**-3):
                break
            count+=1
        self.savehdf5(heart,casename,meshname,savename,saveFolder=toRunCountFolder)
        heartModelTranslator.changeFiberAngles(self.heartModelStr,runParameters+{'path of folder with case':casename,'name of mesh':savename,'subfolder to save files':toRunCountFolder,'subfolder to load files':toRunCountFolder,'LV geometrical fiber direction variation file':temp_fiber_vec_file,'LV geometrical fiber delay variation file':temp_fiber_delay_file,'LV geometrical fiber strength variation file':temp_fiber_strength_file,'LV geometrical overall stiffness constant variation file':temp_fiber_stffvar_file})
        return heart
    def writeRunParameters(self,folder,runParameters):
        runParameters.writeParameters(folder+"/runParameters.txt")
    def readRunParameters(self,folder):
        runParameters=heartParameters.heartParameters(defaultAge=self.LVage)
        runParameters.readParameters(folder+"/runParameters.txt")
        return runParameters
    def getLastcycleCircuitResults(self,BCL,path=None,filename=None):
        if filename is None:
            filename='circuit_results'
        if path is None:
            path=self.casename+"/"+str(self.runCount)
        ngspice_py.getLastcycleCircuitResults(BCL,path,filename)
        return
    def __call__(self,editParameters=None,**kwargs):
        if self.defaultRun_kwargs is None:
            run_kwargs={}
        else:
            run_kwargs=self.defaultRun_kwargs
        for key in kwargs:
            run_kwargs[key]=kwargs[key]
        if self.defaultRunMode=='iterativeRun':
            runParameters=self.iterativeRun(editParameters=editParameters,**run_kwargs)
        elif self.defaultRunMode=='LVbehaviorRun':
            runParameters=self.LVbehaviorRun(editParameters=editParameters,**run_kwargs)
        elif self.defaultRunMode=='ESPVRRun':
            runParameters=self.LVbehaviorRun(editParameters=editParameters,getEDPVR=True,**run_kwargs)
        elif self.defaultRunMode=='unloadedGeometryRun':
            runParameters=self.unloadedGeometryRun(editParameters=editParameters,**run_kwargs)
        elif self.defaultRunMode=='fullWindkesselRun':
            runParameters=self.fullWindkesselRun(editParameters=editParameters,**run_kwargs)
        elif self.defaultRunMode=='flowrateWindkesselRun':
            runParameters=self.flowrateWindkesselRun(editParameters=editParameters,**run_kwargs)
        return runParameters
    def getEDPVR(self,editParameters=None,meshname=None,meshFolder=None,volstep=50,extendVol=0.1,minESV=None,maxEDV=None,saveaddstr='',voladj=0.1,toRunCountFolder=False):
        runParameters=self.defaultParameters.duplicate()
        if editParameters is not None:
            for key in runParameters:
                if key in editParameters:
                    runParameters[key]=editParameters[key]
        else:
            editParameters={}
        self.getLVbehavior(runParameters=runParameters,meshname=meshname,meshFolder=meshFolder,volstep=volstep,extendVol=extendVol,minESV=minESV,maxEDV=maxEDV,saveaddstr=saveaddstr,voladj=voladj,starttime=0,endtime=0,toRunCountFolder=toRunCountFolder)
    def getESPVR(self,editParameters=None,meshname=None,meshFolder=None,volstep=50,extendVol=0.1,minESV=None,maxEDV=None,saveaddstr='',voladj=0.1,toRunCountFolder=False):
        runParameters=self.defaultParameters.duplicate()
        if editParameters is not None:
            for key in runParameters:
                if key in editParameters:
                    runParameters[key]=editParameters[key]
        else:
            editParameters={}
        self.getLVbehavior(runParameters=runParameters,meshname=meshname,meshFolder=meshFolder,volstep=volstep,extendVol=extendVol,minESV=minESV,maxEDV=maxEDV,saveaddstr=saveaddstr,voladj=voladj,starttime=float(int(runParameters['time to maximum LV fiber tension in ms'])),endtime=float(int(runParameters['time to maximum LV fiber tension in ms'])),toRunCountFolder=toRunCountFolder)
    def LVbehaviorMultiCPUTransfer(self,transferFrom,editParameters=None,unloadGeo=True,folderToLVbehavior=None,minESV=None,maxEDV=None,volstep=50):
        #if using displacement to get unloaded geometry, set unloadGeo='displacement'
        if isinstance(folderToLVbehavior,str):
            if folderToLVbehavior[0]!='/':
                folderToLVbehavior='/'+folderToLVbehavior
        runParameters=self.defaultParameters.duplicate()
        if editParameters is not None:
            for key in runParameters:
                if key in editParameters:
                    runParameters[key]=editParameters[key]
        else:
            editParameters={}
        if unloadGeo:
            meshname=self.meshname+'_unloadedmesh'
        else:
            meshname=self.meshname
        if self.heartModelStr=='heArt.BiV':
            cavityNames=['LV','RV']
        elif self.heartModelStr=='AVHeart':
            if runParameters['Mitral valve channel control'] == False:
                cavityNames=['LV','LA']
            else:
                cavityNames=['LV','LA','Mitral']
        else:
            cavityNames=['LV']
            
        os.umask(0)
        os.makedirs(self.casename+"/"+str(self.runCount)+"/multipcu_LVBehavior",mode=0o777,exist_ok=True)
        if folderToLVbehavior is None or min(self.defaultParameters.getParameterRelation(editParameters.keys())+[2])<=1:
            if min(self.defaultParameters.getParameterRelation(editParameters.keys())+[2])<1:
                multipcu_kwargs={'meshname':meshname,'meshFolder':str(self.runCount),'volstep':volstep,'minESV':minESV,'maxEDV':maxEDV,'toRunCountFolder':True}
            elif unloadGeo:
                multipcu_kwargs={'meshname':meshname,'meshFolder':str(self.runCount),'volstep':volstep,'minESV':minESV,'maxEDV':maxEDV,'toRunCountFolder':True}
            else:
                multipcu_kwargs={'meshname':meshname,'volstep':volstep,'minESV':minESV,'maxEDV':maxEDV,'toRunCountFolder':True}
        elif not(os.path.isfile(self.casename+folderToLVbehavior+'/'+meshname+"_Press_VolTime.txt")):
            if unloadGeo:
                multipcu_kwargs={'meshname':meshname,'meshFolder':folderToLVbehavior,'volstep':volstep,'minESV':minESV,'maxEDV':maxEDV,'toRunCountFolder':folderToLVbehavior}
            else:
                multipcu_kwargs={'meshname':meshname,'volstep':volstep,'minESV':minESV,'maxEDV':maxEDV,'toRunCountFolder':folderToLVbehavior}
                
        volumeSpace,timeSpace=self.getLVbehavior(**multipcu_kwargs,returnVolspaceOnly=True,saveFile=True)
        total_volume_to_calculate=volumeSpace.shape[1]**volumeSpace.shape[0]
        
        insert_path_for_heartFEM = sys.path
        calculated_volume=np.loadtxt(transferFrom+"/"+meshname+"_Press_volumeSpace.txt")
        calculated_volume_array=[]
        for curr_volN in range(calculated_volume.shape[1]**calculated_volume.shape[0]):
             urInd=np.unravel_index(curr_volN, [calculated_volume.shape[1]]*calculated_volume.shape[0],order='F')
             calculated_volume_array.append([])
             for m in range(len(urInd)):
                calculated_volume_array[-1].append(calculated_volume[m,urInd[m]])
        calculated_volume_array=np.array(calculated_volume_array)
        logger.info("Commence transfer of "+repr(total_volume_to_calculate)+" volumes from '"+transferFrom+"/multipcu_LVBehavior'"+" to '"+self.casename+"/"+str(self.runCount)+"/multipcu_LVBehavior'")
        os.makedirs(self.casename+"/"+str(self.runCount)+"/multipcu_LVBehavior",mode=0o777,exist_ok=True)
        for curr_volN in range(total_volume_to_calculate):
            urInd=np.unravel_index(curr_volN, [volumeSpace.shape[1]]*volumeSpace.shape[0],order='F')
            volume_toget=[]
            for m in range(len(urInd)):
                volume_toget.append(volumeSpace[m,urInd[m]])
            volume_toget=np.array(volume_toget)
            minInd=np.argmin(np.sum(((calculated_volume_array-volume_toget)/np.mean(volumeSpace,axis=1).reshape((1,-1)))**2.,axis=1))
            if np.all(np.abs(calculated_volume_array[minInd]-volume_toget)<0.0001*np.mean(volumeSpace,axis=1)):
                shutil.copy(transferFrom+"/multipcu_LVBehavior/"+str(minInd)+".log", self.casename+"/"+str(self.runCount)+"/multipcu_LVBehavior/"+str(curr_volN)+".log")
                shutil.copy(transferFrom+"/multipcu_LVBehavior/"+str(minInd)+".npy", self.casename+"/"+str(self.runCount)+"/multipcu_LVBehavior/"+str(curr_volN)+".npy")
            else:
                logger.info("Previous Solution of "+repr(calculated_volume_array[minInd])+" not transfered.\n              Nearest required volume is "+repr(volume_toget))
    def unloadedGeometryRun(self,editParameters=None,inverseHeart=True,tryPressureRatio=0.2,unloadedVolume=None):
        runParameters=self.defaultParameters.duplicate()
        if editParameters is not None:
            for key in runParameters:
                if key in editParameters:
                    runParameters[key]=editParameters[key]
        else:
            editParameters={}
        runParameters['maximum LV fiber tension in Pa']=0.0
        self.writeRunParameters(self.casename+"/"+str(self.runCount),runParameters)
        
        #if self.runCount<=4:
        #    return runParameters
        self.generateMesh(runParameters,toRunCountFolder=True)
        heartRef=heartModelTranslator.HeartModel(self.heartModelStr,self.defaultParameters+{'path of folder with case':self.casename+"/"+str(self.runCount),'name of mesh':self.meshname})
        volRef=heartRef.get_CavityVolume()
        heart=self.getUnloadedGeometry(editParameters=runParameters,casename=self.casename+"/"+str(self.runCount),savename=self.meshname+'_unloadedmesh',toRunCountFolder=False,inverseHeart=inverseHeart,tryPressureRatio=tryPressureRatio,unloadedVolume=unloadedVolume)
        if self.heartModelStr=='heArt.BiV':
            self.solveVolume(heart,np.array([runParameters['LV end diastolic volume in mL'],runParameters['RV end diastolic volume in mL']]),voladj=runParameters['solver volume adjustment'])
        elif self.heartModelStr=='AVHeart':
            if runParameters['Mitral valve channel control'] == False:
                self.solveVolume(heart,np.array([runParameters['LV end diastolic volume in mL'],runParameters['LA end diastolic volume in mL']]),voladj=runParameters['solver volume adjustment'])
            else:
                self.solveVolume(heart,np.array([runParameters['LV end diastolic volume in mL'],runParameters['LA end diastolic volume in mL'],runParameters['Mitral end diastolic volume in mL']]),voladj=runParameters['solver volume adjustment'])
        else:
            self.solveVolume(heart,runParameters['LV end diastolic volume in mL'],voladj=runParameters['solver volume adjustment'])
        if isinstance(runParameters['LV geometrical fiber direction variation file'],str):
            temp_fiber_vec_file=runParameters['LV geometrical fiber direction variation file'][:-4]+'_tempadjmesh'+runParameters['LV geometrical fiber direction variation file'][-4:]
            fiber_temp=heart.get_StretchTensorFunction()
            temp_u=heart.get_DisplacementFunction()
            old_fiber=np.loadtxt(runParameters['LV geometrical fiber direction variation file'])
            temp_probes_f=np.empty((old_fiber.shape[0],9),dtype=float)
            temp_probes_u=np.empty((old_fiber.shape[0],3),dtype=float)
            fiber_temp=fenics.project(fiber_temp,fenics.TensorFunctionSpace(heart.get_Mesh(),'CG',2))
            fiber_temp.set_allow_extrapolation(True)
            temp_u=fenics.project(temp_u,fenics.VectorFunctionSpace(heart.get_Mesh(),'CG',2))
            temp_u.set_allow_extrapolation(True)
            for n in range(old_fiber.shape[0]):
                fiber_temp.eval(temp_probes_f[n],old_fiber[n,:3])
                temp_u.eval(temp_probes_u[n],old_fiber[n,:3])
            temp_probes_f=temp_probes_f.reshape((-1,3,3))
            temp_probes_f=np.matmul(temp_probes_f,old_fiber[:,:3].reshape((-1,3,1))).reshape((-1,3))
            temp_probes_f=temp_probes_f/np.linalg.norm(temp_probes_f,axis=-1,keepdims=True)
            old_fiber[:,:3]=old_fiber[:,:3]+temp_probes_u
            old_fiber[:,3:]=temp_probes_f[:]
            np.savetxt(temp_fiber_vec_file,old_fiber)
        else:
            temp_fiber_vec_file=runParameters['LV geometrical fiber direction variation file']
        if isinstance(runParameters['LV geometrical fiber delay variation file'],str):
            temp_u=heart.get_DisplacementFunction()
            old_fiber=np.loadtxt(runParameters['LV geometrical fiber delay variation file'])
            temp_probes_u=np.empty((old_fiber.shape[0],3),dtype=float)
            temp_u=fenics.project(temp_u,fenics.VectorFunctionSpace(heart.get_Mesh(),'CG',2))
            temp_u.set_allow_extrapolation(True)
            for n in range(old_fiber.shape[0]):
                temp_u.eval(temp_probes_u[n],old_fiber[n,:3])
            old_fiber[:,:3]=old_fiber[:,:3]+temp_probes_u
            temp_fiber_delay_file=runParameters['LV geometrical fiber delay variation file'][:-4]+'_unloadedmesh'+runParameters['LV geometrical fiber delay variation file'][-4:]
            np.savetxt(temp_fiber_delay_file,old_fiber)
        else:
            temp_fiber_delay_file=runParameters['LV geometrical fiber delay variation file']
        if isinstance(runParameters['LV geometrical fiber strength variation file'],str):
            temp_u=heart.get_DisplacementFunction()
            old_fiber=np.loadtxt(runParameters['LV geometrical fiber strength variation file'])
            temp_probes_u=np.empty((old_fiber.shape[0],3),dtype=float)
            temp_u=fenics.project(temp_u,fenics.VectorFunctionSpace(heart.get_Mesh(),'CG',2))
            temp_u.set_allow_extrapolation(True)
            for n in range(old_fiber.shape[0]):
                temp_u.eval(temp_probes_u[n],old_fiber[n,:3])
            old_fiber[:,:3]=old_fiber[:,:3]+temp_probes_u
            temp_fiber_strength_file=runParameters['LV geometrical fiber strength variation file'][:-4]+'_tempadjmesh'+runParameters['LV geometrical fiber strength variation file'][-4:]
            np.savetxt(temp_fiber_strength_file,old_fiber)
        else:
            temp_fiber_strength_file=runParameters['LV geometrical fiber strength variation file']
        if isinstance(runParameters['LV geometrical overall stiffness constant variation file'],str):
            temp_u=heart.get_DisplacementFunction()
            old_fiber=np.loadtxt(runParameters['LV geometrical overall stiffness constant variation file'])
            temp_probes_u=np.empty((old_fiber.shape[0],3),dtype=float)
            temp_u=fenics.project(temp_u,fenics.VectorFunctionSpace(heart.get_Mesh(),'CG',2))
            temp_u.set_allow_extrapolation(True)
            for n in range(old_fiber.shape[0]):
                temp_u.eval(temp_probes_u[n],old_fiber[n,:3])
            old_fiber[:,:3]=old_fiber[:,:3]+temp_probes_u
            temp_fiber_stffvar_file=runParameters['LV geometrical overall stiffness constant variation file'][:-4]+'_tempadjmesh'+runParameters['LV geometrical overall stiffness constant variation file'][-4:]
            np.savetxt(temp_fiber_stffvar_file,old_fiber)
        else:
            temp_fiber_stffvar_file=runParameters['LV geometrical overall stiffness constant variation file']
        fenics.ALE.move(heart.get_Mesh(),heart.get_DisplacementResult())
        np.savetxt(self.casename+"/"+str(self.runCount)+"/coordsDiaMesh.txt",heart.get_MeshCoordinates())
        self.savehdf5(heart,self.casename+"/"+str(self.runCount),self.meshname+'_unloadedmesh',self.meshname+'_unloadedmesh_atDiastole')
        heartModelTranslator.changeFiberAngles(self.heartModelStr,runParameters+{'path of folder with case':self.casename+"/"+str(self.runCount),'name of mesh':self.meshname+'_unloadedmesh_atDiastole','subfolder to save files':'','subfolder to load files':'','LV geometrical fiber direction variation file':temp_fiber_vec_file,'LV geometrical fiber delay variation file':temp_fiber_delay_file,'LV geometrical fiber strength variation file':temp_fiber_strength_file,'LV geometrical overall stiffness constant variation file':temp_fiber_stffvar_file})
        heart=heartModelTranslator.HeartModel(self.heartModelStr,runParameters+{'path of folder with case':self.casename+"/"+str(self.runCount),'name of mesh':self.meshname+'_unloadedmesh'})
        self.solveVolume(heart,volRef,voladj=runParameters['solver volume adjustment'])
        if isinstance(runParameters['LV geometrical fiber direction variation file'],str):
            fiber_temp=heart.get_StretchTensorFunction()
            temp_u=heart.get_DisplacementFunction()
            old_fiber=np.loadtxt(runParameters['LV geometrical fiber direction variation file'])
            temp_probes_f=np.empty((old_fiber.shape[0],9),dtype=float)
            temp_probes_u=np.empty((old_fiber.shape[0],3),dtype=float)
            fiber_temp=fenics.project(fiber_temp,fenics.TensorFunctionSpace(heart.get_Mesh(),'CG',2))
            fiber_temp.set_allow_extrapolation(True)
            temp_u=fenics.project(temp_u,fenics.VectorFunctionSpace(heart.get_Mesh(),'CG',2))
            temp_u.set_allow_extrapolation(True)
            for n in range(old_fiber.shape[0]):
                fiber_temp.eval(temp_probes_f[n],old_fiber[n,:3])
                temp_u.eval(temp_probes_u[n],old_fiber[n,:3])
            temp_probes_f=temp_probes_f.reshape((-1,3,3))
            temp_probes_f=np.matmul(temp_probes_f,old_fiber[:,:3].reshape((-1,3,1))).reshape((-1,3))
            old_fiber[:,:3]=old_fiber[:,:3]+temp_probes_u
            old_fiber[:,3:]=temp_probes_f[:]
            temp_fiber_vec_file=runParameters['LV geometrical fiber direction variation file'][:-4]+'_unloadedmesh'+runParameters['LV geometrical fiber direction variation file'][-4:]
            np.savetxt(temp_fiber_vec_file,old_fiber)
        else:
            temp_fiber_vec_file=runParameters['LV geometrical fiber direction variation file']
        if isinstance(runParameters['LV geometrical fiber delay variation file'],str):
            temp_u=heart.get_DisplacementFunction()
            old_fiber=np.loadtxt(runParameters['LV geometrical fiber delay variation file'])
            temp_probes_u=np.empty((old_fiber.shape[0],3),dtype=float)
            temp_u=fenics.project(temp_u,fenics.VectorFunctionSpace(heart.get_Mesh(),'CG',2))
            temp_u.set_allow_extrapolation(True)
            for n in range(old_fiber.shape[0]):
                temp_u.eval(temp_probes_u[n],old_fiber[n,:3])
            old_fiber[:,:3]=old_fiber[:,:3]+temp_probes_u
            temp_fiber_delay_file=runParameters['LV geometrical fiber delay variation file'][:-4]+'_unloadedmesh'+runParameters['LV geometrical fiber delay variation file'][-4:]
            np.savetxt(temp_fiber_delay_file,old_fiber)
        else:
            temp_fiber_delay_file=runParameters['LV geometrical fiber delay variation file']
        if isinstance(runParameters['LV geometrical fiber strength variation file'],str):
            temp_u=heart.get_DisplacementFunction()
            old_fiber=np.loadtxt(runParameters['LV geometrical fiber strength variation file'])
            temp_probes_u=np.empty((old_fiber.shape[0],3),dtype=float)
            temp_u=fenics.project(temp_u,fenics.VectorFunctionSpace(heart.get_Mesh(),'CG',2))
            temp_u.set_allow_extrapolation(True)
            for n in range(old_fiber.shape[0]):
                temp_u.eval(temp_probes_u[n],old_fiber[n,:3])
            old_fiber[:,:3]=old_fiber[:,:3]+temp_probes_u
            temp_fiber_strength_file=runParameters['LV geometrical fiber strength variation file'][:-4]+'_tempadjmesh'+runParameters['LV geometrical fiber strength variation file'][-4:]
            np.savetxt(temp_fiber_strength_file,old_fiber)
        else:
            temp_fiber_strength_file=runParameters['LV geometrical fiber strength variation file']
        if isinstance(runParameters['LV geometrical overall stiffness constant variation file'],str):
            temp_u=heart.get_DisplacementFunction()
            old_fiber=np.loadtxt(runParameters['LV geometrical overall stiffness constant variation file'])
            temp_probes_u=np.empty((old_fiber.shape[0],3),dtype=float)
            temp_u=fenics.project(temp_u,fenics.VectorFunctionSpace(heart.get_Mesh(),'CG',2))
            temp_u.set_allow_extrapolation(True)
            for n in range(old_fiber.shape[0]):
                temp_u.eval(temp_probes_u[n],old_fiber[n,:3])
            old_fiber[:,:3]=old_fiber[:,:3]+temp_probes_u
            temp_fiber_stffvar_file=runParameters['LV geometrical overall stiffness constant variation file'][:-4]+'_tempadjmesh'+runParameters['LV geometrical overall stiffness constant variation file'][-4:]
            np.savetxt(temp_fiber_stffvar_file,old_fiber)
        else:
            temp_fiber_stffvar_file=runParameters['LV geometrical overall stiffness constant variation file']
        fenics.ALE.move(heart.get_Mesh(),heart.get_DisplacementResult())
        self.savehdf5(heart,self.casename+"/"+str(self.runCount),self.meshname+'_unloadedmesh',self.meshname+'_unloadedmesh_at'+self.meshname)
        heartModelTranslator.changeFiberAngles(self.heartModelStr,runParameters+{'path of folder with case':self.casename+"/"+str(self.runCount),'name of mesh':self.meshname+'_unloadedmesh_at'+self.meshname,'subfolder to save files':'','subfolder to load files':'','LV geometrical fiber direction variation file':temp_fiber_vec_file,'LV geometrical fiber delay variation file':temp_fiber_delay_file,'LV geometrical fiber strength variation file':temp_fiber_strength_file,'LV geometrical overall stiffness constant variation file':temp_fiber_stffvar_file})
        return runParameters
    def LVbehaviorRun(self,editParameters=None,unloadGeo=True,tryPressureRatio=0.2,folderToLVbehavior=None,runCycles=10,minESV=None,maxEDV=None,volstep=50,notAlignEDV=None,multicpu=1,skipKnownErrorStepList=None,reuseLVtable=False,getEDPVR=False,unloadedVolume=None,outputResultList=None,behavior_only=False):
        #if using displacement to get unloaded geometry, set unloadGeo='displacement'
        if isinstance(folderToLVbehavior,str):
            if folderToLVbehavior[0]!='/':
                folderToLVbehavior='/'+folderToLVbehavior
        runParameters=self.defaultParameters.duplicate()
        if editParameters is not None:
            for key in runParameters:
                if key in editParameters:
                    runParameters[key]=editParameters[key]
        else:
            editParameters={}
        self.writeRunParameters(self.casename+"/"+str(self.runCount),runParameters)
        if not(reuseLVtable):
            if min(self.defaultParameters.getParameterRelation(editParameters.keys())+[2])<1:
                self.generateMesh(runParameters,toRunCountFolder=True)
            elif not(os.path.isfile(self.casename+'/'+self.meshname+'.hdf5')):
                if not((folderToLVbehavior is not None and min(self.defaultParameters.getParameterRelation(editParameters.keys())+[2])>0) and os.path.isfile(self.casename+folderToLVbehavior+'/'+self.meshname+'_unloadedmesh.hdf5')):
                    self.generateMesh(runParameters,toRunCountFolder=False)
        if unloadGeo:
            meshname=self.meshname+'_unloadedmesh'
        else:
            meshname=self.meshname
        if not(reuseLVtable):
            if unloadGeo:
                if unloadGeo=='displacement' or self.heartModelStr in ['heArt.BiV']:#,'AVHeart']:
                    inverseHeart=False
                else:
                    inverseHeart=True
                if folderToLVbehavior is not None and min(self.defaultParameters.getParameterRelation(editParameters.keys())+[2])>0:
                    if not(os.path.isfile(self.casename+folderToLVbehavior+'/'+self.meshname+'_unloadedmesh.hdf5')):
                        self.getUnloadedGeometry(editParameters=runParameters,savename=self.meshname+'_unloadedmesh',tryPressureRatio=tryPressureRatio,toRunCountFolder=folderToLVbehavior,inverseHeart=inverseHeart,unloadedVolume=unloadedVolume)
                elif min(self.defaultParameters.getParameterRelation(editParameters.keys())+[2])<1:
                    self.getUnloadedGeometry(editParameters=runParameters,casename=self.casename+'/'+str(self.runCount),savename=self.meshname+'_unloadedmesh',tryPressureRatio=tryPressureRatio,inverseHeart=inverseHeart,unloadedVolume=unloadedVolume)
                elif not(os.path.isfile(self.casename+"/"+str(self.runCount)+'/'+self.meshname+'_unloadedmesh.hdf5')):
                    self.getUnloadedGeometry(editParameters=runParameters,savename=self.meshname+'_unloadedmesh',tryPressureRatio=tryPressureRatio,toRunCountFolder=True,inverseHeart=inverseHeart,unloadedVolume=unloadedVolume)
        if isinstance(runParameters['LV geometrical fiber direction variation file'],str):
            shutil.copy(runParameters['LV geometrical fiber direction variation file'][:-4]+'_tempadjmesh'+runParameters['LV geometrical fiber direction variation file'][-4:],runParameters['LV geometrical fiber direction variation file'][:-4]+'_unloadedmesh'+runParameters['LV geometrical fiber direction variation file'][-4:])
            runParameters['LV geometrical fiber direction variation file']=runParameters['LV geometrical fiber direction variation file'][:-4]+'_unloadedmesh'+runParameters['LV geometrical fiber direction variation file'][-4:]
        if isinstance(runParameters['LV geometrical fiber delay variation file'],str):
            shutil.copy(runParameters['LV geometrical fiber delay variation file'][:-4]+'_tempadjmesh'+runParameters['LV geometrical fiber delay variation file'][-4:],runParameters['LV geometrical fiber delay variation file'][:-4]+'_unloadedmesh'+runParameters['LV geometrical fiber delay variation file'][-4:])
            runParameters['LV geometrical fiber delay variation file']=runParameters['LV geometrical fiber delay variation file'][:-4]+'_unloadedmesh'+runParameters['LV geometrical fiber delay variation file'][-4:]
        if isinstance(runParameters['LV geometrical fiber strength variation file'],str):
            shutil.copy(runParameters['LV geometrical fiber strength variation file'][:-4]+'_tempadjmesh'+runParameters['LV geometrical fiber strength variation file'][-4:],runParameters['LV geometrical fiber strength variation file'][:-4]+'_unloadedmesh'+runParameters['LV geometrical fiber strength variation file'][-4:])
            runParameters['LV geometrical fiber strength variation file']=runParameters['LV geometrical fiber strength variation file'][:-4]+'_unloadedmesh'+runParameters['LV geometrical fiber strength variation file'][-4:]
        if isinstance(runParameters['LV geometrical overall stiffness constant variation file'],str):
            shutil.copy(runParameters['LV geometrical overall stiffness constant variation file'][:-4]+'_tempadjmesh'+runParameters['LV geometrical overall stiffness constant variation file'][-4:],runParameters['LV geometrical overall stiffness constant variation file'][:-4]+'_unloadedmesh'+runParameters['LV geometrical overall stiffness constant variation file'][-4:])
            runParameters['LV geometrical overall stiffness constant variation file']=runParameters['LV geometrical overall stiffness constant variation file'][:-4]+'_unloadedmesh'+runParameters['LV geometrical overall stiffness constant variation file'][-4:]
        if self.heartModelStr=='heArt.BiV':
            cavityNames=['LV','RV']
        elif self.heartModelStr=='AVHeart':
            if runParameters['Mitral valve channel control'] == False:
                cavityNames=['LV','LA']
            else:
                cavityNames=['LV','LA','Mitral']
        else:
            cavityNames=['LV']
        if not(reuseLVtable):
            if multicpu>1:
                if getEDPVR:
                    raise Exception("MultiGPU unsupported for ESPVR")
                if skipKnownErrorStepList is None:
                    skipKnownErrorStepList=[]
                elif isinstance(skipKnownErrorStepList,int):
                    skipKnownErrorStepList=[skipKnownErrorStepList]
                skipKnownErrorStepList=np.array(skipKnownErrorStepList).astype(int)
                os.umask(0)
                os.makedirs(self.casename+"/"+str(self.runCount)+"/multipcu_LVBehavior",mode=0o777,exist_ok=True)
                if folderToLVbehavior is None or min(self.defaultParameters.getParameterRelation(editParameters.keys())+[2])<=1:
                    if min(self.defaultParameters.getParameterRelation(editParameters.keys())+[2])<1:
                        multipcu_kwargs={'meshname':meshname,'meshFolder':str(self.runCount),'volstep':volstep,'minESV':minESV,'maxEDV':maxEDV,'toRunCountFolder':True}
                    elif unloadGeo:
                        multipcu_kwargs={'meshname':meshname,'meshFolder':str(self.runCount),'volstep':volstep,'minESV':minESV,'maxEDV':maxEDV,'toRunCountFolder':True}
                    else:
                        multipcu_kwargs={'meshname':meshname,'volstep':volstep,'minESV':minESV,'maxEDV':maxEDV,'toRunCountFolder':True}
                elif not(os.path.isfile(self.casename+folderToLVbehavior+'/'+meshname+"_Press_VolTime.txt")):
                    if unloadGeo:
                        multipcu_kwargs={'meshname':meshname,'meshFolder':folderToLVbehavior,'volstep':volstep,'minESV':minESV,'maxEDV':maxEDV,'toRunCountFolder':folderToLVbehavior}
                    else:
                        multipcu_kwargs={'meshname':meshname,'volstep':volstep,'minESV':minESV,'maxEDV':maxEDV,'toRunCountFolder':folderToLVbehavior}
                
                volumeSpace,timeSpace=self.getLVbehavior(**multipcu_kwargs,returnVolspaceOnly=True)
                total_volume_to_calculate=volumeSpace.shape[1]**volumeSpace.shape[0]
                
                insert_path_for_heartFEM = sys.path
                if multipcu_kwargs['toRunCountFolder']:
                    if isinstance(multipcu_kwargs['toRunCountFolder'],str):
                        if multipcu_kwargs['toRunCountFolder'][0]!="/":
                            toRunCountFolder_getLVbehavior="/"+multipcu_kwargs['toRunCountFolder']
                        else:
                            toRunCountFolder_getLVbehavior=multipcu_kwargs['toRunCountFolder']
                    else:
                        toRunCountFolder_getLVbehavior="/"+str(self.runCount)
                else:
                    toRunCountFolder_getLVbehavior=""
                rerunall=True
                if os.path.isfile(self.casename+toRunCountFolder_getLVbehavior+"/"+meshname+"_Press_volumeSpace.txt") and os.path.isfile(self.casename+toRunCountFolder_getLVbehavior+"/"+meshname+"_Press_timeSpace.txt"):
                    match_volume=np.loadtxt(self.casename+toRunCountFolder_getLVbehavior+"/"+meshname+"_Press_volumeSpace.txt")
                    if np.all(np.abs(match_volume-volumeSpace)<0.0001*np.mean(volumeSpace,axis=1).reshape((-1,1))):
                        rerunall=False
                    else:
                        temp = input("Volumespace does not match. Rerun all? (y/n)")
                        if not(temp in ["y"," y","y "," y ","yes"," yes","yes "," yes "]):
                            rerunall=False
                para= self.paraToRecreateObj()
                for key in para:
                    multipcu_kwargs['LVclosed_'+key]=para[key]
                import pickle
                print(multipcu_kwargs)
                with open(self.casename+"/"+str(self.runCount)+"/multipcu_LVBehavior/"+'getLVBehavior_dictionary.pkl', 'wb') as f:
                    pickle.dump(multipcu_kwargs, f)
                with open(self.casename+"/"+str(self.runCount)+"/multicpurun.py",'w') as f:
                    fileLines=list(std_multicpu_getLVbehavior_file)
                    if rerunall:
                        fileLines[5]="if False:\n"
                    f.writelines(std_multicpu_getLVbehavior_file)
                multicpu_cmd_func=multicpu_cmd(self.casename+"/"+str(self.runCount),insert_path_for_heartFEM)
                if isinstance(multicpu,int):
                    num = multicpu
                else:
                    num = None  # set to the number of workers you want (it defaults to the cpu count of your machine)
                pool=Pool(num)
                range_total_volume_to_calculate=np.array(range(total_volume_to_calculate))
                range_total_volume_to_calculate=range_total_volume_to_calculate[np.logical_not(np.in1d(range_total_volume_to_calculate,skipKnownErrorStepList))]
                if (len(range_total_volume_to_calculate)+len(skipKnownErrorStepList))!=len(range(total_volume_to_calculate)):
                    raise Exception("Unable to skip "+repr(skipKnownErrorStepList)+" within "+repr(total_volume_to_calculate)+" steps")
                pool.map(multicpu_cmd_func,range_total_volume_to_calculate)
                pool.close()
                pool.join()
                resultdata=[]
                for n in range(total_volume_to_calculate):#goes thru the LV first then RV
                    if os.path.isfile(self.casename+"/"+str(self.runCount)+"/multipcu_LVBehavior/"+str(n)+".npy"):
                        resultdata.append(np.load(self.casename+"/"+str(self.runCount)+"/multipcu_LVBehavior/"+str(n)+".npy"))
                    else:
                        resultdata.append(np.zeros((len(cavityNames),len(timeSpace)))*float('nan'))
                resultdata=np.array(resultdata)
                if np.any(np.isnan(resultdata)):
                    if len(cavityNames)>1:
                        import itertools
                        volume_x=[]
                        for r in itertools.product(*volumeSpace[::-1]):
                            volume_x.append(np.array(r)[::-1])
                        volume_x=np.array(volume_x)
                    else:
                        volume_x=volumeSpace.T
                    train_bool=np.logical_not(np.isnan(resultdata[:,0,0]))
                    volume_xtrain=np.concatenate((np.tile(volume_x[train_bool],(len(timeSpace),1)),np.repeat(timeSpace,sum(train_bool)).reshape((-1,1))),axis=1)
                    volume_xtest=np.concatenate((np.tile(volume_x[np.logical_not(train_bool)],(len(timeSpace),1)),np.repeat(timeSpace,sum(np.logical_not(train_bool))).reshape((-1,1))),axis=1)
                    
                for cavityN in range(len(cavityNames)):
                    temp_resultdata=resultdata[:,cavityN]
                    if np.any(np.isnan(temp_resultdata)):
                        #temp_resultdata[np.logical_not(train_bool)]=scipyinterpolate.RBFInterpolator(volume_xtrain,temp_resultdata[train_bool].reshape(-1,order='F'))(volume_xtest).reshape((-1,len(timeSpace)),order='F')
                        nearest_result_temp=scipyinterpolate.griddata(volume_xtrain, temp_resultdata[train_bool].reshape(-1,order='F'), volume_xtest, method='nearest').reshape((-1,len(timeSpace)),order='F')
                        temp_resultdata[np.logical_not(train_bool)]=scipyinterpolate.griddata(volume_xtrain, temp_resultdata[train_bool].reshape(-1,order='F'), volume_xtest, method='linear').reshape((-1,len(timeSpace)),order='F')
                        temp_resultdata[np.logical_not(train_bool)][np.isnan(temp_resultdata[np.logical_not(train_bool)])] = nearest_result_temp[np.isnan(temp_resultdata[np.logical_not(train_bool)])]
                    press_volTime_base=temp_resultdata[:,0:1]
                    press_volTime=temp_resultdata-press_volTime_base
                    np.savetxt(self.casename+"/"+str(self.runCount)+"/"+meshname+"_Press_VolTime"+"_"+cavityNames[cavityN]+".txt",press_volTime,header=str(total_volume_to_calculate)+'solved rowVolm: '+str(total_volume_to_calculate)+'\n via multicpu. ')
                    np.savetxt(self.casename+"/"+str(self.runCount)+"/"+meshname+"_Press_VolTime"+"_"+cavityNames[cavityN]+"_base"+".txt",press_volTime_base[...,0])
                
            else:
                if getEDPVR:
                    endtime=0
                else:
                    endtime=None
                if folderToLVbehavior is None or min(self.defaultParameters.getParameterRelation(editParameters.keys())+[2])<=1:
                    if min(self.defaultParameters.getParameterRelation(editParameters.keys())+[2])<1:
                        self.getLVbehavior(runParameters=runParameters,meshname=meshname,meshFolder=str(self.runCount),volstep=volstep,minESV=minESV,maxEDV=maxEDV,endtime=endtime,outputResultList=outputResultList,toRunCountFolder=True)
                    elif unloadGeo:
                        self.getLVbehavior(runParameters=runParameters,meshname=meshname,meshFolder=str(self.runCount),volstep=volstep,minESV=minESV,maxEDV=maxEDV,endtime=endtime,outputResultList=outputResultList,toRunCountFolder=True)
                    else:
                        self.getLVbehavior(runParameters=runParameters,meshname=meshname,volstep=volstep,minESV=minESV,maxEDV=maxEDV,endtime=endtime,outputResultList=outputResultList,toRunCountFolder=True)
                elif not(os.path.isfile(self.casename+folderToLVbehavior+'/'+meshname+"_Press_VolTime.txt")):
                    if unloadGeo:
                        self.getLVbehavior(runParameters=runParameters,meshname=meshname,meshFolder=folderToLVbehavior,volstep=volstep,minESV=minESV,maxEDV=maxEDV,endtime=endtime,outputResultList=outputResultList,toRunCountFolder=folderToLVbehavior)
                    else:
                        self.getLVbehavior(runParameters=runParameters,meshname=meshname,volstep=volstep,minESV=minESV,maxEDV=maxEDV,endtime=endtime,outputResultList=outputResultList,toRunCountFolder=folderToLVbehavior)
            if getEDPVR:
                return runParameters
            if folderToLVbehavior is None:
                for funcName in cavityNames:#["LV","RV","LA","RA"]:
                    func=getattr(ngspice_py,"generate"+funcName+"table")
                    func(self.casename+"/"+str(self.runCount)+'/'+meshname,runParameters['duration of one cardiac cycle in ms'],timetopeak=runParameters['time to maximum LV fiber tension in ms'],loading_casename=self.casename+"/"+str(self.runCount)+'/'+meshname,scale_maximum_fiber_tension=runParameters['scale Windkessel '+funcName+' active pressure'])
                    #try:
                    #    func(self.casename+"/"+str(self.runCount)+'/'+meshname,runParameters['duration of one cardiac cycle in ms'],timetopeak=runParameters['time to maximum LV fiber tension in ms'],loading_casename=self.casename+"/"+str(self.runCount)+'/'+meshname,scale_maximum_fiber_tension=runParameters['scale Windkessel active pressure'])
                    #except Exception as e:
                    #    logger.warning(repr(func)+" : "+repr(e))
            else:
                for funcName in cavityNames:#["LV","RV","LA","RA"]:
                    func=getattr(ngspice_py,"generate"+funcName+"table")
                    func(self.casename+"/"+str(self.runCount)+'/'+meshname,runParameters['duration of one cardiac cycle in ms'],timetopeak=runParameters['time to maximum LV fiber tension in ms'],loading_casename=self.casename+folderToLVbehavior+'/'+meshname,scale_maximum_fiber_tension=runParameters['scale Windkessel '+funcName+' active pressure'])
                    #try:
                    #    func(self.casename+"/"+str(self.runCount)+'/'+meshname,runParameters['duration of one cardiac cycle in ms'],timetopeak=runParameters['time to maximum LV fiber tension in ms'],loading_casename=self.casename+folderToLVbehavior+'/'+meshname,scale_maximum_fiber_tension=runParameters['scale Windkessel active pressure'])
                    #except Exception as e:
                    #    logger.warning(repr(func)+" : "+repr(e))
        if behavior_only:
            return runParameters
        if 'timebased' in self.defaultParameters['Windkessel RV source function']:
            cmd = "cp "+self.casename+'/'+self.meshname+'_rvflowrate.txt'+" " + self.casename+"/"+str(self.runCount)+'/'+meshname+'_rvflowrate.txt'
            os.system(cmd)
        
        if runParameters['duration of LV diastole in ms'] is None:
            ngspice_py.createLVcircuit(self.casename+"/"+str(self.runCount)+'/'+meshname,runParameters)
        else:
            ngspice_py.createLVcircuit(self.casename+"/"+str(self.runCount)+'/'+meshname,runParameters,skipVariableList=["timetopeaktension"])
        if notAlignEDV is None:
            if self.heartModelStr=='heArt.BiV':
                ngspice_py.simLVcircuit_align_EDvol_and_EStime(self.casename+"/"+str(self.runCount)+'/'+meshname,runParameters['duration of one cardiac cycle in ms']*runCycles,runParameters+{'Windkessel LV source file':self.casename+"/"+str(self.runCount)+'/'+meshname+'_lvcirtable.txt','Windkessel RV source file':self.casename+"/"+str(self.runCount)+'/'+meshname+'_rvcirtable.txt'},runParameters['duration of one cardiac cycle in ms'],runParameters['duration of LV diastole in ms'],runParameters['time to maximum LV fiber tension in ms'],initLVvol=runParameters['LV end diastolic volume in mL'],try_initRVvol=True,initRVvol=runParameters['RV end diastolic volume in mL'],vla0=runParameters['vla0'],vra0=runParameters['vra0'])
            elif self.heartModelStr=='AVHeart':
                ngspice_py.simLVcircuit_align_EDvol_and_EStime(self.casename+"/"+str(self.runCount)+'/'+meshname,runParameters['duration of one cardiac cycle in ms']*runCycles,runParameters+{'Windkessel LV source file':self.casename+"/"+str(self.runCount)+'/'+meshname+'_lvcirtable.txt','Windkessel LA source file':self.casename+"/"+str(self.runCount)+'/'+meshname+'_rvcirtable.txt'},runParameters['duration of one cardiac cycle in ms'],runParameters['duration of LV diastole in ms'],runParameters['time to maximum LV fiber tension in ms'],initLVvol=runParameters['LV end diastolic volume in mL'],vla0=runParameters['vla0'],vra0=runParameters['vra0'])
            else:
                ngspice_py.simLVcircuit_align_EDvol_and_EStime(self.casename+"/"+str(self.runCount)+'/'+meshname,runParameters['duration of one cardiac cycle in ms']*runCycles,runParameters+{'Windkessel LV source file':self.casename+"/"+str(self.runCount)+'/'+meshname+'_lvcirtable.txt'},runParameters['duration of one cardiac cycle in ms'],runParameters['duration of LV diastole in ms'],runParameters['time to maximum LV fiber tension in ms'],initLVvol=runParameters['LV end diastolic volume in mL'],vla0=runParameters['vla0'],vra0=runParameters['vra0'])
        else:
            ngspice_py.simLVcircuit(self.casename+"/"+str(self.runCount)+'/'+meshname,runParameters['duration of one cardiac cycle in ms']*runCycles,runParameters+{'Windkessel LV source file':self.casename+"/"+str(self.runCount)+'/'+meshname+'_lvcirtable.txt'},initLVvol=notAlignEDV,vla0=runParameters['vla0'],vra0=runParameters['vra0'])
        self.getLastcycleCircuitResults(runParameters['duration of one cardiac cycle in ms'])
        return runParameters
        
    def flowrateWindkesselRun(self,editParameters=None,lvufile=None,lvinputvar=None,stepTime=10):
        #if using displacement to get unloaded geometry, set unloadGeo='displacement'
        runParameters=self.defaultParameters.duplicate()
        if editParameters is not None:
            for key in runParameters:
                if key in editParameters:
                    runParameters[key]=editParameters[key]
        else:
            editParameters={}
        self.writeRunParameters(self.casename+"/"+str(self.runCount),runParameters)
        if lvufile is None:
            lvufile=self.casename+"/"+str(self.runCount)+'/'+self.meshname+'_lvu_outflowrate.txt'
        runParameters+={'Windkessel LV source function':'time based current','Windkessel LV source file':lvufile}
        if not(os.path.isfile(lvufile)):
            ref_time=np.arange(-1,self.defaultParameters['duration of one cardiac cycle in ms']*10.1+1)
            ref_volume=self.manualVol(ref_time)
            ref_outflowrate=(ref_volume[2:]-ref_volume[:-2])/2.
            tempdata=np.array([ref_time[1:-1],ref_outflowrate]).T
            np.savetxt(lvufile,tempdata)
        ngspice_py.createLVcircuit(self.casename+"/"+str(self.runCount)+'/'+self.meshname,runParameters,stepTime=stepTime)
        ngspice_py.simLVcircuit(self.casename+"/"+str(self.runCount)+'/'+self.meshname,runParameters['duration of one cardiac cycle in ms']*10,runParameters,initLVvol=runParameters['LV end diastolic volume in mL'],vla0=runParameters['vla0'],vra0=runParameters['vra0'])
        self.getLastcycleCircuitResults(runParameters['duration of one cardiac cycle in ms'])
        return runParameters
    def fullWindkesselRun(self,editParameters=None,unloadGeo=True,folderToLVbehavior=None,runCycles=10,notAlignEDV=None):
        #if using displacement to get unloaded geometry, set unloadGeo='displacement'
        if isinstance(folderToLVbehavior,str):
            if folderToLVbehavior[0]!='/':
                folderToLVbehavior='/'+folderToLVbehavior
        runParameters=self.defaultParameters.duplicate()
        if editParameters is not None:
            for key in runParameters:
                if key in editParameters:
                    runParameters[key]=editParameters[key]
        else:
            editParameters={}
        self.writeRunParameters(self.casename+"/"+str(self.runCount),runParameters)
        if unloadGeo:
            meshname=self.meshname+'_unloadedmesh'
        else:
            meshname=self.meshname
        
        if folderToLVbehavior is None:
            ngspice_py.generateLVtable(self.casename+"/"+str(self.runCount)+'/'+meshname,runParameters['duration of one cardiac cycle in ms'],timetopeak=runParameters['time to maximum LV fiber tension in ms'],loading_casename=self.casename+"/"+str(self.runCount)+'/'+meshname,scale_maximum_fiber_tension=runParameters['scale Windkessel active pressure'])
        else:
            ngspice_py.generateLVtable(self.casename+"/"+str(self.runCount)+'/'+meshname,runParameters['duration of one cardiac cycle in ms'],timetopeak=runParameters['time to maximum LV fiber tension in ms'],loading_casename=self.casename+folderToLVbehavior+'/'+meshname,scale_maximum_fiber_tension=runParameters['scale Windkessel active pressure'])
        cmd = "cp "+self.casename+'/'+self.meshname+'_rvflowrate.txt'+" " + self.casename+"/"+str(self.runCount)+'/'+meshname+'_rvflowrate.txt'
        os.system(cmd)
        
        if runParameters['duration of LV diastole in ms'] is None:
            ngspice_py.createLVcircuit(self.casename+"/"+str(self.runCount)+'/'+meshname,runParameters)
        else:
            ngspice_py.createLVcircuit(self.casename+"/"+str(self.runCount)+'/'+meshname,runParameters,skipVariableList=["timetopeaktension"])
        if notAlignEDV is None:
            ngspice_py.simLVcircuit_align_EDvol_and_EStime(self.casename+"/"+str(self.runCount)+'/'+meshname,runParameters['duration of one cardiac cycle in ms']*runCycles,runParameters+{'Windkessel LV source file':self.casename+"/"+str(self.runCount)+'/'+meshname+'_lvcirtable.txt'},runParameters['duration of one cardiac cycle in ms'],runParameters['duration of LV diastole in ms'],runParameters['time to maximum LV fiber tension in ms'],initLVvol=runParameters['LV end diastolic volume in mL'],vla0=runParameters['vla0'],vra0=runParameters['vra0'])
        else:
            ngspice_py.simLVcircuit(self.casename+"/"+str(self.runCount)+'/'+meshname,runParameters['duration of one cardiac cycle in ms']*runCycles,runParameters+{'Windkessel LV source file':self.casename+"/"+str(self.runCount)+'/'+meshname+'_lvcirtable.txt'},initLVvol=notAlignEDV,vla0=runParameters['vla0'],vra0=runParameters['vra0'])
        self.getLastcycleCircuitResults(runParameters['duration of one cardiac cycle in ms'])
        return runParameters
    def iterativeRun(self,editParameters=None,runTimeList=None,endTime=None,manualPhaseTimeInput=False,unloadGeo=True,setHeart=None,outputResultList=None,tryPressureRatio=0.2,trackphase=False,electroSim=False,unloadedVolume=None,unload_loop_num=float('inf'),continue_from_break=False,multicpu=1):
        '''
        setHeart: list or self.heartModel object: list [folder name, meshname]
        '''
        if electroSim and continue_from_break:
            raise Exception("'continue_from_break' is not implemented for electro-Simulation.")
        if self.heartModelStr=='heArt.BiV':
            project_form_compiler_parameters={"representation":"uflacs","quadrature_degree":4}
        else:
            project_form_compiler_parameters={"representation":"uflacs"}
        runParameters=self.defaultParameters.duplicate()
        if editParameters is not None:
            if "read_from_file" in editParameters:
                runParameters=self.readRunParameters(editParameters["read_from_file"])
            for key in runParameters:
                if key in editParameters:
                    runParameters[key]=editParameters[key]
        else:
            editParameters={}
        if continue_from_break:
            runParameters=self.readRunParameters(self.casename+"/"+str(self.runCount))
        else:
            self.writeRunParameters(self.casename+"/"+str(self.runCount),runParameters)
        if runTimeList is not None or manualPhaseTimeInput:
            if trackphase:
                logger.warning("trackphase is not compatible with runTimeList set or manualPhaseTimeInput, turning trackphase off")
                trackphase=False
        if setHeart is None and not(continue_from_break):
            if min(self.defaultParameters.getParameterRelation(editParameters.keys())+[2])<1:
                self.generateMesh(runParameters,toRunCountFolder=True)
            elif not(os.path.isfile(self.casename+'/'+self.meshname+'.hdf5')):
                self.generateMesh(runParameters,toRunCountFolder=False)
        if outputResultList is None:
            outputResultList=['displacement','stress']
        elif isinstance(outputResultList,str):
            outputResultList=[outputResultList]
        if electroSim:
            if not("activation_delay" in outputResultList):
                outputResultList+=["activation_delay"]
        
        resultWriter=fenicsResultWriter(self.casename+"/"+str(self.runCount),outputResultList,project_form_compiler_parameters=project_form_compiler_parameters)
        #displacementfile = fenics.File(self.casename+"/"+str(self.runCount)+"/deformation/u_disp.pvd")
        #stress_File = fenics.File(self.casename+"/"+str(self.runCount)+"/stress/_stress.pvd") #Joy Changed here

        
        ########################### Fenics's Newton  #########################################################

        # Closed loop cycle
        BCL = runParameters['duration of one cardiac cycle in ms']
  
        cycle = 0.
        dt=1.
        if runTimeList is None:
            tstep = 0
            t_array=[0]
        else:
            tstep = runTimeList[0]
            t_array=[runTimeList[0]]
        
        

        ####### Loading phase for LV ####################################################
        if self.manualVol is None:
            if self.heartModelStr=='heArt.BiV':
                startVolume=np.array([runParameters['LV end diastolic volume in mL'],runParameters['RV end diastolic volume in mL']])
            elif self.heartModelStr=='AVHeart':
                if runParameters['Mitral valve channel control'] == False:
                    startVolume=np.array([runParameters['LV end diastolic volume in mL'],runParameters['LA end diastolic volume in mL']])
                else:
                    startVolume=np.array([runParameters['LV end diastolic volume in mL'],runParameters['LA end diastolic volume in mL'],runParameters['Mitral end diastolic volume in mL']])
            else:
                startVolume=np.array([runParameters['LV end diastolic volume in mL']])
        elif runTimeList is None:
            startVolume=self.manualVol(0.)
        else:
            startVolume=self.manualVol(runTimeList[0])
        
        #if Cavityvol.vol>(1.1*runParameters['LV end diastolic volume in mL']):
        # 	raise Exception('Starting cavity volume,'+str(Cavityvol.vol)+', is more than target EDV_LV,'+str(runParameters['LV end diastolic volume in mL']))
        
        if setHeart is not None:
            if isinstance(setHeart,list):
                heart=heartModelTranslator.HeartModel(self.heartModelStr,runParameters+{'path of folder with case':setHeart[0],'name of mesh':setHeart[1],'fiber relaxation based on phase during Windkessel':trackphase})
                if electroSim:
                    electroPHeart=heartModelTranslator.HeartModel('heArt.Electro',runParameters+{'path of folder with case':setHeart[0],'name of mesh':setHeart[1]})
            else:
                heart=setHeart
            
        elif not(unloadGeo):
            if min(self.defaultParameters.getParameterRelation(editParameters.keys())+[2])<1:
                heart=heartModelTranslator.HeartModel(self.heartModelStr,runParameters+{'path of folder with case':self.casename+'/'+str(self.runCount),'name of mesh':self.meshname,'fiber relaxation based on phase during Windkessel':trackphase})
                if electroSim:
                    electroPHeart=heartModelTranslator.HeartModel('heArt.Electro',runParameters+{'path of folder with case':self.casename+'/'+str(self.runCount),'name of mesh':self.meshname})
            else:
                heart=heartModelTranslator.HeartModel(self.heartModelStr,runParameters+{'path of folder with case':self.casename,'name of mesh':self.meshname,'fiber relaxation based on phase during Windkessel':trackphase})
                if electroSim:
                    electroPHeart=heartModelTranslator.HeartModel('heArt.Electro',runParameters+{'path of folder with case':self.casename,'name of mesh':self.meshname})
        else:
            copy_file_for_ref=False
            if min(self.defaultParameters.getParameterRelation(editParameters.keys())+[2])<1:
                if not(continue_from_break):
                    self.getUnloadedGeometry(editParameters=runParameters,casename=self.casename+'/'+str(self.runCount),savename=self.meshname+'_unloadedmesh',tryPressureRatio=tryPressureRatio,unloadedVolume=unloadedVolume,unload_loop_num=unload_loop_num)
                    copy_file_for_ref=True
                heart=heartModelTranslator.HeartModel(self.heartModelStr,runParameters+{'path of folder with case':self.casename+'/'+str(self.runCount),'name of mesh':self.meshname+'_unloadedmesh','fiber relaxation based on phase during Windkessel':trackphase})
                if electroSim:
                    electroPHeart=heartModelTranslator.HeartModel('heArt.Electro',runParameters+{'path of folder with case':self.casename+'/'+str(self.runCount),'name of mesh':self.meshname+'_unloadedmesh'})
            else:
                if not(os.path.isfile(self.casename+"/"+self.meshname+'_unloadedmesh.hdf5')) and not(continue_from_break):
                    self.getUnloadedGeometry(editParameters=runParameters,savename=self.meshname+'_unloadedmesh',toRunCountFolder=False,tryPressureRatio=tryPressureRatio,unloadedVolume=unloadedVolume,unload_loop_num=unload_loop_num)
                    copy_file_for_ref=True
                heart=heartModelTranslator.HeartModel(self.heartModelStr,runParameters+{'path of folder with case':self.casename,'name of mesh':self.meshname+'_unloadedmesh','fiber relaxation based on phase during Windkessel':trackphase})
                if electroSim:
                    electroPHeart=heartModelTranslator.HeartModel('heArt.Electro',runParameters+{'path of folder with case':self.casename,'name of mesh':self.meshname+'_unloadedmesh'})
            if copy_file_for_ref:   
                if isinstance(runParameters['LV geometrical fiber direction variation file'],str):
                    shutil.copy(runParameters['LV geometrical fiber direction variation file'][:-4]+'_tempadjmesh'+runParameters['LV geometrical fiber direction variation file'][-4:],runParameters['LV geometrical fiber direction variation file'][:-4]+'_unloadedmesh'+runParameters['LV geometrical fiber direction variation file'][-4:])
                    runParameters['LV geometrical fiber direction variation file']=runParameters['LV geometrical fiber direction variation file'][:-4]+'_unloadedmesh'+runParameters['LV geometrical fiber direction variation file'][-4:]
                if isinstance(runParameters['LV geometrical fiber delay variation file'],str):
                    shutil.copy(runParameters['LV geometrical fiber delay variation file'][:-4]+'_tempadjmesh'+runParameters['LV geometrical fiber delay variation file'][-4:],runParameters['LV geometrical fiber delay variation file'][:-4]+'_unloadedmesh'+runParameters['LV geometrical fiber delay variation file'][-4:])
                    runParameters['LV geometrical fiber delay variation file']=runParameters['LV geometrical fiber delay variation file'][:-4]+'_unloadedmesh'+runParameters['LV geometrical fiber delay variation file'][-4:]
                if isinstance(runParameters['LV geometrical fiber strength variation file'],str):
                    shutil.copy(runParameters['LV geometrical fiber strength variation file'][:-4]+'_tempadjmesh'+runParameters['LV geometrical fiber strength variation file'][-4:],runParameters['LV geometrical fiber strength variation file'][:-4]+'_unloadedmesh'+runParameters['LV geometrical fiber strength variation file'][-4:])
                    runParameters['LV geometrical fiber strength variation file']=runParameters['LV geometrical fiber strength variation file'][:-4]+'_unloadedmesh'+runParameters['LV geometrical fiber strength variation file'][-4:]
                if isinstance(runParameters['LV geometrical overall stiffness constant variation file'],str):
                    shutil.copy(runParameters['LV geometrical overall stiffness constant variation file'][:-4]+'_tempadjmesh'+runParameters['LV geometrical overall stiffness constant variation file'][-4:],runParameters['LV geometrical overall stiffness constant variation file'][:-4]+'_unloadedmesh'+runParameters['LV geometrical overall stiffness constant variation file'][-4:])
                    runParameters['LV geometrical overall stiffness constant variation file']=runParameters['LV geometrical overall stiffness constant variation file'][:-4]+'_unloadedmesh'+runParameters['LV geometrical overall stiffness constant variation file'][-4:]
            del copy_file_for_ref
        heart.set_dTime(1.)
        if multicpu>1:
            if runTimeList is not None:
                runTimeList=list(runTimeList)
            else:
                t = 0
                runTimeList=[0.]
                while(runTimeList[-1] < self.numofCycle*BCL):
                    t = runTimeList[-1] - math.floor(runTimeList[-1]/BCL)*BCL
                    if(t >= 0.0 and t < 4.0):
                        runTimeList.append(runTimeList[-1]+0.50)
                    elif (t >= (runParameters['time to maximum LV fiber tension in ms']-0.5) and t < (runParameters['time to maximum LV fiber tension in ms']+0.5)):
                        runTimeList.append(runTimeList[-1]+3)
                    else :
                        runTimeList.append(runTimeList[-1]+1.0)
            splitruntimelist_num=math.ceil(len(runTimeList)/multicpu)
            splitruntimelist=[]
            for runtime in range(multicpu):
                splitruntimelist.append(runTimeList[:splitruntimelist_num])
                del runTimeList[:splitruntimelist_num]
            os.umask(0)
            os.makedirs(self.casename+"/"+str(self.runCount)+"/multipcu_iterative",mode=0o777,exist_ok=True)
            
            if electroSim:
                raise Exception("electroSim not implemented.")
            if not(unloadGeo):
                if min(self.defaultParameters.getParameterRelation(editParameters.keys())+[2])<1:
                    heartModel_set=[self.casename+'/'+str(self.runCount),self.meshname]
                else:
                    heartModel_set=[self.casename,self.meshname]
                    
            else:
                if min(self.defaultParameters.getParameterRelation(editParameters.keys())+[2])<1:
                    heartModel_set=[self.casename+'/'+str(self.runCount),self.meshname+'_unloadedmesh']
                else:
                    heartModel_set=[self.casename,self.meshname+'_unloadedmesh']
            insert_path_for_heartFEM = sys.path
            
            rerunall=False
            for n in range(multicpu):
                os.makedirs(self.casename+"/"+str(self.runCount)+"/multipcu_iterative/"+str(n),mode=0o777,exist_ok=True)
                self.writeRunParameters(self.casename+"/"+str(self.runCount)+"/multipcu_iterative/"+str(n),runParameters)
                shutil.copy(heartModel_set[0]+"/"+heartModel_set[1]+".hdf5",self.casename+"/"+str(self.runCount)+"/multipcu_iterative/"+str(n)+"/"+heartModel_set[1]+".hdf5")
                if heartModel_set[1]!=self.meshname:
                    shutil.copy(heartModel_set[0]+"/"+self.meshname+".hdf5",self.casename+"/"+str(self.runCount)+"/multipcu_iterative/"+str(n)+"/"+self.meshname+".hdf5")
                if os.path.isfile(self.casename+"/"+str(self.runCount)+"/multipcu_iterative/"+str(n)+"/runtimelist.txt"):
                    match_runTimeList=np.loadtxt(self.casename+"/"+str(self.runCount)+"/multipcu_iterative/"+str(n)+"/runtimelist.txt")
                    if not(np.all(match_runTimeList==np.array(splitruntimelist[n]))):
                        rerunall=True   
                else:
                    rerunall=True
            if rerunall:
                for n in range(multicpu):
                    np.savetxt(self.casename+"/"+str(self.runCount)+"/multipcu_iterative/"+str(n)+"/runtimelist.txt",splitruntimelist[n])

            para= self.paraToRecreateObj()
            multipcu_kwargs={}
            for key in para:
                multipcu_kwargs['LVclosed_'+key]=para[key]
            multipcu_kwargs['LVclosed_casename']=self.casename+"/"+str(self.runCount)+"/multipcu_iterative"
            multipcu_kwargs['outputResultList']=outputResultList
            print(multipcu_kwargs)
            import pickle
            with open(self.casename+"/"+str(self.runCount)+"/multipcu_iterative/"+'iterativerun_dictionary.pkl', 'wb') as f:
                pickle.dump(multipcu_kwargs, f)
            with open(self.casename+"/"+str(self.runCount)+"/multicpurun.py",'w') as f:
                fileLines=list(std_multicpu_iterativerun_file)
                if rerunall:
                    fileLines[5]="if False:\n"
                f.writelines(std_multicpu_iterativerun_file)
            multicpu_cmd_func=multicpu_cmd(self.casename+"/"+str(self.runCount),insert_path_for_heartFEM)
            
            
            pool=Pool(multicpu)
            pool.map(multicpu_cmd_func,np.array(range(multicpu)))
            pool.close()
            pool.join()
            #combine results here!!
            os.makedirs(self.casename+"/"+str(self.runCount)+"/multipcu_iterative/combined",exist_ok=True)
            countstart_temp=[0]
            for n in range(multicpu):
                countstart_temp.append(len(splitruntimelist[n])+countstart_temp[-1])
            temp_outputs={}
            PV_all=np.loadtxt(self.casename+"/"+str(self.runCount)+"/multipcu_iterative/0/"+str(self.runCount)+"/PV_.txt")
            for n in range(1,multicpu):
                PV_all=np.concatenate((PV_all,np.loadtxt(self.casename+"/"+str(self.runCount)+"/multipcu_iterative/"+str(n)+"/"+str(self.runCount)+"/PV_.txt")),axis=0)
            np.savetxt(self.casename+"/"+str(self.runCount)+"/multipcu_iterative/combined/PV_all.txt",PV_all)
            dir_names=os.listdir(self.casename+"/"+str(self.runCount)+"/multipcu_iterative/0/"+str(self.runCount))
            outputResultDir=[]
            for temp_name in dir_names:
                if os.path.isdir(self.casename+"/"+str(self.runCount)+"/multipcu_iterative/0/"+str(self.runCount)+"/"+temp_name):
                    outputResultDir.append(temp_name)
            for n in range(multicpu):
                for temp_name in outputResultDir:
                    if n==0:
                        if os.path.isdir(self.casename+"/"+str(self.runCount)+"/multipcu_iterative/combined/"+temp_name):
                            shutil.rmtree(self.casename+"/"+str(self.runCount)+"/multipcu_iterative/combined/"+temp_name)
                        shutil.copytree(self.casename+"/"+str(self.runCount)+"/multipcu_iterative/0/"+str(self.runCount)+"/"+temp_name,self.casename+"/"+str(self.runCount)+"/multipcu_iterative/combined/"+temp_name)
                        temp_list=os.listdir(self.casename+"/"+str(self.runCount)+"/multipcu_iterative/0/"+str(self.runCount)+"/"+temp_name)
                        temp_outputs[temp_name]=[]
                        for filelist_temp in temp_list:
                            if filelist_temp[-4:]==".pvd":
                                temp_outputs[temp_name].append(filelist_temp[:-4])
                                with open(self.casename+"/"+str(self.runCount)+"/multipcu_iterative/0/"+str(self.runCount)+"/"+temp_name+"/"+filelist_temp,'r') as f:
                                    temp_pvd=f.readlines()
                                for fcount in range(countstart_temp[1],countstart_temp[-1]):
                                    temp_pvd.insert(-2,std_pvd_line.format(filelist_temp[-4:],fcount))
                                with open(self.casename+"/"+str(self.runCount)+"/multipcu_iterative/combined/"+temp_name+"/"+filelist_temp,'w') as f:
                                    f.writelines(temp_pvd)
                    else:
                        for filelist_temp in temp_outputs[temp_name]:
                            for m in range(len(splitruntimelist[n])):
                                if not(os.path.isfile(self.casename+"/"+str(self.runCount)+"/multipcu_iterative/"+str(n)+"/"+str(self.runCount)+"/"+temp_name+"/"+filelist_temp+"{0:06d}.vtu".format(m))):
                                    break
                                shutil.copy(self.casename+"/"+str(self.runCount)+"/multipcu_iterative/"+str(n)+"/"+str(self.runCount)+"/"+temp_name+"/"+filelist_temp+"{0:06d}.vtu".format(m),self.casename+"/"+str(self.runCount)+"/multipcu_iterative/combined/"+temp_name+"/"+filelist_temp+"{0:06d}.vtu".format(m+countstart_temp[n]))
            
            sys.exit()
        elif not(continue_from_break):
            if runTimeList is None:
                if manualPhaseTimeInput and self.manualPhaseTime is not None:
                    heart.set_activeTime(self.manualPhaseTime(0.))
                else:
                    heart.set_activeTime(0.)
            else:
                t=runTimeList[0]
                if manualPhaseTimeInput and self.manualPhaseTime is not None:
                    heart.set_activeTime(self.manualPhaseTime(runTimeList[0]))
                else:
                    heart.set_activeTime(runTimeList[0])
            self.solveLoadActive(heart)
            self.solveVolume(heart,startVolume)
            LVcav_array = [heart.get_CavityVolume()]
            Pcav_array = [heart.get_CavityPressure() ]
        else:
            with open(self.casename+"/"+str(self.runCount)+"/PV_.txt",'r') as f:
                PVfile=f.readlines()
            prev_PV=np.array(PVfile[0].replace('[','').replace(']','').split()).astype(float).reshape((1,-1))
            for n in range(1,len(PVfile)):
                prev_PV=np.concatenate((prev_PV,np.array(PVfile[n].replace('[','').replace(']','').split()).astype(float).reshape((1,-1))),axis=0)
            last_time=float(prev_PV[-1,0])
            num_of_chamber=int((prev_PV.shape[1]-1)/2)
            t_array=list(prev_PV[:,0])
            Pcav_array=list(prev_PV[:,1:(1+num_of_chamber)])
            LVcav_array=list(prev_PV[:,(1+num_of_chamber):(2*(1+num_of_chamber))])
            t=t_array[-1]
            if manualPhaseTimeInput and self.manualPhaseTime is not None:
                t2=self.manualPhaseTime(t)
                heart.set_activeTime(t2)
                if(t2 >= 0.0 and t2 < 4.0):
                    heart.set_dTime(0.50)
                elif (t2 >= (runParameters['time to maximum LV fiber tension in ms']-0.5) and t2 < (runParameters['time to maximum LV fiber tension in ms']+0.5)):
                    heart.set_dTime(3)
                else :
                    heart.set_dTime(1.0)
            else:
                heart.set_dTime(dt)#runTimeList[0]-t
                heart.set_activeTime(t)  
            if(t >= 0.0 and t < 4.0):
                dt = 0.50
            elif (t >= (runParameters['time to maximum LV fiber tension in ms']-0.5) and t < (runParameters['time to maximum LV fiber tension in ms']+0.5)):
                dt = 3
            else :
                dt = 1.0
            logger.info("t_a="+repr(t))
            V_cav=LVcav_array[-1]
            self.solveLoadActive(heart)
            self.solveVolume(heart,startVolume)
        #Joy Changed from here
        #Stress calculation 
        #if os.path.isfile(self.casename+"/"+str(self.runCount)+"/PV_.txt") and not(continue_from_break):
        #    os.remove(self.casename+"/"+str(self.runCount)+"/PV_.txt") 
        ###stress_File << cauchy ## wei xuan remove
        if(fenics.MPI.rank(heart.get_Comm()) == 0): #added to output the volume, stroke volume data before the while loop
            fdataPV =self.casename+"/"+str(self.runCount)+"/PV_.txt"
            if continue_from_break:
                fdatals = open(self.casename+"/"+str(self.runCount)+"/ls.txt", "a")
                fdataSV = open(self.casename+"/"+str(self.runCount)+"/SVp.txt", "a")
                fdataSVt = open(self.casename+"/"+str(self.runCount)+"/SVpt.txt", "a")
                fdata_stress = open( self.casename+"/"+str(self.runCount)+"/_stress_.txt", "a")
                fdataWork = open(self.casename+"/"+str(self.runCount)+"/work.txt", "a")
            else:
                fdatals = open(self.casename+"/"+str(self.runCount)+"/ls.txt", "w")
                fdataSV = open(self.casename+"/"+str(self.runCount)+"/SVp.txt", "w")
                fdataSVt = open(self.casename+"/"+str(self.runCount)+"/SVpt.txt", "w")
                fdata_stress = open( self.casename+"/"+str(self.runCount)+"/_stress_.txt", "w")
                fdataWork = open(self.casename+"/"+str(self.runCount)+"/work.txt", "w")
        if not(continue_from_break):
            p_cav = heart.get_CavityPressure()
            V_cav = heart.get_CavityVolume()
            if(fenics.MPI.rank(heart.get_Comm()) == 0):
            	logger.info("Cycle number = "+repr(cycle)+ " cell time = "+repr(t)+ " tstep = "+repr(tstep)+" dt = "+repr(heart.get_dTime()))
            	p_cav_str=np.array2string(p_cav)
            	if p_cav_str[0]=='[':
                    p_cav_str=p_cav_str[1:-1]
            	V_cav_str=np.array2string(V_cav)
            	if V_cav_str[0]=='[':
                    V_cav_str=V_cav_str[1:-1]
            	with open(fdataPV,'w') as f:
                    f.write(str(tstep)+' '+ p_cav_str  +' '+ V_cav_str)
        
            #if(fenics.MPI.rank(heart.get_Comm()) == 0):
            #    displacementfile << heart.get_DisplacementResult()
            #Joy changed  from here
            #cauchy1 =  heart.uflforms.Cauchy1() + heart.activeforms.cauchy()
            #cauchy = fenics.project(cauchy1,fenics.TensorFunctionSpace(heart.get_Mesh(), "DG", 1), form_compiler_parameters=project_form_compiler_parameters)
            #cauchy.rename("cauchy_stress","cauchy_stress")
            #stress_File << cauchy
            
            resultWriter(heart)
            sigma_fiber_LV = heart.get_VolumeAverageFiberStress()
       
       	#Joy chnaged till here
        #save files
       
        # Closed-loop phase
        if self.manualVol is None:
            ngspice_py.createLVcircuit(self.casename+'/'+self.meshname,runParameters)
        if runTimeList is not None:
            runTimeList=list(runTimeList)
            runTimeList.pop(0)
            if continue_from_break:
                while len(runTimeList)>0 and runTimeList[0]<=last_time:
                    runTimeList.pop(0)
        if endTime is None:
            endTime = BCL*self.numofCycle
        if electroSim:
            electroPHeart.replace_State_object(heart.get_State_object())
            electroPHeart.reset_Phase()
        if continue_from_break:
            tstep = last_time
        while(cycle < self.numofCycle):
        	
        	#################################################################
            p_cav = heart.get_CavityPressure()
            V_cav = heart.get_CavityVolume()
            
        	
            if runTimeList is None:
                tstep = tstep + dt
            elif len(runTimeList)<=0:
                break
            else:
                tstep =runTimeList.pop(0)
            if tstep>endTime:
                break
            cycle = math.floor(tstep/BCL)
            logger.info("cycle="+repr(cycle))
            t = tstep - cycle*BCL
            #these lines added as between the two timesteps if there is a large volume it will crash, therefore below increases volume gradually
            if manualPhaseTimeInput and self.manualPhaseTime is not None:
                t2=self.manualPhaseTime(t)
                heart.set_activeTime(t2)
                if(t2 >= 0.0 and t2 < 4.0):
                    heart.set_dTime(0.50)
                elif (t2 >= (runParameters['time to maximum LV fiber tension in ms']-0.5) and t2 < (runParameters['time to maximum LV fiber tension in ms']+0.5)):
                    heart.set_dTime(3)
                else :
                    heart.set_dTime(1.0)
            else:
                heart.set_dTime(dt)#runTimeList[0]-t
                heart.set_activeTime(t)   
            if(t >= 0.0 and t < 4.0):
                dt = 0.50
            elif (t >= (runParameters['time to maximum LV fiber tension in ms']-0.5) and t < (runParameters['time to maximum LV fiber tension in ms']+0.5)):
                dt = 3
            else :
                dt = 1.0
                
            logger.info("t_a="+repr(t))
        
        	
            if self.manualVol is not None:
                V_cav=self.manualVol(tstep)
            else:
                V_cav=self.runWindkessel(heart.get_Comm(),tstep,dt)
            t_array.append(t)
            LVcav_array.append(V_cav)
            Pcav_array.append(p_cav )
            
            self.solveVolume(heart,V_cav)
            #self.solveTa(heart,t,peaktime=runParameters['time to maximum LV fiber tension in ms'])
            if electroSim:
                if t<dt:
                    electroPHeart.reset_Phase()
                electroPHeart.solve()
                heart.update_Phase()
                electroPHeart.update_Phase()
                potential_ref=electroPHeart.get_interpolate_potential_ep2me_phi(heart.get_Mesh())
                heart.update_activationTime(potential_ref)
            p_cav = heart.get_CavityPressure()
            V_cav = heart.get_CavityVolume()
            
            if(fenics.MPI.rank(heart.get_Comm()) == 0):
                logger.info("Cycle number = "+repr(cycle)+ " cell time = "+repr(t)+ " tstep = "+repr(tstep)+" dt = "+repr(heart.get_dTime()))
                p_cav_str=np.array2string(p_cav)
                if p_cav_str[0]=='[':
                    p_cav_str=p_cav_str[1:-1]
                V_cav_str=np.array2string(V_cav)
                if V_cav_str[0]=='[':
                    V_cav_str=V_cav_str[1:-1]
                with open(fdataPV,'a') as f:
                    f.write('\n'+str(tstep)+' '+ p_cav_str  +' '+ V_cav_str)
            ls0 = runParameters['LV sarcomere length threshold where no tension develops in um']
            ls = fenics.sqrt(fenics.dot(heart.get_FiberDirectionVector(), heart.get_RightCauchyGreenDeformationTensorFunction()*heart.get_FiberDirectionVector()))*ls0
            ls1 = fenics.project(ls,heart.get_ScalarFunctionSpace(),form_compiler_parameters=project_form_compiler_parameters).vector().get_local()[:]
            eca = fenics.project(heart.get_ActiveCalciumStrain(), heart.get_ScalarFunctionSpace(),form_compiler_parameters=project_form_compiler_parameters).vector().get_local()[:]
            t_r = fenics.project(heart.get_RelaxationTimeLength(), heart.get_ScalarFunctionSpace(),form_compiler_parameters=project_form_compiler_parameters).vector().get_local()[:]
            if(fenics.MPI.rank(heart.get_Comm()) == 0):
                print(heart.get_activeTime(), min(ls1), max(ls1), min(eca), max(eca), min(t_r), max(t_r), file=fdatals)
                resultWriter(heart)
            sigma_fiber_LV = heart.get_VolumeAverageFiberStress()
            
            if(fenics.MPI.rank(heart.get_Comm()) == 0):
                print (fdata_stress, tstep, sigma_fiber_LV, file=fdata_stress )
                print(tstep, heart.get_TotalStrainEnergy() , file=fdataWork)
        	#Joy changed Till here         
            
            if trackphase:
                heart.update_Phase()
        

        if not(os.path.isfile(self.casename+"/"+str(self.runCount)+"/PV_.txt")):
            np.savetxt(self.casename+"/"+str(self.runCount)+"/PV_.txt",np.array([t_array,Pcav_array,LVcav_array]))
        if (tstep >= runParameters['duration of one cardiac cycle in ms']):
        	t_array=np.array(t_array)
        	if (tstep >= 2*runParameters['duration of one cardiac cycle in ms']):
        		temp_Ind=np.nonzero((t_array[2:]-t_array[:-2])>0)[0][-1]
        		temp_ind=np.nonzero(np.logical_and(t_array[:temp_Ind]>t_array[temp_Ind],t_array[:temp_Ind]<t_array[temp_Ind+2]))[0][-1]
    
        		LVcav_arrayT=LVcav_array[(temp_ind-1):]
        		Pcav_arrayT=Pcav_array[(temp_ind-1):]
        	else:
        		LVcav_arrayT=LVcav_array
        		Pcav_arrayT=Pcav_array
        	plt.plot(LVcav_arrayT, Pcav_arrayT)
        	plt.xlabel("Volume")
        	plt.ylabel("Pressure")
        	plt.legend()
        	plt.savefig(self.casename+"/PV.png")
        	
        	print(np.nanmax(LVcav_arrayT),np.nanmin(LVcav_arrayT),np.nanmax(LVcav_arrayT)-np.nanmin(LVcav_arrayT), np.nanmax(Pcav_arrayT),np.nanmin(Pcav_arrayT),np.nanmax(Pcav_arrayT)-np.nanmin(Pcav_arrayT), file=fdataSV) 
        	print(LVcav_arrayT,Pcav_arrayT, file=fdataSVt)
        
        if(fenics.MPI.rank(heart.get_Comm()) == 0):
        	fdataSV.close()
        	fdataSVt.close()
        	fdata_stress.close() #Joy Changed here	
        	fdataWork.close()
        ######################################################################################################
        return runParameters
class fenicsResultWriter:
    def __new__(cls, *args, **kwargs):
        return super().__new__(cls)
    def __init__(self,savePath,outputResultList,project_form_compiler_parameters=None,tag=None):
        self.savePath=savePath
        self.displacementFile=None
        self.deformationgradienttensorFile=None
        self.stressFile=None
        self.fiberStressFile=None
        self.stresstraceFile=None
        self.stressvonmisesFile=None
        self.longiStrainFile=None
        self.fiberStrainFile=None
        self.workFile=None
        self.fiberFile=None
        self.lcrstretchFile=None
        self.activationdelayFile=None
        if project_form_compiler_parameters is None:
            self.project_form_compiler_parameters={"representation":"uflacs"}
        else:
            self.project_form_compiler_parameters=project_form_compiler_parameters
        if tag is not None:
            self.tag=[np.zeros((0,len(tag))),tag]
        else:
            self.tag=None
        if "displacement" in outputResultList or "deformation" in outputResultList:
            os.makedirs(self.savePath+"/deformation",exist_ok=True)
            self.displacementFile = fenics.File(self.savePath+"/deformation/u_disp.pvd")
        if "deformationgradienttensor" in outputResultList:
            os.makedirs(self.savePath+"/deformationgradienttensor",exist_ok=True)
            self.deformationgradienttensorFile = fenics.File(self.savePath+"/deformationgradienttensor/_deformationgradienttensor.pvd") #Joy Changed here
        if "stress" in outputResultList:
            os.makedirs(self.savePath+"/stress",exist_ok=True)
            self.stressFile = fenics.File(self.savePath+"/stress/_stress.pvd") #Joy Changed here
        if "fiberstress" in outputResultList:
            os.makedirs(self.savePath+"/fiberstress",exist_ok=True)
            self.fiberStressFile = []
            for n in ['ff','ss','nn','fs','fn','ns']:
                self.fiberStressFile.append(fenics.File(self.savePath+"/fiberstress/_fiberstress_"+n+".pvd")) #Joy Changed here
        if "stresstrace" in outputResultList:
            os.makedirs(self.savePath+"/stresstrace",exist_ok=True)
            self.stresstraceFile = fenics.File(self.savePath+"/stresstrace/_stresstrace.pvd") #mf added
        if "stressvonmises" in outputResultList:
            os.makedirs(self.savePath+"/stressvonmises",exist_ok=True)
            self.stressvonmisesFile = fenics.File(self.savePath+"/stressvonmises/_stressvonmises.pvd") #mf added
        if "longistrain" in outputResultList:
            os.makedirs(self.savePath+"/strain",exist_ok=True)
            self.longiStrainFile = []
            for n in ['ll','cc','rr','lc','lr','rc']:
                self.longiStrainFile.append(fenics.File(self.savePath+"/strain/_longistrain_"+n+".pvd"))
        if "fiberstrain" in outputResultList:
            os.makedirs(self.savePath+"/strain",exist_ok=True)
            self.fiberStrainFile = []
            for n in ['ff','ss','nn','fs','fn','ns']:
                self.fiberStrainFile.append(fenics.File(self.savePath+"/strain/_fiberstrain_"+n+".pvd"))
        if "strain_energy_density" in outputResultList:
            os.makedirs(self.savePath+"/strain_energy_density",exist_ok=True)
            self.workFile = fenics.File(self.savePath+"/strain_energy_density/_strain_energy_density.pvd")
        if "fiber_stretch" in outputResultList:
            os.makedirs(self.savePath+"/fiber_stretch",exist_ok=True)
            self.fiberFile=[]
            for n in ['ff','ss','nn']:
                self.fiberFile.append(fenics.File(self.savePath+"/fiber_stretch/_fiber_"+n+".pvd"))
        if "lcr_stretch" in outputResultList:
            os.makedirs(self.savePath+"/lcr_stretch",exist_ok=True)
            self.lcrstretchFile=[]
            for n in ['longitudinal','circumferential','radial']:
                self.lcrstretchFile.append(fenics.File(self.savePath+"/lcr_stretch/_"+n+".pvd"))
        if "activation_delay" in outputResultList:
            os.makedirs(self.savePath+"/electro",exist_ok=True)
            self.activationdelayFile = fenics.File(self.savePath+"/electro/activation_delay.pvd")
    def __call__(self,heart,tag=None):
        if(fenics.MPI.rank(heart.get_Comm()) == 0):
            if self.displacementFile is not None:
                self.displacementFile << heart.get_DisplacementResult()
            if self.deformationgradienttensorFile is not None:
                Fmat1=heart.get_StretchTensorFunction()
                Fmat=fenics.project(Fmat1,fenics.TensorFunctionSpace(heart.get_Mesh(), "DG", 1), form_compiler_parameters=self.project_form_compiler_parameters)
                Fmat.rename("deformation_gradient_tensor","deformation_gradient_tensor")
                self.deformationgradienttensorFile << Fmat
            if self.stressFile is not None:
                cauchy = heart.get_CauchyFiberStressTensor()
                cauchy.rename("cauchy_stress","cauchy_stress")
                self.stressFile << cauchy
            if self.fiberStressFile is not None:
                cauchy = heart.get_CauchyFiberStressTensor()
                f0 = heart.get_FiberDirectionVector()
                s0 = heart.get_FiberSheetletDirectionVector()
                n0 = heart.get_FiberNormalDirectionVector()
                cauchyff = fenics.project(fenics.inner(f0, cauchy*f0),fenics.FunctionSpace(heart.get_Mesh(), "CG", 1), form_compiler_parameters=self.project_form_compiler_parameters)
                cauchyff.rename("ff_stress","ff_stress")
                self.fiberStressFile[0] << cauchyff
                cauchyss = fenics.project(fenics.inner(s0, cauchy*s0),fenics.FunctionSpace(heart.get_Mesh(), "CG", 1), form_compiler_parameters=self.project_form_compiler_parameters)
                cauchyss.rename("ss_stress","ss_stress")
                self.fiberStressFile[1] << cauchyss
                cauchynn = fenics.project(fenics.inner(n0, cauchy*n0),fenics.FunctionSpace(heart.get_Mesh(), "CG", 1), form_compiler_parameters=self.project_form_compiler_parameters)
                cauchynn.rename("nn_stress","nn_stress")
                self.fiberStressFile[2] << cauchynn
                cauchyfs = fenics.project(fenics.inner(f0, cauchy*s0),fenics.FunctionSpace(heart.get_Mesh(), "CG", 1), form_compiler_parameters=self.project_form_compiler_parameters)
                cauchyfs.rename("fs_stress","fs_stress")
                self.fiberStressFile[3] << cauchyfs
                cauchyfn = fenics.project(fenics.inner(f0, cauchy*n0),fenics.FunctionSpace(heart.get_Mesh(), "CG", 1), form_compiler_parameters=self.project_form_compiler_parameters)
                cauchyfn.rename("fn_stress","fn_stress")
                self.fiberStressFile[4] << cauchyfn
                cauchyns = fenics.project(fenics.inner(n0, cauchy*s0),fenics.FunctionSpace(heart.get_Mesh(), "CG", 1), form_compiler_parameters=self.project_form_compiler_parameters)
                cauchyns.rename("ns_stress","ns_stress")
                self.fiberStressFile[5] << cauchyns
            if self.stresstraceFile is not None:#mf added
                cauchy = heart.get_CauchyFiberStressTensor()
                trace= fenics.project(fenics.tr(cauchy),fenics.FunctionSpace(heart.get_Mesh(), "CG", 1), form_compiler_parameters=self.project_form_compiler_parameters)
                trace.rename("stresstrace","stresstrace")
                self.stresstraceFile << trace
            if self.stressvonmisesFile is not None:#mf added
                cauchy = heart.get_CauchyFiberStressTensor()
                stress_von_mises =fenics.sqrt(0.5*((cauchy[0,0]-cauchy[1,1])**2+(cauchy[1,1]-cauchy[2,2])**2+(cauchy[0,0]-cauchy[2,2])**2)+3*(cauchy[0,1]**2+cauchy[1,2]**2+cauchy[2,0]**2))
                von_mises= fenics.project(stress_von_mises,fenics.FunctionSpace(heart.get_Mesh(), "CG", 1), form_compiler_parameters=self.project_form_compiler_parameters)
                von_mises.rename("stressvonmises","stressvonmises")
                self.stressvonmisesFile << von_mises
            if self.longiStrainFile is not None:
                Ea=heart.get_StrainTensorFunction()
                long0 = heart.get_LongitudinalDirectionVector()
                cir0 = heart.get_CircumferentialDirectionVector()
                rad0 = heart.get_RadialDirectionVector()
                Ell = fenics.project(fenics.inner(long0, Ea*long0),fenics.FunctionSpace(heart.get_Mesh(), "CG", 1), form_compiler_parameters=self.project_form_compiler_parameters)
                Ell.rename("Ell_strain","Ell_strain")
                self.longiStrainFile[0] << Ell
                Ecc = fenics.project(fenics.inner(cir0, Ea*cir0),fenics.FunctionSpace(heart.get_Mesh(), "CG", 1), form_compiler_parameters=self.project_form_compiler_parameters)
                Ecc.rename("Ecc_strain","Ecc_strain")
                self.longiStrainFile[1] << Ecc
                Err = fenics.project(fenics.inner(rad0, Ea*rad0),fenics.FunctionSpace(heart.get_Mesh(), "CG", 1), form_compiler_parameters=self.project_form_compiler_parameters)
                Err.rename("Err_strain","Err_strain")
                self.longiStrainFile[2] << Err
                Elc = fenics.project(fenics.inner(long0, Ea*cir0),fenics.FunctionSpace(heart.get_Mesh(), "CG", 1), form_compiler_parameters=self.project_form_compiler_parameters)
                Elc.rename("Elc_strain","Elc_strain")
                self.longiStrainFile[3] << Elc
                Elr = fenics.project(fenics.inner(long0, Ea*rad0),fenics.FunctionSpace(heart.get_Mesh(), "CG", 1), form_compiler_parameters=self.project_form_compiler_parameters)
                Elr.rename("Elr_strain","Elr_strain")
                self.longiStrainFile[4] << Elr
                Erc = fenics.project(fenics.inner(rad0, Ea*cir0),fenics.FunctionSpace(heart.get_Mesh(), "CG", 1), form_compiler_parameters=self.project_form_compiler_parameters)
                Erc.rename("Erc_strain","Erc_strain")
                self.longiStrainFile[5] << Erc
            if self.fiberStrainFile is not None:
                Ea=heart.get_StrainTensorFunction()
                f0 = heart.get_FiberDirectionVector()
                s0 = heart.get_FiberSheetletDirectionVector()
                n0 = heart.get_FiberNormalDirectionVector()
                Eff = fenics.project(fenics.inner(f0, Ea*f0),fenics.FunctionSpace(heart.get_Mesh(), "CG", 1), form_compiler_parameters=self.project_form_compiler_parameters)
                Eff.rename("Eff_strain","Eff_strain")
                self.fiberStrainFile[0] << Eff
                Ess = fenics.project(fenics.inner(s0, Ea*s0),fenics.FunctionSpace(heart.get_Mesh(), "CG", 1), form_compiler_parameters=self.project_form_compiler_parameters)
                Ess.rename("Ess_strain","Ess_strain")
                self.fiberStrainFile[1] << Ess
                Enn = fenics.project(fenics.inner(n0, Ea*n0),fenics.FunctionSpace(heart.get_Mesh(), "CG", 1), form_compiler_parameters=self.project_form_compiler_parameters)
                Enn.rename("Enn_strain","Enn_strain")
                self.fiberStrainFile[2] << Enn
                Efs = fenics.project(fenics.inner(f0, Ea*s0),fenics.FunctionSpace(heart.get_Mesh(), "CG", 1), form_compiler_parameters=self.project_form_compiler_parameters)
                Efs.rename("Efs_strain","Efs_strain")
                self.fiberStrainFile[3] << Efs
                Efn = fenics.project(fenics.inner(f0, Ea*n0),fenics.FunctionSpace(heart.get_Mesh(), "CG", 1), form_compiler_parameters=self.project_form_compiler_parameters)
                Efn.rename("Efn_strain","Efn_strain")
                self.fiberStrainFile[4] << Efn
                Ens = fenics.project(fenics.inner(n0, Ea*s0),fenics.FunctionSpace(heart.get_Mesh(), "CG", 1), form_compiler_parameters=self.project_form_compiler_parameters)
                Ens.rename("Ens_strain","Ens_strain")
                self.fiberStrainFile[5] << Ens
            if self.workFile is not None:
                work=heart.get_StrainEnergyDensity()
                work.rename("strain_energy_density","strain_energy_density")
                self.workFile << work
            if self.fiberFile is not None:
                Fmat1=heart.get_StretchTensorFunction()
                f0 = heart.get_FiberDirectionVector()
                s0 = heart.get_FiberSheetletDirectionVector()
                n0 = heart.get_FiberNormalDirectionVector()
                fiber_ff = fenics.project(Fmat1*f0,fenics.VectorFunctionSpace(heart.get_Mesh(), "CG", 1), form_compiler_parameters=self.project_form_compiler_parameters)
                fiber_ff.rename("fiber_ff","fiber_ff")
                self.fiberFile[0] << fiber_ff
                fiber_ss = fenics.project(Fmat1*s0,fenics.VectorFunctionSpace(heart.get_Mesh(), "CG", 1), form_compiler_parameters=self.project_form_compiler_parameters)
                fiber_ss.rename("fiber_ss","fiber_ss")
                self.fiberFile[1] << fiber_ss
                fiber_nn = fenics.project(Fmat1*n0,fenics.VectorFunctionSpace(heart.get_Mesh(), "CG", 1), form_compiler_parameters=self.project_form_compiler_parameters)
                fiber_nn.rename("fiber_nn","fiber_nn")
                self.fiberFile[2] << fiber_nn
            if self.lcrstretchFile is not None:
                Fmat1=heart.get_StretchTensorFunction()
                long0 = heart.get_LongitudinalDirectionVector()
                cir0 = heart.get_CircumferentialDirectionVector()
                rad0 = heart.get_RadialDirectionVector()
                stretch_ll = fenics.project(Fmat1*long0,fenics.VectorFunctionSpace(heart.get_Mesh(), "CG", 1), form_compiler_parameters=self.project_form_compiler_parameters)
                stretch_ll.rename("stretch_ll","stretch_ll")
                self.lcrstretchFile[0] << stretch_ll
                stretch_cc = fenics.project(Fmat1*cir0,fenics.VectorFunctionSpace(heart.get_Mesh(), "CG", 1), form_compiler_parameters=self.project_form_compiler_parameters)
                stretch_cc.rename("stretch_cc","stretch_cc")
                self.lcrstretchFile[1] << stretch_cc
                stretch_rr = fenics.project(Fmat1*rad0,fenics.VectorFunctionSpace(heart.get_Mesh(), "CG", 1), form_compiler_parameters=self.project_form_compiler_parameters)
                stretch_rr.rename("stretch_rr","stretch_rr")
                self.lcrstretchFile[2] << stretch_rr
            if self.activationdelayFile is not None:
                activationdelay_Function = heart.get_ActivationDelay()
                activationdelay = fenics.project(activationdelay_Function,fenics.VectorFunctionSpace(heart.get_Mesh(), "CG", 1), form_compiler_parameters=self.project_form_compiler_parameters)
                activationdelay.rename("activation_delay","activation_delay")
                self.activationdelayFile << activationdelay
            if self.tag is not None:
                temp_tag_val=np.array([[float('nan')]*len(self.tag[1])])
                if tag is not None:
                    for n in range(len(self.tag[1])):
                        if self.tag[1][n] in tag:
                            temp_tag_val[0,n]=tag[self.tag[1][n]]
                self.tag[0]=np.concatenate((self.tag[0],temp_tag_val),axis=0)
                np.savetxt(self.savePath+"/result_tags.txt",np.concatenate((np.arange(self.tag[0].shape[0]).reshape((-1,1)),self.tag[0]),axis=1),comments="#ID ",header=" ".join(self.tag[1]))
class tanh_limit:
    def __init__(self,low,high,inverse=False,log=0):
        self.log=log
        if log>0:
            self.low=np.log(low)/np.log(self.log)
            self.high=np.log(high)/np.log(self.log)
        else:
            self.low=low
            self.high=high
        self.inverse=inverse
    def __call__(self,y,inverse=None):
        if inverse is None:
            inverse=self.inverse
        if inverse:
            if self.log>0:
                x=np.log(y)/np.log(self.log)
            else:
                x=y
            if self.low==float('-inf') and self.high==float('inf'):
                return x
            elif self.high==float('inf'):
                return np.log(x-self.low)
            elif self.low==float('-inf'):
                return -np.log(self.high-x)
            else:
                if x<=self.low:
                    return float('-inf')
                elif x>=self.high:
                    return float('inf')
                else:
                    #return np.arctanh(min(1.,max(-1.,(x-self.low)*2./(self.high-self.low)-1.)))*(self.high-self.low)**2.
                    return -np.log((self.high-self.low)/(x-self.low)-1.)*(self.high-self.low)**2.
        else:
            if self.low==float('-inf') and self.high==float('inf'):
                result=y
            elif self.high==float('inf'):
                result=np.exp(y)+self.low
            elif self.low==float('-inf'):
                result=-np.exp(-y)+self.high
            else:
                #return (np.tanh(x/(self.high-self.low)**2.)+1.)/2.*(self.high-self.low)+self.low
                result=(self.high-self.low)/(1+np.exp(-y/(self.high-self.low)**2.))+self.low
            if self.log>0:
                return self.log**result
            else:
                return result
def plotLVBehavior(mainPath,cavity='LV',labelList=None,timeSlice=None,volumeSlice=None,fmt='png',show=False):
    if not(isinstance(mainPath,list)):
        mainPath=[mainPath]
    removeVolumeratio=0.1
    VolTime_total=[]
    timeSpace=[]
    volumeSpace=[]
    VolTime_base=[]
    VolTime=[]
    for n in range(len(mainPath)):
        timeSpace.append(np.loadtxt(mainPath[n]+'_Press_timeSpace.txt'))
        volumeSpace.append(np.loadtxt(mainPath[n]+'_Press_volumeSpace.txt'))
        VolTime_base.append(np.loadtxt(mainPath[n]+'_Press_VolTime_'+cavity+'_base.txt'))
        VolTime.append(np.loadtxt(mainPath[n]+'_Press_VolTime_'+cavity+'.txt'))
        
        remove_volumeSpace=int(len(volumeSpace[-1])*removeVolumeratio)
        remove_timeSpace=int(len(timeSpace[-1])*0.5)
        volumeSpace[-1]=volumeSpace[-1][remove_volumeSpace:-remove_volumeSpace]
        VolTime_base[-1]=VolTime_base[-1][remove_volumeSpace:-remove_volumeSpace]
        VolTime[-1]=VolTime[-1][remove_volumeSpace:-remove_volumeSpace]
        timeSpace[-1]=timeSpace[-1][:-remove_timeSpace]
        VolTime[-1]=VolTime[-1][:,:-remove_timeSpace]
        
        VolTime_total.append(VolTime_base[-1].reshape((-1,1))+VolTime[-1])
    import matplotlib.pyplot as plt
    from matplotlib import cm
    import matplotlib.ticker as mticker
    import matplotlib.colors as colors
    cList=['k','b','r','g','c','m','y']
    if timeSlice is not None:
        fig,ax=plt.subplots(1,1)
        for m in range(len(mainPath)):
            data=[]
            for n in range(len(volumeSpace[m])):
                spl=scipyinterpolate.splrep(timeSpace[m], VolTime_total[m][n,:])
                data.append(scipyinterpolate.splev(np.array([timeSlice]), spl)[0])
            if labelList is not None:
                plt.plot(volumeSpace[m],data,color=cList[m],label=labelList[m])
            else:
                plt.plot(volumeSpace[m],data,color=cList[m])
        ax.set_title('LV Behavior at {0:.1f}ms'.format(timeSlice))
        ax.set_xlabel('Volume (mL)')
        ax.set_ylabel('Pressure (mmHg)')
        if labelList is not None:
            ax.legend(loc=0)
        plt.savefig(mainPath[0]+'_LVBehaviorPlot_time_{0:.1f}ms.png'.format(timeSlice),dpi=fig.dpi,format=fmt,bbox_inches='tight')
        fig.clf()
    if volumeSlice is not None:
        fig,ax=plt.subplots(1,1)
        for m in range(len(mainPath)):
            data=[]
            for n in range(len(timeSpace[m])):
                spl=scipyinterpolate.splrep(volumeSpace[m][::-1], VolTime_total[m][::-1,n])
                data.append(scipyinterpolate.splev(np.array([volumeSlice]), spl)[0])
            if labelList is not None:
                plt.plot(timeSpace[m],data,color=cList[m],label=labelList[m])
            else:
                plt.plot(timeSpace[m],data,color=cList[m])
        ax.set_title('LV Behavior at {0:.2f}mL'.format(volumeSlice))
        ax.set_xlabel('Time (ms)')
        ax.set_ylabel('Pressure (mmHg)')
        if labelList is not None:
            ax.legend(loc=0)
        plt.savefig(mainPath[0]+'_LVBehaviorPlot_volume_{0:.1f}mL.png'.format(volumeSlice),dpi=fig.dpi,format=fmt,bbox_inches='tight')
        fig.clf()
        
    fig,ax=plt.subplots(1,1)
    for m in range(len(mainPath)):
        plt.plot(volumeSpace[m],VolTime_base[m],color=cList[m],label=mainPath[m])
    ax.set_title('passive LV Behavior - Pressure (mmHg)')
    ax.set_xlabel('Volume (mL)')
    ax.set_ylabel('Passive Pressure (mmHg)')
    if len(mainPath)>1:
        ax.legend(loc=0)
    plt.savefig(mainPath[0]+'_passiveLVBehaviorPlot.png'.format(volumeSlice),dpi=fig.dpi,format=fmt,bbox_inches='tight')
    fig.clf()
    for m in range(len(mainPath)):
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
        X, Y = np.meshgrid(timeSpace, volumeSpace)
        
        # Customize the z axis.
        vmax=VolTime[m].max()
        vmin=VolTime[m].min()
        surf = ax.plot_surface(X, Y, VolTime[m],cmap=cm.coolwarm,linewidth=0, antialiased=False)
    
        ax.zaxis.set_major_locator(mticker.LinearLocator(10))
        ax.zaxis.set_major_formatter('{x:.01f}')
        ax.set_zlim(vmin,vmax)
            
        
        # Add a color bar which maps values to colors.
        fig.colorbar(surf, shrink=0.5, aspect=5)
        ax.set_title('active LV Behavior - Pressure (mmHg)')
        ax.set_xlabel('Time (ms)')
        ax.set_ylabel('Volume (mL)')
        ax.set_zlabel('Pressure (mmHg)')
        plt.savefig(mainPath[0]+'_activeLVBehaviorPlot.png',dpi=fig.dpi,format=fmt,bbox_inches='tight')
        fig.clf()
        
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
        X, Y = np.meshgrid(timeSpace, volumeSpace)
        
        # Customize the z axis.
        vmax=VolTime_total[m].max()
        vmin=VolTime_total[m].min()
        surf = ax.plot_surface(X, Y, VolTime_total[m],cmap=cm.coolwarm,linewidth=0, antialiased=False)
    
        ax.zaxis.set_major_locator(mticker.LinearLocator(10))
        ax.zaxis.set_major_formatter('{x:.01f}')
        ax.set_zlim(vmin,vmax)
            
        
        # Add a color bar which maps values to colors.
        fig.colorbar(surf, shrink=0.5, aspect=5)
        ax.set_title('LV Behavior - Pressure (mmHg)')
        ax.set_xlabel('Time (ms)')
        ax.set_ylabel('Volume (mL)')
        ax.set_zlabel('Pressure (mmHg)')
        plt.savefig(mainPath[0]+'_LVBehaviorPlot.png',dpi=fig.dpi,format=fmt,bbox_inches='tight')
        if show:
            plt.show()
        fig.clf()

def exactFunc(x,inverse=None):
    return x
class optimiser_linker:
    def __new__(cls, *args, **kwargs):
        return super().__new__(cls)
    def __init__(self,runSimulationClass,calculateCostClass):
        self.runSimulation=runSimulationClass
        self.calculateCost=calculateCostClass
        self.variableList=[]
        self.variableOperatorList=[]
        self.minCost=float('inf')
        self.kwargs={}
    def addVariable(self,variables,variableOperator=None,log=0):
        if isinstance(variables,str):
            variables=[variables]
        if variableOperator is None:
            variableOperator=[exactFunc]*len(variables)
        elif not(isinstance(variableOperator,(list,tuple))):
            variableOperator=[variableOperator]*len(variables)
        elif isinstance(variableOperator[0],(float,int)):
                variableOperator=[variableOperator]*len(variables)
        if not(isinstance(log,(list,tuple))):
            logList=[log]*len(variables)
        else:
            logList=log
        self.variableList+=variables
        for n in range(len(variables)):
            if variableOperator[n] is None:
                self.variableOperatorList+=[exactFunc]
            if isinstance(variableOperator[n],(list,tuple)):
                self.variableOperatorList+=[tanh_limit(variableOperator[n][0],variableOperator[n][1],log=logList[n])]
            else:
                self.variableOperatorList+=[variableOperator[n]]
    def invOperate_para(self,para):
        result=[]
        for n in range(len(para)):
            result.append(self.variableOperatorList[n](para[n],inverse=True))
        return np.array(result)
    def set_kwargs(self,**kwargs):
        if len(kwargs)>0:
            self.kwargs=kwargs
    def __call__(self,para):
        if len(self.variableList)==0:
            raise Exception('Select variable to optimise with optimiser.addVariable(VARIABLE_STRING[s])')
        if isinstance(self.runSimulation.runCount,str):
            try:
                removecounter=self.runSimulation.runCount[::-1].index('/')
                subFolderName=self.runSimulation.runCount[:-removecounter]
            except:
                subFolderName=''
        os.makedirs(self.runSimulation.casename+'/'+subFolderName+'best_fit',exist_ok=True)
        editParameters={}
        for n in range(len(self.variableList)):
            editParameters[self.variableList[n]]=self.variableOperatorList[n](para[n])
        try:
            if isinstance(self.runSimulation.runCount,int):
                self.runSimulation.runCount+=1
            elif isinstance(self.runSimulation.runCount,str):
                if re.match('.*?([0-9]+)$', self.runSimulation.runCount)== None:
                    last_digits = None
                else:
                    last_digits = re.match('.*?([0-9]+)$', self.runSimulation.runCount).group(1)
                if last_digits is None:
                    self.runSimulation.runCount+='0'
                else:
                    self.runSimulation.runCount=self.runSimulation.runCount[:-len(last_digits)]+str(int(last_digits)+1)
            kwargs=self.kwargs
            runParameters=self.runSimulation(editParameters=editParameters,**kwargs)  
        except Exception as e_inst:
            print(repr(e_inst))
            with open(self.runSimulation.casename+'/'+str(self.runSimulation.runCount)+"/runParameters.txt", "r") as f:
                 old = f.readlines() # read everything in the file
            old.insert(0,"#Status , FAILED\n# cost , "+str(type(e_inst))+'  '+str(type(e_inst))+"\n")
            with open(self.runSimulation.casename+'/'+str(self.runSimulation.runCount)+"/runParameters.txt", "w") as f:
                 f.writelines(old) # write the new line before
            logger.warning("FAILED RUN")
            return float('inf')
        else:
            component_Cost=self.calculateCost(self.runSimulation.casename+'/'+str(self.runSimulation.runCount),runParameters)
            cost_current=np.sum(component_Cost)
            addline=''
            for n in range(len(self.calculateCost.costs)):
                addline+=' , '+self.calculateCost.costs[n]+'='+str(component_Cost[n])+'/'+str(self.calculateCost.weights[n])
            with open(self.runSimulation.casename+'/'+str(self.runSimulation.runCount)+"/runParameters.txt", "r") as f:
                 old = f.readlines() # read everything in the file
            old.insert(0,"#Status , SUCCESS\n# cost , "+str(cost_current)+addline+"\n")
            with open(self.runSimulation.casename+'/'+str(self.runSimulation.runCount)+"/runParameters.txt", "w") as f:
                 f.writelines(old) # write the new line before
                 
            if cost_current<self.minCost:
                self.minCost=cost_current
                
                    
                cmd = "rm -r " + self.runSimulation.casename+'/'+subFolderName+'best_fit'
                os.system(cmd)
                os.makedirs(self.runSimulation.casename+'/'+subFolderName+'best_fit',exist_ok=True)
                cmd = "cp -r " + self.runSimulation.casename+'/'+str(self.runSimulation.runCount)+'/* '+ self.runSimulation.casename+'/'+subFolderName+'best_fit'
                os.system(cmd)
            logger.info("cost="+repr(cost_current))
            return cost_current
def MeanSquare(x,ref):
    result=(x-ref)**2.
    if isinstance(result,np.ndarray):
        result.mean()
    return result
def Exact(x,ref):
    if isinstance(x,np.ndarray):
        result=x.mean()
    else:
        result=x
    return result
def getCircuitResult(folder,ngspice_varname,last_cycle=False):
    filename='/circuit_results'
    if last_cycle:
        filename=filename+'_lastcycle'
    filename=filename+'.txt'
    with open(folder+filename,'r') as f:
        header=f.readline()
    if header[-1]=='\n':
        header=header[:-1]
    headerList=header.split()
    while not('time' in headerList[0]):
        headerList.pop(0)
    if last_cycle:
        result=np.loadtxt(folder+filename)
    else:
        result=np.loadtxt(folder+filename,skiprows=1)
    if result.shape[1]!=len(headerList):
        logger.warning("Circuit results shape "+repr(result.shape)+' does not match header list with length '+str(len(headerList)))
    return result[:,headerList.index(ngspice_varname)]
        
class cost_function:
    def __new__(cls, *args, **kwargs):
        return super().__new__(cls)
    def __init__(self,runMode,matchingFolder='ref'):
        self.runMode=runMode
        self.costs=[]
        self.weights=[]
        if os.path.isdir(matchingFolder):
            self.matchingFolder=matchingFolder
        else:
            self.matchingFolder=os.getcwd()+'/'+matchingFolder
        self.collateNameSwitch={'MeanSquare':MeanSquare,
                                'Exact':Exact
                            }
        self.varNameSwitch={'LeftVentriclePressure':'p_cav',
                         'LeftVentricleVolume':'V_cav',
                         'StrokeVolume':'V_sv',
                         'MaxLVLAPressure':'pmax_LVLA',
                         'MaxLAVALVPressure':'pmax_LAVALVregurge',
                         'MaxLVAAPressure':'pmax_LVAA',
                         'MaxLVVALVPressure':'pmax_LVVALV',
                         'MinLVVALVPressure':'pmin_LVVALV',
                         'IsovolumetricContractionTime':'t_isocontract',
                         'IsovolumetricRelaxationTime':'t_isorelax',
                         't_la':'t_la', 
                         'VeinPressure':'Pven',
                         'ArterialPressure':'Part', 
                         'LeftAtriumPressure':'PLA', 
                         'Qla':'Qla', 
                         'AortaQ':'Qao', 
                         'perQ':'Qper', 
                         'mvQ':'Qmv', 
                         'VeinVolume':'V_ven', 
                         'ArterialVolume':'V_art', 
                         'LeftAtriumVolume':'V_LA',
                         'RightVentriclePressure':'PRV',
                         'RightAtriumPressure':'PRA',
                         'PulmonaryArteryPressure':'Ppa1',
                         'AscendingAortaPressure':'Paa',
                         'DiastolicMeshCoordinates':'coordsDiaMesh',
                         'SystolicMeshCoordinates':'coordsSysMesh'                        
                         }
    def __call__(self,folderName,runParameters):
        f=[]
        for n in range(len(self.costs)):
            func,var=self.decodeCost(self.costs[n])
            f.append(self.getCost(folderName,runParameters,var,func)*self.weights[n])
        return f
    def decodeCost(self,costLabel):
        for key in self.collateNameSwitch:
            if len(costLabel)>len(key):
                if costLabel[-len(key):]==key:
                    return (self.collateNameSwitch[key],self.varNameSwitch[costLabel[:-len(key)]])
        logger.warning("No Cost Found for "+costLabel)
        logger.warning("Collating Cost with "+repr(self.collateNameSwitch.keys()))
        logger.warning("Variable avaliable "+repr(self.varNameSwitch.keys()))
    def addCost(self,costString,weights=None):
        #put weight to negative to cost which is better if higher
        if weights is None:
            weights=1.
        if isinstance(costString,str):
            costString=[costString]
        if isinstance(weights,(float,int)):
            weights=np.ones(len(costString))*weights
        for n in range(len(costString)):
            self.costs.append(costString[n])
            self.weights.append(weights[n])
    def removeCost(self,costString):
        if isinstance(costString,str):
            costString=[costString]
        for string in costString:
            if string in self.costs:
                temp_ind=self.costs.index(string)
                self.costs.pop(temp_ind)
                self.weights.pop(temp_ind)
    def getspline(self,folderName,varName,returnArray=False):
        if returnArray:
            def returnFunc(*xx):
                return np.array(xx).T
        else:
            returnFunc=scipyinterpolate.splrep
        if self.runMode=='LVbehaviorRun' or self.runMode=='fullWindkesselRun' or self.runMode=='flowrateWindkesselRun':
            fdatacl = np.loadtxt(folderName+"/circuit_results.txt",skiprows=1)
            name=[' ', 'p_cav',' ', 'V_cav',' ', 'PLA',' ','PRV',' ','PRA', ' ','Ppa1',' ','Paa' ,' ','Plavalv' ,' ','Plvvalv' ,' ','QLV_in' , ' ','PhaseTime' , ' ','V_RV' , ' ','QLALV' , ' ','QLVAA' , ' ','QLVLA']
            return returnFunc(fdatacl[:,0], fdatacl[:,name.index(varName)])
        elif self.runMode=='iterativeRun':
            try:
                fdatacl = np.loadtxt(folderName+"/circuit_results.txt",skiprows=1)
                name=[' ', 'p_cav',' ', 'V_cav',' ', 'PLA',' ','PRV',' ','PRA', ' ','Ppa1',' ','Paa',' ','Plvvalv' ]
                return returnFunc(fdatacl[:,0], fdatacl[:,name.index(varName)])
            except:
                try:
                    fdatacl = np.loadtxt(folderName+"/cl.txt")
                    name=['tstep', 'p_cav', 'V_cav', 't_la', 'Pven', 'PLV', 'Part', 'PLA', 'Qla', 'Qao', 'Qper', 'Qmv', 'V_ven', 'V_cav', 'V_art', 'V_LA']
                except:
                    fdatacl = np.loadtxt(folderName+"/PV_.txt")
                    name=['tstep', 'p_cav', 'V_cav']
                return returnFunc(fdatacl[:,0], fdatacl[:,name.index(varName)])
        elif self.runMode=='unloadedGeometry':
            fdatacl = np.loadtxt(folderName+"/"+varName+".txt")
            return fdatacl
    def getCost(self,folderName,runParameters,varName,func):
        if varName[:6]=="coords":
            try:
                ref=np.loadtxt(self.matchingFolder+"/"+varName+".txt")
            except:
                ref=None
            simulated=self.getspline(folderName,varName)
            return func(simulated,ref)
        elif varName=='V_sv':
            try:
                ref=float(np.loadtxt(self.matchingFolder+"/"+varName+".txt"))
            except:
                ref=None
            valueArray=self.getspline(folderName,'V_cav',returnArray=True)
            valueArray=valueArray[valueArray[:,0]>=(valueArray[:,0].max()-runParameters['duration of one cardiac cycle in ms'])]
            cycletime=valueArray[:,0]%runParameters['duration of one cardiac cycle in ms']
            diastolictimeInd=np.argmin(np.minimum(cycletime,runParameters['duration of one cardiac cycle in ms']-cycletime))
            minval=np.min(valueArray[:,1])
            maxval=valueArray[diastolictimeInd,1]#np.max(valueArray[:,1])#
            ratio=1.#2.*(valueArray[diastolictimeInd,1]-minval)/(maxval-minval)-1
            return func(ratio*(maxval-minval),ref)
        elif varName=='pmax_LVLA':  
            try:
                ref=float(np.loadtxt(self.matchingFolder+"/"+varName+".txt"))
            except:
                ref=None
            valueArray_p_cav=self.getspline(folderName,'p_cav',returnArray=True)
            valueArray=valueArray_p_cav.copy()
            valueArray_PLA=self.getspline(folderName,'PLA',returnArray=True)
            valueArray[:,1]=valueArray_p_cav[:,1]-valueArray_PLA[:,1]
            valueArray=valueArray[valueArray_p_cav[:,0]>=(valueArray_p_cav[:,0].max()-runParameters['duration of one cardiac cycle in ms'])]
            return func(valueArray[:,1].max(),ref)
        elif varName=='pmax_LAVALVregurge':  
            try:
                ref=float(np.loadtxt(self.matchingFolder+"/"+varName+".txt"))
            except:
                ref=None
            valueArray_p_cav=self.getspline(folderName,'p_cav',returnArray=True)
            valueArray=valueArray_p_cav.copy()
            valueArray_PLA=self.getspline(folderName,'Plavalv',returnArray=True)
            valueArray[:,1]=valueArray_p_cav[:,1]-valueArray_PLA[:,1]
            valueArray=valueArray[valueArray_p_cav[:,0]>=(valueArray_p_cav[:,0].max()-runParameters['duration of one cardiac cycle in ms'])]
            return func(valueArray[:,1].max(),ref)
        elif varName=='pmax_LVAA':
            try:
                ref=float(np.loadtxt(self.matchingFolder+"/"+varName+".txt"))
            except:
                ref=None
            valueArray_p_cav=self.getspline(folderName,'p_cav',returnArray=True)
            valueArray=valueArray_p_cav.copy()
            valueArray_PAA=self.getspline(folderName,'Paa',returnArray=True)
            valueArray[:,1]=valueArray_p_cav[:,1]-valueArray_PAA[:,1]
            valueArray=valueArray[valueArray_p_cav[:,0]>=(valueArray_p_cav[:,0].max()-runParameters['duration of one cardiac cycle in ms'])]
            return func(valueArray[:,1].max(),ref)
        elif varName=='pmax_LVVALV':
            try:
                ref=float(np.loadtxt(self.matchingFolder+"/"+varName+".txt"))
            except:
                ref=None
            valueArray_p_cav=self.getspline(folderName,'p_cav',returnArray=True)
            valueArray=valueArray_p_cav.copy()
            valueArray_PAA=self.getspline(folderName,'Plvvalv',returnArray=True)
            valueArray[:,1]=valueArray_p_cav[:,1]-valueArray_PAA[:,1]
            valueArray=valueArray[valueArray_p_cav[:,0]>=(valueArray_p_cav[:,0].max()-runParameters['duration of one cardiac cycle in ms'])]
            return func(valueArray[:,1].max(),ref)
        elif varName=='pmin_LVVALV':
            try:
                ref=float(np.loadtxt(self.matchingFolder+"/"+varName+".txt"))
            except:
                ref=None
            valueArray_p_cav=self.getspline(folderName,'p_cav',returnArray=True)
            valueArray=valueArray_p_cav.copy()
            valueArray_PAA=self.getspline(folderName,'Plvvalv',returnArray=True)
            valueArray[:,1]=valueArray_p_cav[:,1]-valueArray_PAA[:,1]
            valueArray=valueArray[valueArray_p_cav[:,0]>=(valueArray_p_cav[:,0].max()-runParameters['duration of one cardiac cycle in ms'])]
            return func(valueArray[:,1].min(),ref)
        elif varName =='t_isorelax':
            try:
                ref=np.loadtxt(self.matchingFolder+"/"+varName+".txt")[0]
            except:
                ref=None
            valueArray=self.getspline(folderName,'V_cav')
            valueArray=valueArray[valueArray[:,0]>=(valueArray[:,0].max()-runParameters['duration of one cardiac cycle in ms'])]
            minVolindmax=np.argmin(valueArray[:,1])
            minVolindmin=minVolindmax
            minVol=valueArray[minVolindmax,1]
            for n in range(minVolindmax+1,valueArray.shape[0]):
                if valueArray[n,1]<=(minVol*1.01):
                    minVolindmax+=1
                else:
                    break
            for n in range(minVolindmin-1,-1,-1):
                if valueArray[n,1]<=(minVol*1.01):
                    minVolindmin-=1
                else:
                    break
            relaxtime=valueArray[minVolindmax,0]-valueArray[minVolindmin,0]
            return func(relaxtime,ref)
        elif varName=='t_isocontract':
            try:
                ref=np.loadtxt(self.matchingFolder+"/"+varName+".txt")[0]
            except:
                ref=None
            valueArray=self.getspline(folderName,'V_cav')
            valueArray=valueArray[valueArray[:,0]>=(valueArray[:,0].max()-runParameters['duration of one cardiac cycle in ms'])]
            maxVolindmax=np.argmax(valueArray[:,1])
            maxVolindmin=maxVolindmax
            maxVol=valueArray[maxVolindmax,1]
            for n in range(maxVolindmax+1,valueArray.shape[0]):
                if valueArray[n,1]<=(maxVol*0.99):
                    maxVolindmax+=1
                else:
                    break
            for n in range(maxVolindmin-1,-1,-1):
                if valueArray[n,1]<=(maxVol*0.99):
                    maxVolindmin-=1
                else:
                    break
            contracttime=valueArray[maxVolindmax,0]-valueArray[maxVolindmin,0]
            return func(contracttime,ref)
        else:
            try:
                ref=np.loadtxt(self.matchingFolder+"/"+varName+".txt")
            except:
                ref=None
            spl=self.getspline(folderName,varName)
            simulated = scipyinterpolate.splev(ref[:,0], spl)
            return func(simulated,ref[:,1])
def optimiseEDP(LVclosed_object,maxLALV_pressure,folder=None,try_stlPressure=4.):
    if folder is None:
        folder='optimiseEDP'
    os.makedirs(LVclosed_object.casename+'/'+folder+'/best_fit',exist_ok=True)
    LVclosed_object.runCount=folder
    LVclosed_object.generateMesh(LVclosed_object.defaultParameters,toRunCountFolder=True)
    tune=1.
    end_start_Diastolic_pressure=0.
    stlPressure_adj=try_stlPressure*0.3
    try_stlPressure_current=try_stlPressure*0.7
    while end_start_Diastolic_pressure*tune<maxLALV_pressure*tune:
        if abs(end_start_Diastolic_pressure-maxLALV_pressure)<0.001 or stlPressure_adj<0.001:
            break
        try_stlPressure_current+=stlPressure_adj*tune
        LVclosed_object.runCount=folder+'/{0:.3f}'.format(try_stlPressure_current)
        heart=LVclosed_object.adjustMeshToPressureWithInverse(LVclosed_object.casename+'/'+folder,LVclosed_object.meshname,try_stlPressure_current,savename=LVclosed_object.meshname+'_unloadedmesh',iterationNumber=0,saveFolder='/{0:.3f}'.format(try_stlPressure_current))
        heart=heartModelTranslator.HeartModel(LVclosed_object.heartModelStr,LVclosed_object.defaultParameters+{'path of folder with case':LVclosed_object.casename+'/'+LVclosed_object.runCount,'name of mesh':LVclosed_object.meshname+'_unloadedmesh'})
        LVclosed_object.solveVolume(heart,LVclosed_object.defaultParameters['LV end systolic volume in mL'])
        temp_ESP=heart.get_CavityPressure() 
        LVclosed_object.solveVolume(heart,LVclosed_object.defaultParameters['LV end diastolic volume in mL'])
        temp_EDP=heart.get_CavityPressure() 
        end_start_Diastolic_pressure=temp_EDP-temp_ESP
        LVclosed_object.defaultParameters['LV end diastolic pressure in mmHg']=temp_EDP
        LVclosed_object.writeRunParameters(LVclosed_object.casename+'/'+LVclosed_object.runCount,LVclosed_object.defaultParameters)
        if end_start_Diastolic_pressure*tune>maxLALV_pressure*tune:
            tune*=-1.
            stlPressure_adj*=0.3
    
    cmd = "cp -r " + LVclosed_object.casename+'/'+folder+'/'+LVclosed_object.meshname+'* '+ LVclosed_object.casename+'/'+folder+'/best_fit'
    os.system(cmd)
    cmd = "cp -r " + LVclosed_object.casename+'/'+LVclosed_object.runCount+'/* '+ LVclosed_object.casename+'/'+folder+'/best_fit'
    os.system(cmd)
def optimiseStiffness(LVclosed_object,LAtime_Vol,maxLALV_pressure,folder=None,relaxationPressureCorrection=1.,tryPressureRatio=0.5):
    if folder is None:
        folder='optimiseStiffness'
    os.makedirs(LVclosed_object.casename+'/'+folder+'/best_fit',exist_ok=True)
    
    LAtime_Vol_spl=scipyinterpolate.splrep(LAtime_Vol[:,0], LAtime_Vol[:,1])
    ageInWeeks=float(LVclosed_object.LVage[5:])
    LA_unloadedVolume=LVclosed_object.defaultParameters.getFetalPopulation_LA_unloadedVolume(ageInWeeks)
    
    count=0
    LVclosed_object.runCount=folder+'/'+str(count)
    if not(os.path.isfile(LVclosed_object.casename+'/'+LVclosed_object.runCount+'/'+LVclosed_object.meshname+'_unloadedmesh.hdf5')):
        LVclosed_object.unloadedGeometryRun(tryPressureRatio=tryPressureRatio)
    LVclosed_object.iterativeRun(endTime=LVclosed_object.defaultParameters["BCL"],setHeart=[LVclosed_object.casename+'/'+LVclosed_object.runCount,LVclosed_object.meshname+'_unloadedmesh'],outputResultList=[],trackphase=True)
    LV_tPV=np.loadtxt(LVclosed_object.casename+"/"+str(LVclosed_object.runCount)+"/PV_.txt")
    LV_tPV[:,1]*=relaxationPressureCorrection
    LA_Vol = scipyinterpolate.splev(LV_tPV[:,0], LAtime_Vol_spl)
    LA_pressures=(LA_Vol-LA_unloadedVolume)/LVclosed_object.defaultParameters['lac']
    
    LALV_pressure=(LA_pressures-LV_tPV[:,1]).max()
    with open(LVclosed_object.casename+'/'+LVclosed_object.runCount+"/runParameters.txt", "r") as f:
         old = f.readlines() # read everything in the file
    old.insert(0,"#current maxLALV pressure= "+repr(LALV_pressure)+"/"+repr(relaxationPressureCorrection)+" , target LALV pressure "+repr(maxLALV_pressure)+" \n")
    with open(LVclosed_object.casename+'/'+LVclosed_object.runCount+"/runParameters.txt", "w") as f:
         f.writelines(old)
    tune=1.
    if maxLALV_pressure>LALV_pressure:
        relation=1
    else:
        relation=-1
    try_strainCoef=LVclosed_object.defaultParameters["StrainEnergyDensityFunction_Coef"]
    adj_strainCoef=try_strainCoef*0.1
    while maxLALV_pressure*tune*relation>LALV_pressure*tune*relation:
        if abs(maxLALV_pressure-LALV_pressure)<0.001:
            break
        try_strainCoef+=tune*adj_strainCoef
        count+=1
        LVclosed_object.defaultParameters["StrainEnergyDensityFunction_Coef"]=try_strainCoef
        LVclosed_object.runCount=folder+'/'+str(count)
        LVclosed_object.unloadedGeometryRun(tryPressureRatio=tryPressureRatio)
        LVclosed_object.iterativeRun(endTime=LVclosed_object.defaultParameters["BCL"],setHeart=[LVclosed_object.casename+'/'+LVclosed_object.runCount,LVclosed_object.meshname+'_unloadedmesh'],outputResultList=[],trackphase=True)
        LV_tPV=np.loadtxt(LVclosed_object.casename+"/"+str(LVclosed_object.runCount)+"/PV_.txt")
        LV_tPV[:,1]*=relaxationPressureCorrection
        LALV_pressure=(LA_pressures-LV_tPV[:,1]).max()
        with open(LVclosed_object.casename+'/'+LVclosed_object.runCount+"/runParameters.txt", "r") as f:
             old = f.readlines() # read everything in the file
        old.insert(0,"#current maxLALV pressure= "+repr(LALV_pressure)+"/"+repr(relaxationPressureCorrection)+" , target LALV pressure "+repr(maxLALV_pressure)+" \n")
        with open(LVclosed_object.casename+'/'+LVclosed_object.runCount+"/runParameters.txt", "w") as f:
             f.writelines(old)
        if maxLALV_pressure*tune*relation<LALV_pressure*tune*relation:
            tune*=-1.
            adj_strainCoef*=0.3
    
    cmd = "cp -r " + LVclosed_object.casename+'/'+LVclosed_object.runCount+'/* '+ LVclosed_object.casename+'/'+folder+'/best_fit'
    os.system(cmd)
def optimiseFiberAngle(LVclosed_object,optVarList=None,folder=None,tryPressureRatio=0.2):
    if optVarList is None:
        optVarList=['endo_angle','epi_angle']
    if folder is None:
        folder='optimiseFiberAngle'
    bound=[]
    for var in optVarList:
        if var=='endo_angle':
            bound.append((0.,90.))
        elif var=='epi_angle':
            bound.append((-90.,0))
        elif var=='fiberSheetletAngle':
            bound.append((-90.,90.))
        elif var=='fiberSheetletWidth':
            bound.append((0.,float('inf')))
        elif var=='radialFiberAngle':
            bound.append((-90.,90.))
        else:
            bound.append((float('-inf'),float('inf')))
    LVclosed_object.defaultRunMode='iterativeRun'
    LVclosed_object.runCount=folder+'/0'
    os.makedirs(LVclosed_object.casename+'/'+folder+'/ref',exist_ok=True)
    np.savetxt(LVclosed_object.casename+'/'+folder+'/ref/p_cav.txt',np.array([[0.,LVclosed_object.defaultParameters['LV end diastolic pressure in mmHg']],[LVclosed_object.defaultParameters['time to maximum LV fiber tension in ms'],500.]]))
    cost=cost_function(LVclosed_object.defaultRunMode,matchingFolder=LVclosed_object.casename+'/'+folder+'/ref')
    cost.addCost('LeftVentriclePressureMeanSquare')
    linker=optimiser_linker(LVclosed_object,cost)
    linker.addVariable(optVarList)
    linker.set_kwargs(endTime=LVclosed_object.defaultParameters['time to maximum LV fiber tension in ms']+3,outputResultList=[],tryPressureRatio=tryPressureRatio)
    initPara=[]
    for n in range(len(optVarList)):
        initPara.append(LVclosed_object.defaultParameters[optVarList[n]])
    iniPara=np.array(initPara)
    #optimize.fmin(linker,iniPara,xtol=0.001, ftol=0.001)
    optimize.minimize(linker,iniPara,bounds=bound,tol=0.001)

def optimiseTact(LVclosed_object,time_in_ms,pressure_in_mmHg,folder=None,setHeart=None):
    if folder is None:
        folder='optimiseTact'
    LVclosed_object.defaultRunMode='iterativeRun'
    LVclosed_object.runCount=folder+'/0'
    os.makedirs(LVclosed_object.casename+'/'+folder+'/ref',exist_ok=True)
    np.savetxt(LVclosed_object.casename+'/'+folder+'/ref/p_cav.txt',np.array([[0.,LVclosed_object.defaultParameters['LV end diastolic pressure in mmHg']],[time_in_ms,pressure_in_mmHg]]))
    cost=cost_function(LVclosed_object.defaultRunMode,matchingFolder=LVclosed_object.casename+'/'+folder+'/ref')
    cost.addCost('LeftVentriclePressureMeanSquare')
    linker=optimiser_linker(LVclosed_object,cost)
    linker.addVariable(['maximum LV fiber tension in Pa'])
    linker.set_kwargs(endTime=time_in_ms+3,outputResultList=[],setHeart=setHeart)
    iniPara=np.array([LVclosed_object.defaultParameters['maximum LV fiber tension in Pa']])
    optimize.minimize(linker,iniPara,bounds=[(1e4,LVclosed_object.defaultParameters['maximum LV fiber tension in Pa']*1.1)],tol=0.001)

def optimiseWinkesselParameters(LVclosed_object,strokeVolume,maxLVLApressureDiff,maxLVAApressureDiff,minLVAApressureDiff=None,step=0,stopstep=3,Q_regurge_flowratio=0.1,mitralinflow_pressure=None,optVarList=None,costList=None,baseFolder=None,folder=None,folderToLVbehavior=None,runCycles=20):
    if folder is None:
        folder='optimiseWinkesselParameters'
    if baseFolder is not None:
        folder=baseFolder+'/'+folder
    if optVarList is None:
        optVarList=['lvlvvalvk','scale Windkessel active pressure']
        logs=[3.,0]
    else:
        optVarList=['lvlvvalvk','scale Windkessel active pressure']+optVarList
        logs=[3.,0]+list([0]*len(optVarList))
    if costList is None:
        costList=['StrokeVolumeMeanSquare','MaxLVVALVPressureMeanSquare']
    if maxLVLApressureDiff is not None:
        optVarList=optVarList+['lvregurgevalveratio']
        logs=logs+[3.]
        costList=costList+['MaxLAVALVPressureMeanSquare']
    if minLVAApressureDiff is not None:
        optVarList=optVarList+['aaregurgevalveratio']
        logs=logs+[3.]
        costList=costList+['MinLVVALVPressureMeanSquare']
    os.makedirs(LVclosed_object.casename+'/'+folder+'/ref',exist_ok=True)
    weights=[]
    for n in range(len(costList)):
        if 'StrokeVolume' in costList[n]:
            weights.append(10./LVclosed_object.defaultParameters['LV end systolic volume in mL'])
            np.savetxt(LVclosed_object.casename+'/'+folder+'/ref/V_sv.txt',np.array([strokeVolume]))
        if 'MaxLAVALVPressure' in costList[n] and maxLVLApressureDiff is not None:
            weights.append(1./maxLVLApressureDiff)
            np.savetxt(LVclosed_object.casename+'/'+folder+'/ref/pmax_LAVALVregurge.txt',np.array([maxLVLApressureDiff]))
        if 'MaxLVVALVPressure' in costList[n]:
            if maxLVLApressureDiff is not None:
                weights.append(2./(maxLVLApressureDiff+maxLVAApressureDiff))
            else:
                weights.append(2./(maxLVAApressureDiff))
            np.savetxt(LVclosed_object.casename+'/'+folder+'/ref/pmax_LVVALV.txt',np.array([maxLVAApressureDiff]))
        if 'MinLVVALVPressure' in costList[n]:
            np.savetxt(LVclosed_object.casename+'/'+folder+'/ref/pmin_LVVALV.txt',np.array([minLVAApressureDiff]))
    LVclosed_object.defaultRunMode='flowrateWindkesselRun'
    LVLA_Q=LVclosed_object.manualVol(np.arange(LVclosed_object.defaultParameters['duration of one cardiac cycle in ms']*10.)/10.)
    LVLA_Q=(LVLA_Q[1:]-LVLA_Q[:-1])/0.1
    maxLVLA_Q=-LVLA_Q.min()
    if mitralinflow_pressure is not None:
        guess_lalavalvk=(mitralinflow_pressure/(LVLA_Q.max()**2.))
        LVclosed_object.defaultParameters['lalavalvk']=guess_lalavalvk
    ref_value=[maxLVLApressureDiff,maxLVAApressureDiff]
    if step<=0 and stopstep>=0:
        os.makedirs(LVclosed_object.casename+'/'+folder+'/init_Qratio/best_fit',exist_ok=True)
        value_count=0
        
        if Q_regurge_flowratio<0:
            forwardguess=-1
        else:
            forwardguess=1
        adj_Q_regurge_flowratio=forwardguess*Q_regurge_flowratio*0.3
        try_Q_regurge_flowratio=forwardguess*Q_regurge_flowratio-forwardguess*adj_Q_regurge_flowratio
        real_Qratio=1.0
        tune=1
        while try_Q_regurge_flowratio*tune<real_Qratio*tune:
            if abs(real_Qratio-try_Q_regurge_flowratio)<0.0001:
                break
            value_count+=1
            if (try_Q_regurge_flowratio+tune*forwardguess*adj_Q_regurge_flowratio)<=0:
                adj_Q_regurge_flowratio=try_Q_regurge_flowratio*0.3
            elif (try_Q_regurge_flowratio+tune*forwardguess*adj_Q_regurge_flowratio)>=1:
                adj_Q_regurge_flowratio=(1.-try_Q_regurge_flowratio)*0.3
            try_Q_regurge_flowratio+=tune*forwardguess*adj_Q_regurge_flowratio
            LVclosed_object.runCount=folder+'/init_Qratio/'+str(value_count)
            
            if maxLVLApressureDiff is not None:
                if (maxLVAApressureDiff-0.001*(maxLVLA_Q*(1.-try_Q_regurge_flowratio)))<0:
                    guess_lvlvvalvk=((maxLVAApressureDiff)/((maxLVLA_Q*(1.-try_Q_regurge_flowratio))**2.))
                else:
                    guess_lvlvvalvk=((maxLVAApressureDiff-0.001*(maxLVLA_Q*(1.-try_Q_regurge_flowratio)))/((maxLVLA_Q*(1.-try_Q_regurge_flowratio))**2.))
            
                if (maxLVLApressureDiff-0.001*(maxLVLA_Q*(1.-try_Q_regurge_flowratio)))<0:
                    guess_lvregurgevalveratio=LVclosed_object.defaultParameters['lalavalvk']/((maxLVLApressureDiff)/((maxLVLA_Q*try_Q_regurge_flowratio)**2.))
                else:
                    guess_lvregurgevalveratio=LVclosed_object.defaultParameters['lalavalvk']/((maxLVLApressureDiff-0.001*(maxLVLA_Q*(1.-try_Q_regurge_flowratio)))/((maxLVLA_Q*try_Q_regurge_flowratio)**2.))
                LVclosed_object.defaultParameters['lvlvvalvk']=guess_lvlvvalvk
                LVclosed_object.defaultParameters['lvregurgevalveratio']=guess_lvregurgevalveratio
            else:
                LVclosed_object.defaultParameters['lvlvvalvk']=((maxLVAApressureDiff)/(maxLVLA_Q**2.))
                LVclosed_object.defaultParameters['lvregurgevalveratio']=-1
            LVclosed_object.runCount=folder+'/init_Qratio/'+str(value_count)
            LVclosed_object.flowrateWindkesselRun()
            QLV_maxind=np.argmin(getCircuitResult(LVclosed_object.casename+'/'+LVclosed_object.runCount,'i(Vlvprobe)',last_cycle=True))
            QLVLA=getCircuitResult(LVclosed_object.casename+'/'+LVclosed_object.runCount,'i(Vlvregurgitation)',last_cycle=True)[QLV_maxind]
            QLVAA=getCircuitResult(LVclosed_object.casename+'/'+LVclosed_object.runCount,'i(Vlvlvvalv)',last_cycle=True)[QLV_maxind]
            VLV=getCircuitResult(LVclosed_object.casename+'/'+LVclosed_object.runCount,'v(lv)',last_cycle=True)
            VLA=getCircuitResult(LVclosed_object.casename+'/'+LVclosed_object.runCount,'v(lavalvp1)',last_cycle=True)
            VLVVALV=getCircuitResult(LVclosed_object.casename+'/'+LVclosed_object.runCount,'v(lvvalv)',last_cycle=True)
            real_P_value=[(VLV-VLA).max(),(VLV-VLVVALV).max()]
            
            real_Qratio=QLVLA/(QLVLA+QLVAA)
            with open(LVclosed_object.casename+'/'+LVclosed_object.runCount+"/runParameters.txt", "r") as f:
                 old = f.readlines() # read everything in the file
            old.insert(0,"#Real total Q= "+repr(QLVLA+QLVAA)+" , guess total Q "+repr(maxLVLA_Q)+" \n")
            old.insert(0,"#PLVLA,PLVAA= "+repr(real_P_value)+" , to match "+repr(ref_value)+" , try Qratio="+str(try_Q_regurge_flowratio)+' , real_Qratio= '+str(real_Qratio)+" \n")
            with open(LVclosed_object.casename+'/'+LVclosed_object.runCount+"/runParameters.txt", "w") as f:
                 f.writelines(old)
            if try_Q_regurge_flowratio*tune>real_Qratio*tune:
                tune*=-1
                adj_Q_regurge_flowratio*=0.3
            #maxLVLA_Q=0.5*(maxLVLA_Q+QLVLA+QLVAA)
            if maxLVLApressureDiff is None:
                break
        #copy over
        cmd = "cp -r " + LVclosed_object.casename+'/'+LVclosed_object.runCount+'/* '+ LVclosed_object.casename+'/'+folder+'/init_Qratio/best_fit'
        os.system(cmd)
    #initialise lvlvvalvk' with  'MaxLVAAPressure'
    if step<=1 and stopstep>=1:
        #get lvlvvalvk' circuit results with volume data
        if step==1 and os.path.isfile(LVclosed_object.casename+'/'+folder+'/init_lvlvvalvk/best_fit/runParameters.txt'):
            results=LVclosed_object.readRunParameters(LVclosed_object.casename+'/'+folder+'/init_lvlvvalvk/best_fit')
        else:
            results=LVclosed_object.readRunParameters(LVclosed_object.casename+'/'+folder+'/init_Qratio/best_fit')
        LVclosed_object.defaultParameters['lvlvvalvk']=results['lvlvvalvk']
        if 'aaregurgevalveratio' in optVarList and LVclosed_object.defaultParameters['aaregurgevalveratio']<0:
            LVclosed_object.defaultParameters['aaregurgevalveratio']=10.**-5.
        LVclosed_object.defaultParameters['lvregurgevalveratio']=results['lvregurgevalveratio']
        if 'lvlvvalvr' in optVarList:
            LVclosed_object.defaultParameters['lvlvvalvr']=results['lvlvvalvr']
        os.makedirs(LVclosed_object.casename+'/'+folder+'/init_lvlvvalvk/ref',exist_ok=True)
        for n in range(len(costList)):
            if 'MaxLAVALVPressure' in costList[n] and maxLVLApressureDiff is not None:
                np.savetxt(LVclosed_object.casename+'/'+folder+'/init_lvlvvalvk/ref/pmax_LAVALVregurge.txt',np.array([maxLVLApressureDiff]))
            if 'MaxLVVALVPressure' in costList[n]:
                np.savetxt(LVclosed_object.casename+'/'+folder+'/init_lvlvvalvk/ref/pmax_LVVALV.txt',np.array([maxLVAApressureDiff]))
            if 'MinLVVALVPressure' in costList[n]:
                np.savetxt(LVclosed_object.casename+'/'+folder+'/init_lvlvvalvk/ref/pmin_LVVALV.txt',np.array([minLVAApressureDiff]))
        
        for n in [1]:
            LVclosed_object.runCount=folder+'/init_lvlvvalvk/0'
            cost=cost_function(LVclosed_object.defaultRunMode,matchingFolder=LVclosed_object.casename+'/'+folder+'/init_lvlvvalvk/ref')
            if maxLVLApressureDiff is None:
                cost.addCost(['MaxLVVALVPressureMeanSquare'],weights=[2./(maxLVAApressureDiff)])
            else:
                cost.addCost(['MaxLVVALVPressureMeanSquare'],weights=[2./(maxLVLApressureDiff+maxLVAApressureDiff)])
            if 'lvregurgevalveratio' in optVarList:
                cost.addCost(['MaxLAVALVPressureMeanSquare'],weights=[1./maxLVLApressureDiff])
            if 'aaregurgevalveratio' in optVarList:
                cost.addCost(['MinLVVALVPressureMeanSquare'],weights=[-1./minLVAApressureDiff])
            linker=optimiser_linker(LVclosed_object,cost)
            linker.addVariable(['lvlvvalvk'],[(LVclosed_object.defaultParameters['lvlvvalvk']/3.**n,LVclosed_object.defaultParameters['lvlvvalvk']*3.**n)],log=[3.])
            linker_para=[LVclosed_object.defaultParameters['lvlvvalvk']]
            if 'aaregurgevalveratio' in optVarList:
                linker.addVariable(['aaregurgevalveratio'],[(LVclosed_object.defaultParameters['aaregurgevalveratio']/10.**n,LVclosed_object.defaultParameters['aaregurgevalveratio']*10.**n)],log=[10.])
                linker_para.append(LVclosed_object.defaultParameters['aaregurgevalveratio'])
            if 'lvregurgevalveratio' in optVarList:
                linker.addVariable(['lvregurgevalveratio'],[(LVclosed_object.defaultParameters['lvregurgevalveratio']/10.**n,LVclosed_object.defaultParameters['lvregurgevalveratio']*10.**n)],log=[10.])
                linker_para.appends(LVclosed_object.defaultParameters['lvregurgevalveratio'])
            if 'lvlvvalvr' in optVarList:
                linker.addVariable(['lvlvvalvr'],[(LVclosed_object.defaultParameters['lvlvvalvr']/2.**n,LVclosed_object.defaultParameters['lvlvvalvr']*3.**n)],log=[3.])
                linker_para.appends(LVclosed_object.defaultParameters['lvlvvalvr'])
            optimize.minimize(linker,linker.invOperate_para(np.array(linker_para)),jac='3-point',options={'eps':0.00001,'maxiter':8},tol=0.0001)
            results=LVclosed_object.readRunParameters(LVclosed_object.casename+'/'+folder+'/init_lvlvvalvk/best_fit')
            LVclosed_object.defaultParameters['lvlvvalvk']=results['lvlvvalvk']
            if 'aaregurgevalveratio' in optVarList:
                LVclosed_object.defaultParameters['aaregurgevalveratio']=results['aaregurgevalveratio']
            if 'lvregurgevalveratio' in optVarList:
                LVclosed_object.defaultParameters['lvregurgevalveratio']=results['lvregurgevalveratio']
            if 'lvlvvalvr' in optVarList:
                LVclosed_object.defaultParameters['lvlvvalvr']=results['lvlvvalvr']
    if os.path.isfile(LVclosed_object.casename+'/'+folder+'/init_lvlvvalvk/best_fit/runParameters.txt'):
        results=LVclosed_object.readRunParameters(LVclosed_object.casename+'/'+folder+'/init_lvlvvalvk/best_fit')
    else:
        results=LVclosed_object.readRunParameters(LVclosed_object.casename+'/'+folder+'/init_Qratio/best_fit')
    LVclosed_object.defaultParameters['lvlvvalvk']=results['lvlvvalvk']
    if 'aaregurgevalveratio' in optVarList:
        LVclosed_object.defaultParameters['aaregurgevalveratio']=results['aaregurgevalveratio']
    if 'lvregurgevalveratio' in optVarList:
        LVclosed_object.defaultParameters['lvregurgevalveratio']=results['lvregurgevalveratio']
    if 'lvlvvalvr' in optVarList:
        LVclosed_object.defaultParameters['lvlvvalvr']=results['lvlvvalvr']
    LVclosed_object.defaultRunMode='fullWindkesselRun'
    if step<=2 and stopstep>=2:
    
        #get appropriate Windkessel_scale_maximum_fiber_tension
        if step==2 and os.path.isfile(LVclosed_object.casename+'/'+folder+'/init_WindScale/best_fit/runParameters.txt'):
            init_Windkessel_scale_maximum_fiber_tension=LVclosed_object.readRunParameters(LVclosed_object.casename+'/'+folder+'/init_WindScale/best_fit')['scale Windkessel active pressure']
        else:
            LVclosed_object.runCount=folder+'/init_WindScale/0'
            LVclosed_object.fullWindkesselRun(unloadGeo=True,folderToLVbehavior=folderToLVbehavior,runCycles=runCycles)
            results_LVBehavior=getCircuitResult(LVclosed_object.casename+'/'+LVclosed_object.runCount,'v(lvvol)',last_cycle=True)
            temp_SV=results_LVBehavior.max()-results_LVBehavior.min()
            init_Windkessel_scale_maximum_fiber_tension=min(0.9999999,max(0.0000001,LVclosed_object.defaultParameters['scale Windkessel active pressure']*strokeVolume/temp_SV))
        #initialise 'scale Windkessel active pressure' with Stroke Volume
        os.makedirs(LVclosed_object.casename+'/'+folder+'/init_WindScale/ref',exist_ok=True)
        np.savetxt(LVclosed_object.casename+'/'+folder+'/init_WindScale/ref/V_sv.txt',np.array([strokeVolume]))
        LVclosed_object.runCount=folder+'/init_WindScale/0'
        cost=cost_function(LVclosed_object.defaultRunMode,matchingFolder=LVclosed_object.casename+'/'+folder+'/init_WindScale/ref')
        cost.addCost(['StrokeVolumeMeanSquare'])
        linker=optimiser_linker(LVclosed_object,cost)
        linker.addVariable(['scale Windkessel active pressure'],[(init_Windkessel_scale_maximum_fiber_tension*0.5,min(1,init_Windkessel_scale_maximum_fiber_tension*1.5))])
        linker.set_kwargs(folderToLVbehavior=folderToLVbehavior,runCycles=runCycles)
        bounds=[(float('-inf'),float('inf'))]
        optimize.minimize(linker,linker.invOperate_para(np.array([init_Windkessel_scale_maximum_fiber_tension])),jac='3-point',options={'eps':0.000001,'maxiter':1},tol=0.0001)
        results=LVclosed_object.readRunParameters(LVclosed_object.casename+'/'+folder+'/init_WindScale/best_fit')['scale Windkessel active pressure']
        LVclosed_object.defaultParameters['scale Windkessel active pressure']=results
        for n in range(1):
            if abs(results-strokeVolume)/strokeVolume<0.01:
                break
            LVclosed_object.runCount=folder+'/init_WindScale/0'
            cost=cost_function(LVclosed_object.defaultRunMode,matchingFolder=LVclosed_object.casename+'/'+folder+'/init_WindScale/ref')
            cost.addCost(['StrokeVolumeMeanSquare'])
            linker=optimiser_linker(LVclosed_object,cost)
            linker.addVariable(['scale Windkessel active pressure'],[(init_Windkessel_scale_maximum_fiber_tension*0.5,min(1,init_Windkessel_scale_maximum_fiber_tension*1.5))])
            linker.set_kwargs(folderToLVbehavior=folderToLVbehavior,runCycles=runCycles)
            bounds=[(float('-inf'),float('inf'))]
            optimize.minimize(linker,linker.invOperate_para(np.array([init_Windkessel_scale_maximum_fiber_tension])),jac='3-point',options={'eps':0.000001,'maxiter':2},tol=0.0001)
            results=LVclosed_object.readRunParameters(LVclosed_object.casename+'/'+folder+'/init_WindScale/best_fit')['scale Windkessel active pressure']
            LVclosed_object.defaultParameters['scale Windkessel active pressure']=results
    results=LVclosed_object.readRunParameters(LVclosed_object.casename+'/'+folder+'/init_WindScale/best_fit')['scale Windkessel active pressure']
        
    LVclosed_object.defaultParameters['scale Windkessel active pressure']=results
    if step<=3 and stopstep>=3:
        #optimise all
        if step==3 and os.path.isfile(LVclosed_object.casename+'/'+folder+'/best_fit/runParameters.txt'):
            results=LVclosed_object.readRunParameters(LVclosed_object.casename+'/'+folder+'/best_fit')
            for n in range(len(optVarList)):
                LVclosed_object.defaultParameters[optVarList[n]]=results[optVarList[n]]
        for m in range(3):
            LVclosed_object.runCount=folder+'/0'
            iniPara=[]
            bounds=[]
            all_limits=[]
            for n in range(0,len(optVarList)):
                if optVarList[n]=='lvlvvalvk':
                    all_limits.append((LVclosed_object.defaultParameters['lvlvvalvk']/3.,LVclosed_object.defaultParameters['lvlvvalvk']*3.))
                    iniPara.append(LVclosed_object.defaultParameters['lvlvvalvk'])
                    bounds.append((float('-inf'),float('inf')))
                if optVarList[n]=='lvregurgevalveratio':
                    all_limits.append((LVclosed_object.defaultParameters['lvregurgevalveratio']/3.,min(1.,LVclosed_object.defaultParameters['lvregurgevalveratio']*3.)))
                    iniPara.append(LVclosed_object.defaultParameters['lvregurgevalveratio'])
                    bounds.append((float('-inf'),float('inf')))
                if optVarList[n]=='scale Windkessel active pressure':
                    all_limits.append((LVclosed_object.defaultParameters['scale Windkessel active pressure']*0.5,min(1.,LVclosed_object.defaultParameters['scale Windkessel active pressure']*1.5)))
                    iniPara.append(LVclosed_object.defaultParameters['scale Windkessel active pressure'])
                    bounds.append((float('-inf'),float('inf')))
            iniPara=np.array(iniPara)
            
            cost=cost_function(LVclosed_object.defaultRunMode,matchingFolder=LVclosed_object.casename+'/'+folder+'/ref')
            cost.addCost(costList,weights=weights)
            linker=optimiser_linker(LVclosed_object,cost)
            linker.addVariable(optVarList,all_limits,log=logs)
            linker.set_kwargs(folderToLVbehavior=folderToLVbehavior,runCycles=runCycles)
            optimize.minimize(linker,linker.invOperate_para(iniPara),jac='3-point',options={'eps':0.000001,'maxiter':10},tol=0.001)
            
            results=LVclosed_object.readRunParameters(LVclosed_object.casename+'/'+folder+'/best_fit')
            for n in range(len(optVarList)):
                LVclosed_object.defaultParameters[optVarList[n]]=results[optVarList[n]]
            

def optimiseAllWinkesselParameters(LVclosed_object,strokeVolume,maxLVLApressureDiff,maxLVLVVALVpressureDiff,Q_regurge_flowratio=0.1,iteration=2,folderToLVbehavior=None,runCycles=20):
    baseFolder='optimiseWinkesselParameters'
    folder='0'
    EDVol=LVclosed_object.defaultParameters['LV end diastolic volume in mL']
    ESVol=LVclosed_object.defaultParameters['LV end systolic volume in mL']
    os.makedirs(LVclosed_object.casename+'/'+baseFolder+'/'+folder,exist_ok=True)
    if not(os.path.isfile(LVclosed_object.casename+'/'+baseFolder+'/'+folder+'/init_WindScale/best_fit/circuit_results_lastcycle.txt')):
        if not(os.path.isfile(LVclosed_object.casename+'/'+baseFolder+'/'+folder+'/init_lvlvvalvk/best_fit/circuit_results_lastcycle.txt')):
            optimiseWinkesselParameters(LVclosed_object,strokeVolume,maxLVLApressureDiff,maxLVLVVALVpressureDiff,step=0,stopstep=2,Q_regurge_flowratio=Q_regurge_flowratio,baseFolder=baseFolder,folder=folder,folderToLVbehavior=folderToLVbehavior,runCycles=runCycles)
        else:
            optimiseWinkesselParameters(LVclosed_object,strokeVolume,maxLVLApressureDiff,maxLVLVVALVpressureDiff,step=2,stopstep=2,Q_regurge_flowratio=Q_regurge_flowratio,baseFolder=baseFolder,folder=folder,folderToLVbehavior=folderToLVbehavior,runCycles=runCycles)
    for startiter in range(1,10):
        if not(os.path.isfile(LVclosed_object.casename+'/'+baseFolder+'/'+str(startiter)+'/init_WindScale/best_fit/circuit_results_lastcycle.txt')):
            break
    for n in range(startiter,startiter+iteration):
        if n == (startiter+iteration-1):
            stopstep=3
        else:
            stopstep=2
        folder=str(n)
        if not(os.path.isfile(LVclosed_object.casename+'/'+baseFolder+'/'+folder+'/init_lvlvvalvk/best_fit/circuit_results_lastcycle.txt')):
            lastResultsPara=LVclosed_object.readRunParameters(LVclosed_object.casename+'/'+baseFolder+'/'+str(n-1)+'/init_WindScale/best_fit')
            for key in lastResultsPara:
                LVclosed_object.defaultParameters[key] = lastResultsPara[key]
            QLV_maxind=np.argmin(getCircuitResult(LVclosed_object.casename+'/'+baseFolder+'/'+str(n-1)+'/init_WindScale/best_fit','i(Vlvprobe)',last_cycle=True))
            QLVLA=getCircuitResult(LVclosed_object.casename+'/'+baseFolder+'/'+str(n-1)+'/init_WindScale/best_fit','i(Vlvregurgitation)',last_cycle=True)[QLV_maxind]
            QLVAA=getCircuitResult(LVclosed_object.casename+'/'+baseFolder+'/'+str(n-1)+'/init_WindScale/best_fit','i(Vlvlvvalv)',last_cycle=True)[QLV_maxind]
            real_Qratio=QLVLA/(QLVLA+QLVAA)
            
            os.makedirs(LVclosed_object.casename+'/'+baseFolder+'/'+folder,exist_ok=True)
            time=np.loadtxt(LVclosed_object.casename+'/'+baseFolder+'/'+str(n-1)+'/init_WindScale/best_fit/circuit_results_lastcycle.txt')[:,0]
            Vol_LV=getCircuitResult(LVclosed_object.casename+'/'+baseFolder+'/'+str(n-1)+'/init_WindScale/best_fit','v(lvvol)',last_cycle=True)
            while time[0]<0:
                time=time[1:]
                Vol_LV=Vol_LV[1:]
            while time[-1]>LVclosed_object.defaultParameters['duration of one cardiac cycle in ms']:
                time=time[:-1]
                Vol_LV=Vol_LV[:-1]
            if time[0]>0:
                time=np.concatenate(([0.],time))
                Vol_LV=np.concatenate(([Vol_LV[0]],Vol_LV))
            argmin=np.argmin(Vol_LV)
            Vol_LV[:argmin]=(Vol_LV[:argmin]-Vol_LV[argmin])/(Vol_LV[0]-Vol_LV[argmin])*(EDVol-ESVol)+ESVol
            Vol_LV[argmin:]=(Vol_LV[argmin:]-Vol_LV[argmin])/(Vol_LV[-1]-Vol_LV[argmin])*(EDVol-ESVol)+ESVol
            if time[-1]<LVclosed_object.defaultParameters['duration of one cardiac cycle in ms']:
                time=np.concatenate((time,[LVclosed_object.defaultParameters['duration of one cardiac cycle in ms']]))
                Vol_LV=np.concatenate((Vol_LV,[EDVol]))
            #Vol_LV=(Vol_LV-Vol_LV.min())/(Vol_LV[0]-Vol_LV.min())*(EDVol-ESVol)+ESVol
            np.savetxt(LVclosed_object.casename+'/'+baseFolder+'/'+str(n)+'/volume.txt',np.array([time,Vol_LV]).T)
            Q_regurge_flowratio=real_Qratio*0.7+Q_regurge_flowratio*0.3
            LVclosed_object.setManualVolumeFile(LVclosed_object.casename+'/'+baseFolder+'/'+str(n)+'/volume.txt')
            optimiseWinkesselParameters(LVclosed_object,strokeVolume,maxLVLApressureDiff,maxLVLVVALVpressureDiff,step=0,stopstep=stopstep,Q_regurge_flowratio=Q_regurge_flowratio,optVarList=None,baseFolder=baseFolder,folder=folder,folderToLVbehavior=folderToLVbehavior,runCycles=runCycles)
        else:
            lastResultsPara=LVclosed_object.readRunParameters(LVclosed_object.casename+'/'+baseFolder+'/'+folder+'/init_lvlvvalvk/best_fit')
            for key in lastResultsPara:
                LVclosed_object.defaultParameters[key] = lastResultsPara[key]
            LVclosed_object.setManualVolumeFile(LVclosed_object.casename+'/'+baseFolder+'/'+str(n)+'/volume.txt')
            optimiseWinkesselParameters(LVclosed_object,strokeVolume,maxLVLApressureDiff,maxLVLVVALVpressureDiff,step=2,stopstep=stopstep,Q_regurge_flowratio=Q_regurge_flowratio,optVarList=None,baseFolder=baseFolder,folder=folder,folderToLVbehavior=folderToLVbehavior,runCycles=runCycles)
    