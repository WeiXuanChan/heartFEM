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
                                                            - added 'Windkessel_scale_T0_LV' to sclae LVbehavior during table generation for ngspice
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
'''
_version='3.5.1'
import logging
logger = logging.getLogger('heartFEM v'+_version)
logger.info('heartFEM version '+_version)

import sys

import os as os
#import shutil
try:
    import dolfin as fenics
except:
    import fenics as fenics
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
from heartFEM import ngspice_py
from heartFEM import heartParameters
from mpi4py import MPI as pyMPI
try:
    from heartFEM import lcleeHeart
except Exception as e:
    logger.warning('Unable to load Module "lcleeHeart".')
    logger.warning(repr(e))
import motionSegmentation.BsplineFourier as bsf
 
WindkesselComponents=['lv','la','rv','ra','aa','ao1','ao2','ao3','ao4','br','ca','ub','he','inte','ivc','kid','leg','lung','pa1','pa2','plac','svc','uv']
WindkessellinkComponents=['aaao1','ao1ao2','ao2ao3','ao3ao4','pa1pa2','pa2lung','da','ao1ca','cabr','brsvc','ao1ub','ubsvc','ao3he','ao3inte','intehe','ao3kid','kidivc','ao4plac','placuv','ao4leg','legivc','uvhe','heivc','dv','svcra','ivcra','lungla','fo','raravalv','rvrvvalv','lvlvvalv','rvrvvalv']
defaultAgeScalePower={'defaultr':-1.,'pa2lungr':-1.2,'lunglar':-1.2,'cabrr':-1.1,'brsvcr':-1.1,'dvr':-0.55,
                      'defaultl':-0.33,
                      'defaultc':1.33,'brc':1.471,'lungc':1.6,'ra':0.5,'la':0.5,
                      'defaultk':0.,'fok':-0.6,'dak': -2.5,'dvk':-0.88,'raravalvk':-1.33,'rvrvvalvk':-1.33,'lalavalvk':-1.33,'lvlvvalvk':-1.33,
                      'defaultb':0.}

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
class LVclosed:
    def __new__(cls, *args, **kwargs):
        return super().__new__(cls)
    def __init__(self,defaultParameters=None,defaultAge='fetal28'):
        self.defaultRunMode='iterativeRun'
        self.defaultRun_kwargs=None
        #define defaults
        self.meshname = "t0"
        self.casename = os.getcwd()
        self.numofCycle=12
        self.LVage=defaultAge
        self.heartMotionFile=None
        
        self.defaultParameters=heartParameters.heartParameters(defaultParameters=defaultParameters,defaultAge=defaultAge)
        self.manualVol=None
        self.manualPhaseTime=None
        
        self.constrains=['Kspring','translation','rotation']
        self.runCount=0
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
        try:
            data=np.loadtxt(self.manualVolFile)
        except:
            data=np.loadtxt(self.manualVolFile,skiprows=1)
        timedata=data[:,self.manualVolFiletimecol]*self.manualVolFiletimescale
        voldata=data[:,self.manualVolFilevolcol]*self.manualVolFilevolscale
        if isinstance(time,(int,float)):
            singleValued=True
            time=np.array([time])
        else:
            singleValued=False
            time=time.copy()
        for n in range(len(time)):
            while time[n]>timedata[-1]:
                time[n]-=self.defaultParameters['BCL']
            while time[n]<timedata[0]:
                time[n]+=self.defaultParameters['BCL']
        spl=scipyinterpolate.splrep(timedata, voldata)
        if singleValued:
            return scipyinterpolate.splev(time, spl)[0]
        else:
            return scipyinterpolate.splev(time, spl)
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
            time-=self.defaultParameters['BCL']
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
    def setAtrialpulse(self,amp,peaktime,width,atrial_side='rl'):
        #atrial_side = "r" or "l" or 'lr" or "rl"
        if 'l' in atrial_side:
            self.defaultParameters['laamp']=amp
            self.defaultParameters['lapeaktime']=peaktime
            self.defaultParameters['lawidth']=width
        if 'r' in atrial_side:
            self.defaultParameters['raamp']=amp
            self.defaultParameters['rapeaktime']=peaktime
            self.defaultParameters['rawidth']=width
    def setRVcurrentInflowFourier(self,filename,fourierTerms=1,fineAdjustPhase=False):
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
            logger.info("Adjusted RVcurrentInflowFourier by "+str(adjust_x)+" radian")
        else:
            adjust_x=0.
        self.defaultParameters['rvfunc']='fourier'+str(fourierTerms)
        self.defaultParameters['rvfuncarg']=[]
        for n in range(len(popt)):
            self.defaultParameters['rvfuncarg'].append('rvfuncarg'+str(n))
            popt_temp=popt[n]
            if n%2==1:
                if adjust_x!=0:
                    popt_temp=popt_temp+(int(n/2.)+1) *2.* np.pi / period * adjust_x
                while popt_temp<-np.pi:
                    popt_temp+=2.*np.pi
                while popt_temp>np.pi:
                    popt_temp-=2.*np.pi
                popt_temp*=180./np.pi
            self.defaultParameters['rvfuncarg'+str(n)]=popt_temp
        
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
    def getLongAxis(self):#not general, more accurate is long axis is longer, wrong result if spherical
        meshfilename = self.casename+ '/'+self.meshname+ '.stl'
        pdata = lcleeHeart.readSTL(meshfilename)
        pts=[]
        for p in range(pdata.GetNumberOfPoints()):
            pts.append(pdata.GetPoint(p))
        pts=np.array(pts)
        center,boxlengthfromcenter,rotM=obb(pts)
        axis_rotated_max=np.argmax(np.sum(boxlengthfromcenter**2.,axis=1))
        pts_rotated=np.dot(pts,rotM)
        max_r=pts_rotated[:,axis_rotated_max].max()
        min_r=pts_rotated[:,axis_rotated_max].min()
        center_low,boxlengthfromcenter_low,rotM_low=obb(pts_rotated[pts_rotated[:,axis_rotated_max]<(min_r*0.9+max_r*0.1)])
        center_high,boxlengthfromcenter_high,rotM_high=obb(pts_rotated[pts_rotated[:,axis_rotated_max]>(min_r*0.1+max_r*0.9)])
        Laxis=np.zeros(3)
        if np.sum(boxlengthfromcenter_high**2.)>np.sum(boxlengthfromcenter_low**2.):
            Laxis[axis_rotated_max]=1.
        else:
            Laxis[axis_rotated_max]=-1.
        Laxis=np.dot(Laxis.reshape((1,3)),np.transpose(rotM)).reshape(-1)
        self.defaultParameters['Laxis_X']=Laxis[0]
        self.defaultParameters['Laxis_Y']=Laxis[1]
        self.defaultParameters['Laxis_Z']=Laxis[2]
        logger.info('Estimated Laxis as '+repr(Laxis))
    def generateMesh(self,Laxis=None,endo_angle='',epi_angle='',fiberSheetletAngle=0.,fiberSheetletWidth=0.,radialFiberAngle=0.,fiberLength=0.,clipratio=0.95,meshsize=0.6,saveaddstr='',toRunCountFolder=False):
        if Laxis is None:
            Laxis=np.array([np.array([self.defaultParameters['Laxis_X'],self.defaultParameters['Laxis_Y'],self.defaultParameters['Laxis_Z']])])
        if toRunCountFolder:
            if isinstance(toRunCountFolder,str):
                if toRunCountFolder[0]!='/':
                    toRunCountFolder='/'+toRunCountFolder
            else:
                toRunCountFolder='/'+str(self.runCount)
            return lcleeHeart.generateMesh(self.casename,self.meshname,Laxis,endo_angle=endo_angle,epi_angle=epi_angle,fiberSheetletAngle=fiberSheetletAngle,fiberSheetletWidth=fiberSheetletWidth,radialFiberAngle=radialFiberAngle,fiberLength=fiberLength,clipratio=clipratio,meshsize=meshsize,meshname=self.meshname+saveaddstr,saveSubFolder=toRunCountFolder)
        else:
            return lcleeHeart.generateMesh(self.casename,self.meshname,Laxis,endo_angle=endo_angle,epi_angle=epi_angle,fiberSheetletAngle=fiberSheetletAngle,fiberSheetletWidth=fiberSheetletWidth,radialFiberAngle=radialFiberAngle,fiberLength=fiberLength,clipratio=clipratio,meshsize=meshsize,meshname=self.meshname+saveaddstr)
    def runWindkessel(self,comm,tstep,dt_dt,*args):
        if self.LVage=='adult':
            #args=(p_cav,V_cav,V_art,V_ven,V_LA)
        	p_cav,V_cav,V_art,V_ven,V_LA=args
        	BCL = self.runParameters['BCL']
        	AV = self.runParameters['AV']
        	Cao = self.runParameters['Cao']
        	Cven = self.runParameters['Cven']
        	Vart0 = self.runParameters['Vart0']
        	Vven0 = self.runParameters['Vven0']
        	Rao = self.runParameters['Rao']
        	Rven = self.runParameters['Rven']
        	Rper = self.runParameters['Rper']
        	Rmv = self.runParameters['Rmv']
            
        	#### For Calculating P_LA ######################################## 
        	Ees_la = self.runParameters['Ees_la']
        	A_la = self.runParameters['A_la']
        	B_la = self.runParameters['B_la']
        	V0_la = self.runParameters['V0_la']
        	Tmax_la = self.runParameters['Tmax_la']
        	tau_la = self.runParameters['tau_la']
            
            #Time varying elastance function fro LA and RA
        	def et(t, Tmax, tau):
        		if (t <= 1.5*Tmax):
        			out = 0.5*(math.sin((math.pi/Tmax)*t - math.pi/2) + 1);
        		else:
        			out = 0.5*math.exp((-t + (1.5*Tmax))/tau);
        		return out   
        		logger.info("out="+repr(out))
        	Part = 1.0/Cao*(V_art - Vart0);
        	Pven = 1.0/Cven*(V_ven - Vven0);
        	PLV = p_cav;
            
        	if (tstep < (BCL - AV)):
        		t_la = tstep + AV;
        	else: 
        		t_la = tstep - BCL + AV;
        
        	if(fenics.MPI.rank(comm) == 0):
        		logger.info("t_LA = "+repr(t_la))
        
        		PLA = et(t_la,Tmax_la,tau_la)*Ees_la*(V_LA - V0_la) + (1 - et(t_la,Tmax_la,tau_la))*A_la*(math.exp(B_la*(V_LA - V0_la)) - 1);
        	##################################################################################################################################
        
        	if(fenics.MPI.rank(comm) == 0):
        		logger.info("P_ven = "+repr(Pven))
        		logger.info("P_LV = "+repr(PLV))
        		logger.info("P_art = "+repr(Part))	
        		logger.info("P_LA = "+repr(PLA))
        
        	#### conditions for Valves#######################################
        	if(PLV <= Part):
        		Qao = 0.0;
        	else:
        		Qao = 1.0/Rao*(PLV - Part);
        	
        	if(PLV >= PLA):
        		Qla = 0.0;
        	else: 
        		Qla = 1.0/Rmv*(PLA - PLV);
        
        	H = (Part - PLV)*0.0075;	
        	
        	Qper = 1.0/Rper*(Part - Pven);
        	Qmv = 1.0/Rven*(Pven - PLA);
        
        	if(fenics.MPI.rank(comm) == 0):
        		logger.info("Q_LA = "+repr(Qla))
        		logger.info("Q_ao = "+repr(Qao))
        		logger.info("Q_per = "+repr(Qper))
        		logger.info("Q_mv = "+repr(Qmv))
        
        	V_cav_prev = V_cav
        	V_art_prev = V_art
        	V_ven_prev = V_ven
        
        	V_cav = V_cav + dt_dt*(Qla - Qao);
        	V_art = V_art + dt_dt*(Qao - Qper);
        	V_ven = V_ven + dt_dt*(Qper - Qmv);
        	V_LA = V_LA + dt_dt*(Qmv - Qla);
                    
        	if(fenics.MPI.rank(comm) == 0):
        		logger.info("V_ven = "+repr(V_ven))
        		logger.info("V_LV = "+repr(V_cav))
        		logger.info("V_art = "+repr(V_art))
        		logger.info("V_LA = "+repr(V_LA))
                
        	if(fenics.MPI.rank(comm) == 0):
        		with open(self.casename+"/"+str(self.runCount)+"/cl.txt", "a") as fdatacl:
        		    print(tstep, p_cav*0.0075, V_cav, t_la, Pven*0.0075, PLV*0.0075, Part*0.0075, PLA*0.0075, Qla, Qao, Qper, Qmv, V_ven, V_cav, V_art, V_LA,'\n', file=fdatacl)
        	return (V_cav,V_art,V_ven,V_LA )
        else:
            lvufile=self.casename+"/"+str(self.runCount)+'/LVgenerator.txt'
            lvudata=np.loadtxt(self.casename+"/"+str(self.runCount)+"/PV_.txt")
            lvucontroldata=np.loadtxt(self.casename+'/'+self.meshname+'_lvucontrol.txt')
            lvudata2=np.concatenate((lvucontroldata[:-1],lvudata[:,2]+np.array([[lvucontroldata[-1,0],0]])),axis=0)
            np.savetxt(lvufile,lvudata2)
            ngspice_py.simLVcircuit(self.casename+'/'+self.meshname,lvucontroldata[-1,0]+tstep+dt_dt,lvufile)
            self.manualVolFile=self.casename+'/'+self.meshname+'_circuit.txt'
            V_cav=self.getVolumeFromFile(lvucontroldata[-1,0]+tstep+dt_dt)
            return V_cav
    def solveVolume(self,heart,targetVolume,voladj=0.05):
        if heart.uflforms.cavityvol()>(targetVolume):
        	tuneVol=-1
        elif heart.uflforms.cavityvol()<(targetVolume):
            tuneVol=1
        else:
            tuneVol=0
        while (heart.uflforms.cavityvol()*tuneVol)<(targetVolume*tuneVol):
        	pressadjust_voladj=max(0.01,0.75**(np.log(max(1.,abs(heart.uflforms.cavitypressure()*.0075)/10.))/np.log(2.)))
        	if (((1+tuneVol*voladj*pressadjust_voladj)*heart.uflforms.cavityvol())*tuneVol)<(targetVolume*tuneVol):
        		heart.Cavityvol.vol =heart.uflforms.cavityvol()* (1+tuneVol*voladj*pressadjust_voladj)
        	else:
        		heart.Cavityvol.vol = targetVolume
        	heart.solver.solvenonlinear()
        	p_cav = heart.uflforms.cavitypressure()
        	V_cav = heart.uflforms.cavityvol()
        	if isinstance(p_cav,(float,int)):
        		logger.info("PLV = "+repr(p_cav*.0075)+" VLV = "+repr(V_cav))
        	else:
        		raise Exception('Unable to converge to volume '+repr(targetVolume))
        	if abs(heart.uflforms.cavityvol()/targetVolume-1.)<10**-4:
        		break	
            
            
        return 1
    def solveTa(self,heart,targetTa,peaktime=None):
        if heart.t_a.t_a>(targetTa):
        	tuneTa=-1
        elif heart.t_a.t_a<(targetTa):
        	tuneTa=1
        else:
        	tuneTa=0
        while (heart.t_a.t_a*tuneTa)<(targetTa*tuneTa):
        	t=heart.t_a.t_a
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
        		heart.t_a.t_a =t+dt*tuneTa
        	else:
        		heart.t_a.t_a = targetTa
            
        	heart.dt.dt=abs(heart.t_a.t_a-t)
        	logger.info("heart.t_a.t_a = "+repr(heart.t_a.t_a)+" heart.dt.dt = "+repr(heart.dt.dt))
        	heart.solver.solvenonlinear()
        	p_cav = heart.uflforms.cavitypressure()
        	V_cav = heart.uflforms.cavityvol()
        	if isinstance(p_cav,(float,int)):
        		logger.info("PLV = "+repr(p_cav*.0075)+" VLV = "+repr(V_cav))
        	else:
        		raise Exception('Unable to converge to volume '+repr(targetTa))  
        return 1
    def solvePressure(self,heart,targetPressure,voladj=0.05,minVolumeBreak=0.,maxVolumeBreak=float('inf')):
        #target pressure in mmHg
        targetPressure=targetPressure/.0075
        if abs(heart.uflforms.cavitypressure()/targetPressure-1.)<10**-4:
            tuneVol=0
        elif heart.uflforms.cavitypressure()>(targetPressure):
        	tuneVol=-1
        elif heart.uflforms.cavitypressure()<(targetPressure):
            tuneVol=1
        else:
            tuneVol=0
        heart.Cavityvol.vol = heart.uflforms.cavityvol()
        while (heart.uflforms.cavitypressure()*tuneVol)<(targetPressure*tuneVol):
            if heart.uflforms.cavityvol()<minVolumeBreak or heart.uflforms.cavityvol()>maxVolumeBreak:
                return 0
            heart.Cavityvol.vol = heart.Cavityvol.vol*(1+tuneVol*voladj)
            heart.solver.solvenonlinear()
            p_cav = heart.uflforms.cavitypressure()
            V_cav = heart.uflforms.cavityvol()
            if isinstance(p_cav,(float,int)):
                logger.info("PLV = "+repr(p_cav*.0075)+" VLV = "+repr(V_cav))	
            else:
                raise Exception('Unable to converge to pressure '+repr(targetPressure)+' mmHg')
            if (heart.uflforms.cavitypressure()*tuneVol)>(targetPressure*tuneVol):
                voladj=voladj*0.33
                tuneVol=tuneVol*-1
                logging.debug('Change tuneVol to '+repr(tuneVol)+'  voladj='+repr(voladj))
            if abs(heart.uflforms.cavitypressure()/targetPressure-1.)<10**-4:
                break
        return 1
    def getLVbehavior(self,runParameters=None,meshname=None,meshFolder=None,volstep=50,extendVol=0.1,minESV=None,maxEDV=None,saveaddstr='',voladj=0.1,starttime=0,endtime=None,toRunCountFolder=False):
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
        elif isinstance(runParameters,str):
            if runParameters[0]!='/':
                runParameters='/'+runParameters
            runParameters=self.readRunParameters(self.casename+runParameters)
            trycontinueprevious=True
        else:
            trycontinueprevious=False
        if meshname is None:
            meshname=self.meshname
        heart=lcleeHeart.heart(self.casename+meshFolder,meshname,runParameters)
        
        
        ########################### Fenics's Newton  #########################################################
        heart.Cavityvol.vol = heart.uflforms.cavityvol()

        
        # Closed loop cycle
        BCL = runParameters['BCL']
        
        
        #######set volume solve space#############     
        EDV=runParameters['EDV_LV']
        if maxEDV is None:
            maxEDV=EDV*(1+extendVol)
        if 'ESV_LV' in runParameters:
            ESV=runParameters['ESV_LV']
        else:
            ESV=heart.Cavityvol.vol
        if minESV is None:
            minESV=ESV*(1-extendVol)
        if isinstance(volstep,(list,np.ndarray)):
            volumeSpace=np.array(volstep)
        else:
            volumeSpace=minESV*(maxEDV/minESV)**(np.linspace(0,1,num=volstep))[::-1]
        ####### Loading phase to MAX volume ####################################################
        #self.solveVolume(heart,volumeSpace[0])

        timeSpace=[]
        tstep=0
        while(tstep < 2*BCL):
            timeSpace.append(tstep)
            if(tstep >= 0.0 and tstep < 4.0):
           		tstep += 0.50
            elif (tstep >= (runParameters['t0']-0.5) and tstep < (runParameters['t0']+0.5)):
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
        if(fenics.MPI.rank(heart.comm) == 0):
            if len(timeSpace)>1:
                np.savetxt(self.casename+toRunCountFolder+"/"+meshname+"_Press_volumeSpace"+saveaddstr+".txt",np.array(volumeSpace))
                np.savetxt(self.casename+toRunCountFolder+"/"+meshname+"_Press_timeSpace"+saveaddstr+".txt",np.array(timeSpace))
        if len(timeSpace)==1:
            volPress=np.zeros((len(volumeSpace),2))
            volPress[:,0]=volumeSpace.copy()
            if abs(timeSpace[0])<=0.5:
                volPressStr='EDPVR'
            elif abs(timeSpace[0]-runParameters['t0'])<=0.5:
                volPressStr='ESPVR'
            else:
                volPressStr='t'+str(int(timeSpace[0]))+'PVR'
        logger.info('========START============')
        press_volTime_base=np.zeros(len(volumeSpace))
        press_volTime=np.zeros((len(volumeSpace),len(timeSpace)))
        startVolumeIndex=0
        
        if trycontinueprevious:
            try:
                volumeSpace2=np.loadtxt(self.casename+toRunCountFolder+"/"+meshname+"_Press_volumeSpace"+saveaddstr+".txt")
                timeSpace2=np.loadtxt(self.casename+toRunCountFolder+"/"+meshname+"_Press_timeSpace"+saveaddstr+".txt")
                if np.all(np.abs(volumeSpace2-volumeSpace)<(runParameters['EDV_LV']/1000.)) and np.all(np.abs(timeSpace2-timeSpace)<0.01):
                    press_volTime=np.loadtxt(self.casename+toRunCountFolder+"/"+meshname+"_Press_VolTime"+saveaddstr+".txt")
                    press_volTime_base=np.loadtxt(self.casename+toRunCountFolder+"/"+meshname+"_Press_VolTime_base"+saveaddstr+".txt")
                    startVolumeIndex=np.nonzero(press_volTime_base==0.)[0][0]
                    volumeSpace=volumeSpace2
                    timeSpace=timeSpace2
            except:
                pass
        for curr_volN in range(startVolumeIndex,len(volumeSpace)):
            if curr_volN==0 or len(timeSpace)>1:
                heart=lcleeHeart.heart(self.casename+meshFolder,meshname,runParameters)
            press_volTime_temp=[]
            samePressure=0
            self.solveVolume(heart,volumeSpace[curr_volN])
            for curr_timeN in range(len(timeSpace)):
            	t = timeSpace[curr_timeN]
            	logger.info('vol:'+repr(curr_volN)+'/'+repr(len(volumeSpace))+'  time:'+repr(curr_timeN)+'/'+repr(len(timeSpace))+" t_a="+repr(t))
            
            	self.solveTa(heart,t,peaktime=runParameters['t0'])
            	logger.info("PLV = "+repr(heart.uflforms.cavitypressure()*.0075)+" VLV = "+repr(heart.uflforms.cavityvol()))
            	press_volTime_temp.append(heart.uflforms.cavitypressure()*.0075)
            	if len(press_volTime_temp)>1:
                    if abs(press_volTime_temp[-1]-press_volTime_temp[-2])<10**-6:
                        samePressure+=1
                    else:
                        samePressure=0
                    if samePressure>1:
                        break
            press_volTime_base[curr_volN]=min(press_volTime_temp[-1],press_volTime_temp[0])
            press_volTime[curr_volN,:len(press_volTime_temp)]=np.array(press_volTime_temp)-press_volTime_base[curr_volN]
            if len(timeSpace)==1:
                volPress[curr_volN,0]=heart.uflforms.cavityvol()
                volPress[curr_volN,1]=press_volTime_base[curr_volN]+press_volTime[curr_volN,0]
            if(fenics.MPI.rank(heart.comm) == 0):
                if len(timeSpace)>1:
                    np.savetxt(self.casename+toRunCountFolder+"/"+meshname+"_Press_VolTime"+saveaddstr+".txt",press_volTime,header=str(curr_volN)+'solved rowVolm: '+str(volumeSpace)+'\ncolTime: '+str(timeSpace))
                    np.savetxt(self.casename+toRunCountFolder+"/"+meshname+"_Press_VolTime_base"+saveaddstr+".txt",press_volTime_base)
                else:
                    np.savetxt(self.casename+toRunCountFolder+"/"+meshname+"_"+volPressStr+saveaddstr+".txt",volPress,header="Volume, Pressure")

        if(fenics.MPI.rank(heart.comm) == 0):
            if len(timeSpace)>1:
                np.savetxt(self.casename+toRunCountFolder+"/"+meshname+"_Press_VolTime"+saveaddstr+".txt",np.array(press_volTime),header='rowVolm: '+str(volumeSpace)+'\ncolTime: '+str(timeSpace))
                np.savetxt(self.casename+toRunCountFolder+"/"+meshname+"_Press_VolTime_base"+saveaddstr+".txt",press_volTime_base)
            else:
                np.savetxt(self.casename+toRunCountFolder+"/"+meshname+"_"+volPressStr+saveaddstr+".txt",volPress,header="Volume, Pressure")
    def savehdf5(self,heart,casename,meshname,savename,saveFolder=''):
        #mesh = fenics.Mesh()
        f = fenics.HDF5File(fenics.MPI.comm_world, casename+'/'+meshname+".hdf5", 'r') 
        #f.read(mesh, meshname, False)
        facetboundaries = fenics.MeshFunction("size_t", heart.mesh, 2)
        f.read(facetboundaries, meshname+"/"+"facetboundaries")
        
        fiberFS = fenics.VectorFunctionSpace(heart.mesh, 'DG', 0)
        VQuadelem = fenics.VectorElement("Quadrature", heart.mesh.ufl_cell(), degree=4, quad_scheme="default") #segmentation fault happened due to degree = 2. degree = 4 should be implied. I changed it. (Joy) 
        VQuadelem._quad_scheme = 'default'
        for e in VQuadelem.sub_elements():
        	e._quad_scheme = 'default'
        fiberFS = fenics.FunctionSpace(heart.mesh, VQuadelem)
        f0 = fenics.Function(fiberFS)
        s0 = fenics.Function(fiberFS)
        n0 = fenics.Function(fiberFS)
        c0 = fenics.Function(fiberFS)
        l0 = fenics.Function(fiberFS)
        r0 = fenics.Function(fiberFS)
        matid = fenics.MeshFunction('size_t',heart.mesh, 3, heart.mesh.domains())
        facetboundaries = fenics.MeshFunction("size_t", heart.mesh, 2)
        edgeboundaries = fenics.MeshFunction("size_t", heart.mesh, 1)
        
        f.read(f0, meshname+"/"+"eF")
        f.read(s0, meshname+"/"+"eS")
        f.read(n0, meshname+"/"+"eN")
        f.read(c0, meshname+"/"+"eC")
        f.read(l0, meshname+"/"+"eL")
        f.read(r0, meshname+"/"+"eR")
        f.read(matid, meshname+"/"+"matid")
        f.read(facetboundaries, meshname+"/"+"facetboundaries")
        f.read(edgeboundaries, meshname+"/"+"edgeboundaries")
        
        
        matid_filename = casename+'/'+savename + "_matid.pvd"
        fenics.File(matid_filename) << matid
        f = fenics.HDF5File(heart.mesh.mpi_comm(), casename+saveFolder+'/'+savename+".hdf5", 'w')
        f.write(heart.mesh, savename)
        f.close()
       	
        f = fenics.HDF5File(heart.mesh.mpi_comm(), casename+saveFolder+'/'+savename+".hdf5", 'a') 
        f.write(facetboundaries, savename+"/"+"facetboundaries") 
        f.write(edgeboundaries, savename+"/"+"edgeboundaries") 
        f.write(matid, savename+"/"+"matid") 
        f.write(f0, savename+"/"+"eF") 
        f.write(s0, savename+"/"+"eS") 
        f.write(n0, savename+"/"+"eN")
        f.write(c0, savename+"/"+"eC") 
        f.write(l0, savename+"/"+"eL") 
        f.write(r0, savename+"/"+"eR")
       
        f.close()
       
       
        fenics.File(casename+saveFolder+'/'+savename + "_facetboundaries"+".pvd") << facetboundaries
        fenics.File(casename+saveFolder+'/'+savename + "_edgeboundaries"+".pvd") << edgeboundaries
        fenics.File(casename+saveFolder+'/'+savename + "_mesh" + ".pvd") << heart.mesh
        fenics.File(casename+saveFolder+'/'+savename + "_matid" +".pvd") << matid
        
    def adjustMeshToPressure(self,casename,meshname,targetPressure,savename=None,softAdjust=float('inf'),usePrevious=3,iterationNumber=float('inf'),initHeart=None,prevDeform=None,saveFolder='',tolerance=10.**-4.):
        #set usePrevious to False or 0 to not use previous mesh for deformation, set True to use previous mesh only and set a number to use previous mesh a number of times is it doesnt converge
        if not(isinstance(usePrevious,bool)):
            usePreviousiteration=usePrevious
            usePrevious=False
        else:
            usePreviousiteration=0
        runParameters=self.defaultParameters
        refHeart=lcleeHeart.heart(casename,meshname,runParameters)
        if initHeart is None or (usePrevious and (prevDeform is None)):
            heart=lcleeHeart.heart(casename,meshname,runParameters)
            prevDeformVectors=heart.w.sub(0).vector()[:].copy()*0.
        else:
            heart=initHeart
            if prevDeform is None:
                prevDeformVectors=heart.w.sub(0).vector()[:].copy()*0.
            else:
                prevDeformVectors=prevDeform
        targetVolume=refHeart.uflforms.cavityvol()
        
        prevMeshCoords=heart.mesh.coordinates()[:].copy()
        self.solvePressure(heart,targetPressure,voladj=0.05,maxVolumeBreak=heart.uflforms.cavityvol()*softAdjust)
        deformVectors=heart.w.sub(0).vector()[:].copy()
        maxadjustmentsq=None
        newheartVol=heart.uflforms.cavityvol()
        prevHeartVol=heart.uflforms.cavityvol()
        
        count=0
        while abs(heart.uflforms.cavitypressure()*0.0075/targetPressure-1)>tolerance or abs(heart.uflforms.cavityvol()/targetVolume-1)>tolerance:##alternate convergence criteria##(np.mean((heart.mesh.coordinates()[:]-du.vector()[:].reshape((-1,3))-refHeart.mesh.coordinates()[:])**2.).max()**1.5/refHeart.uflforms.cavityvol())<10.**-6.:
            V = fenics.VectorFunctionSpace(heart.mesh, "CG", 1)
            du = fenics.project(heart.w.sub(0),V)
            if maxadjustmentsq is None:
                maxadjustmentsq=np.sum(du.vector()[:].reshape((-1,3))**2.,axis=1).max()**1.5/newheartVol
                ratio=1.
            else:
                ratio=min(1.,maxadjustmentsq/np.sum((heart.mesh.coordinates()[:]+du.vector()[:].reshape((-1,3))-refHeart.mesh.coordinates()[:])**2.,axis=1).max()**1.5/newheartVol)
            
            
            newheartVol=heart.uflforms.cavityvol()
            logging.info('Adjusting mesh... Max adjustment ='+repr(np.sum(du.vector()[:].reshape((-1,3)),axis=1).max()))
            logging.info('ratio='+repr(ratio))
            logging.info('target volume='+repr(targetVolume))
            logging.info('prevheart volume='+repr(prevHeartVol))
            logging.info('newheart volume ='+repr(heart.uflforms.cavityvol()))
            prevHeartVol=heart.uflforms.cavityvol()
            if usePrevious:
                newDeform=heart.w.sub(0).copy(deepcopy=True)
                newDeform.vector()[:]*=-1
                newDeform.vector()[:]+=prevDeformVectors
                fenics.ALE.move(heart.mesh,newDeform)
                self.savehdf5(heart,casename,meshname,'tempadjmesh_'+meshname,saveFolder=saveFolder)
                lcleeHeart.changeFiberAngles(casename,'tempadjmesh_'+meshname,runParameters['endo_angle'],runParameters['epi_angle'],fiberSheetletAngle=runParameters['fiberSheetletAngle'],fiberSheetletWidth=runParameters['fiberSheetletWidth'],radialFiberAngle=runParameters['radialFiberAngle'],fiberLength=runParameters['lr'],saveSubFolder=saveFolder,loadSubFolder=saveFolder)
            else:
                refHeart=lcleeHeart.heart(casename,meshname,runParameters)
                newDeform=refHeart.w.sub(0).copy(deepcopy=True)
                newDeform.vector()[:]*=0.
                newDeform.vector()[:]-=heart.w.sub(0).vector()[:]
                fenics.ALE.move(refHeart.mesh,newDeform)
                self.savehdf5(refHeart,casename,meshname,'tempadjmesh_'+meshname,saveFolder=saveFolder)
                lcleeHeart.changeFiberAngles(casename,'tempadjmesh_'+meshname,runParameters['endo_angle'],runParameters['epi_angle'],fiberSheetletAngle=runParameters['fiberSheetletAngle'],fiberSheetletWidth=runParameters['fiberSheetletWidth'],radialFiberAngle=runParameters['radialFiberAngle'],fiberLength=runParameters['lr'],saveSubFolder=saveFolder,loadSubFolder=saveFolder)
            prevDeformVectors=deformVectors.copy()
            heart=lcleeHeart.heart(casename+saveFolder,'tempadjmesh_'+meshname,runParameters)
            logging.info('newheart starting volume ='+repr(heart.uflforms.cavityvol()))
            self.solvePressure(heart,targetPressure,voladj=0.05,maxVolumeBreak=newheartVol*softAdjust)
            deformVectors=heart.w.sub(0).vector()[:].copy()
            maxcoordchange=np.sum((heart.mesh.coordinates()[:]-prevMeshCoords)**2.,axis=1).max()**1.5
            prevMeshCoords=heart.mesh.coordinates()[:].copy()
            count+=1
            if (maxcoordchange/targetVolume)<10**-8:
                break
            else:
                logging.info('max coord change ratio= '+repr(maxcoordchange/targetVolume))
            if count>iterationNumber:
                logging.warning('Maximum Iteration reached for adjustMeshToPressure')
                break
        logging.info('Final heart volume ='+repr(heart.uflforms.cavityvol()))
        if savename is not None:
            self.savehdf5(heart,casename,meshname,savename,saveFolder=saveFolder)
            lcleeHeart.changeFiberAngles(casename,savename,runParameters['endo_angle'],runParameters['epi_angle'],fiberSheetletAngle=runParameters['fiberSheetletAngle'],fiberSheetletWidth=runParameters['fiberSheetletWidth'],radialFiberAngle=runParameters['radialFiberAngle'],fiberLength=runParameters['lr'],saveSubFolder=saveFolder,loadSubFolder=saveFolder)
        if not(usePrevious) and (usePreviousiteration>0) and (abs(heart.uflforms.cavitypressure()*0.0075/targetPressure-1)>10**-4. or abs(heart.uflforms.cavityvol()/targetVolume-1)>10**-4.):
            newDeform=heart.w.sub(0).copy(deepcopy=True)
            newDeform.vector()[:]*=-1
            newDeform.vector()[:]+=prevDeformVectors
            fenics.ALE.move(heart.mesh,newDeform)
            self.savehdf5(heart,casename,meshname,'tempadjmesh_'+meshname,saveFolder=saveFolder)
            lcleeHeart.changeFiberAngles(casename,'tempadjmesh_'+meshname,runParameters['endo_angle'],runParameters['epi_angle'],fiberSheetletAngle=runParameters['fiberSheetletAngle'],fiberSheetletWidth=runParameters['fiberSheetletWidth'],radialFiberAngle=runParameters['radialFiberAngle'],fiberLength=runParameters['lr'],saveSubFolder=saveFolder,loadSubFolder=saveFolder)
            heart=lcleeHeart.heart(casename+saveFolder,'tempadjmesh_'+meshname,runParameters)
            logging.info('newheart starting volume ='+repr(heart.uflforms.cavityvol()))
            logging.info('trying usePrevious = True')
            heart=self.adjustMeshToPressure(casename,meshname,targetPressure,savename=savename+'_tryusingprev',softAdjust=float('inf'),usePrevious=True,iterationNumber=usePreviousiteration,initHeart=heart,prevDeform=deformVectors,saveFolder=saveFolder)
        return heart
    def adjustMeshToPressure_tryNerror(self,casename,meshname,targetPressure,savename=None,first_tryvolume=None,iterationNumber=float('inf')):
        #use negative pressure to try and error
        runParameters=self.defaultParameters
        refHeart=lcleeHeart.heart(casename,meshname,runParameters)
        heart=lcleeHeart.heart(casename,meshname,runParameters)
        targetVolume=heart.uflforms.cavityvol()
        tryVolume_adj=0.05
        tuneVolume=-1
        if first_tryvolume is not None:
            try_negative=first_tryvolume
        else:
            try_negative=targetVolume*0.5
        self.solveVolume(heart,try_negative,voladj=0.05)
        print(np.max(np.abs(refHeart.mesh.coordinates()[:]-heart.mesh.coordinates()[:])))
        #V = fenics.VectorFunctionSpace(heart.mesh, "CG", 1)
        #du = fenics.project(heart.w.sub(0),V)
        #heart.mesh.coordinates()[:]=refHeart.mesh.coordinates()[:]+du.vector()[:].reshape((-1,3),order='C')
        #heart.mesh.smooth(20)
        fenics.ALE.move(heart.mesh,heart.w.sub(0))
        self.savehdf5(heart,casename,meshname,'tempadjmesh_'+meshname)
        lcleeHeart.changeFiberAngles(casename,'tempadjmesh_'+meshname,runParameters['endo_angle'],runParameters['epi_angle'],fiberSheetletAngle=runParameters['fiberSheetletAngle'],fiberSheetletWidth=runParameters['fiberSheetletWidth'],radialFiberAngle=runParameters['radialFiberAngle'],fiberLength=runParameters['lr'])
        heart=lcleeHeart.heart(casename,'tempadjmesh_'+meshname,runParameters)
        logger.info("PLV = "+repr(heart.uflforms.cavitypressure()*.0075)+" VLV = "+repr(heart.uflforms.cavityvol()))
        self.solvePressure(heart,targetPressure,voladj=0.05)
        if heart.uflforms.cavityvol()*tuneVolume<targetVolume*tuneVolume:
            try_negative=try_negative*(1+tuneVolume*tryVolume_adj)
        else:
            tuneVolume*=-1
            try_negative=try_negative*(1+tuneVolume*tryVolume_adj)
        count=0
        while abs(heart.uflforms.cavityvol()/targetVolume-1)>10**-4.:
            heart=lcleeHeart.heart(casename,meshname,runParameters)
            self.solveVolume(heart,try_negative,voladj=0.05)
            #V = fenics.VectorFunctionSpace(heart.mesh, "CG", 1)
            #du = fenics.project(heart.w.sub(0),V)
            #heart.mesh.coordinates()[:]=refHeart.mesh.coordinates()[:]+du.vector()[:].reshape((-1,3),order='C')
            #heart.mesh.smooth(20)
            fenics.ALE.move(heart.mesh,heart.w.sub(0))
            self.savehdf5(heart,casename,meshname,'tempadjmesh_'+meshname)
            lcleeHeart.changeFiberAngles(casename,'tempadjmesh_'+meshname,runParameters['endo_angle'],runParameters['epi_angle'],fiberSheetletAngle=runParameters['fiberSheetletAngle'],fiberSheetletWidth=runParameters['fiberSheetletWidth'],radialFiberAngle=runParameters['radialFiberAngle'],fiberLength=runParameters['lr'])
            heart=lcleeHeart.heart(casename,'tempadjmesh_'+meshname,runParameters)
        
            self.solvePressure(heart,targetPressure,voladj=0.05)#,maxVolumeBreak=heart.uflforms.cavityvol()*1.031)
            if heart.uflforms.cavityvol()*tuneVolume<targetVolume*tuneVolume:
                try_negative=try_negative*(1+tuneVolume*tryVolume_adj)
            else:
                tuneVolume*=-1
                tryVolume_adj*=0.33
                try_negative=try_negative*(1+tuneVolume*tryVolume_adj)
            count+=1
            if count>iterationNumber:
                break
        if savename is not None:
            self.savehdf5(heart,casename,meshname,savename)
            lcleeHeart.changeFiberAngles(casename,savename,runParameters['endo_angle'],runParameters['epi_angle'],fiberSheetletAngle=runParameters['fiberSheetletAngle'],fiberSheetletWidth=runParameters['fiberSheetletWidth'],radialFiberAngle=runParameters['radialFiberAngle'],fiberLength=runParameters['lr'])
        return heart
    
    def getHeartCoords(self,casename,meshname,time):
        timestep=self.timetotimestepScale*time+self.timetotimestepShift
        heart=lcleeHeart.heart(casename,meshname,self.defaultParameters)
        coords=np.copy(heart.mesh.coordinates()[:])
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
    def adjustMeshToPressureWithInverse(self,casename,meshname,targetPressure,savename=None,softAdjust=float('inf'),iterationNumber=float('inf'),initHeart=None,prevDeform=None,saveFolder='',tolerance=10.**-4.):
        #set usePrevious to False or 0 to not use previous mesh for deformation, set True to use previous mesh only and set a number to use previous mesh a number of times is it doesnt converge
        runParameters=self.defaultParameters
        refHeart=lcleeHeart.heart(casename,meshname,runParameters)
        
        targetVolume=refHeart.uflforms.cavityvol()
        logging.info('Target Volume = '+repr(targetVolume)+', target pressure = '+repr(targetPressure))
        targetMeshCoords=refHeart.mesh.coordinates()[:].copy()
        self.solvePressure(refHeart,targetPressure,voladj=0.05)
        finv=fenics.grad(refHeart.w.sub(0))+fenics.Identity(refHeart.u.geometric_dimension())
        
        heart=lcleeHeart.heart(casename,meshname,runParameters,inverseHeart=True)
        heart.invF.vector()[:]=fenics.project(finv,refHeart.invFspace).vector()[:].copy()

        heart.solver.solvenonlinear()
        deform=fenics.project(heart.w.sub(0),heart.displacementSpace)
        deform_vectors=deform.vector()[:].copy()
        fenics.ALE.move(heart.mesh,deform)
        self.savehdf5(heart,casename,meshname,'tempadjmesh_'+meshname,saveFolder=saveFolder)
        lcleeHeart.changeFiberAngles(casename,'tempadjmesh_'+meshname,runParameters['endo_angle'],runParameters['epi_angle'],fiberSheetletAngle=runParameters['fiberSheetletAngle'],fiberSheetletWidth=runParameters['fiberSheetletWidth'],radialFiberAngle=runParameters['radialFiberAngle'],fiberLength=runParameters['lr'],saveSubFolder=saveFolder,loadSubFolder=saveFolder)
        
        maxadjustmentsq=None
        newheartVol=heart.uflforms.cavityvol()
        prevHeartVol=heart.uflforms.cavityvol()
        
        refHeart=lcleeHeart.heart(casename+saveFolder,'tempadjmesh_'+meshname,runParameters)
        logging.info('New start Volume = '+repr(refHeart.uflforms.cavityvol()))
        prevMeshCoords=refHeart.mesh.coordinates()[:].copy()
        self.solvePressure(refHeart,targetPressure,voladj=0.05)
        meanSqError=pointsMeansqDiff_kabsch(refHeart.mesh.coordinates()[:],targetMeshCoords)
        logging.info('meanSqError= '+repr(meanSqError))
        count=0
        while abs(refHeart.uflforms.cavitypressure()*0.0075/targetPressure-1)>tolerance or abs(refHeart.uflforms.cavityvol()/targetVolume-1)>tolerance:##alternate convergence criteria##(np.mean((heart.mesh.coordinates()[:]-du.vector()[:].reshape((-1,3))-refHeart.mesh.coordinates()[:])**2.).max()**1.5/refHeart.uflforms.cavityvol())<10.**-6.:
            if count>=iterationNumber:
                logging.warning('Maximum Iteration reached for adjustMeshToPressure')
                break
            new_deform_vectors=fenics.project(refHeart.w.sub(0),refHeart.displacementSpace).vector()[:]
            originalHeart=lcleeHeart.heart(casename,meshname,runParameters)
            inv_deform=fenics.project(originalHeart.w.sub(0),originalHeart.displacementSpace)
            inv_deform.vector()[:]=deform_vectors+new_deform_vectors
            finv=fenics.grad(inv_deform)+fenics.Identity(originalHeart.u.geometric_dimension())
            
            
            heart=lcleeHeart.heart(casename+saveFolder,'tempadjmesh_'+meshname,runParameters,inverseHeart=True)
            heart.invF.vector()[:]=fenics.project(finv,originalHeart.invFspace).vector()[:].copy()
            heart.solver.solvenonlinear()
            deform=fenics.project(heart.w.sub(0),heart.displacementSpace)
            deform_vectors=deform_vectors+deform.vector()[:].copy()
            
            deform=fenics.project(originalHeart.w.sub(0),originalHeart.displacementSpace)
            deform.vector()[:]=deform_vectors.copy()
            fenics.ALE.move(originalHeart.mesh,deform)
            self.savehdf5(originalHeart,casename,meshname,'tempadjmesh_'+meshname,saveFolder=saveFolder)
            lcleeHeart.changeFiberAngles(casename,'tempadjmesh_'+meshname,runParameters['endo_angle'],runParameters['epi_angle'],fiberSheetletAngle=runParameters['fiberSheetletAngle'],fiberSheetletWidth=runParameters['fiberSheetletWidth'],radialFiberAngle=runParameters['radialFiberAngle'],fiberLength=runParameters['lr'],saveSubFolder=saveFolder,loadSubFolder=saveFolder)
            
            refHeart=lcleeHeart.heart(casename+saveFolder,'tempadjmesh_'+meshname,runParameters)
            logging.info('New start Volume = '+repr(refHeart.uflforms.cavityvol()))
            self.solvePressure(refHeart,targetPressure,voladj=0.05)
            meanSqError=pointsMeansqDiff_kabsch(refHeart.mesh.coordinates()[:],targetMeshCoords)
            logging.info('meanSqError= '+repr(meanSqError))
            maxcoordchange=pointsMeansqDiff_kabsch(refHeart.mesh.coordinates()[:],prevMeshCoords,averageFunc=np.max)**1.5
            prevMeshCoords=refHeart.mesh.coordinates()[:].copy()
            count+=1
            if (maxcoordchange/targetVolume)<10**-8:
                break
            else:
                logging.info('max coord change ratio= '+repr(maxcoordchange/targetVolume))
        logging.info('Final heart volume ='+repr(refHeart.uflforms.cavityvol()))
        if savename is not None:
            self.savehdf5(refHeart,casename,meshname,savename,saveFolder=saveFolder)
            lcleeHeart.changeFiberAngles(casename,savename,runParameters['endo_angle'],runParameters['epi_angle'],fiberSheetletAngle=runParameters['fiberSheetletAngle'],fiberSheetletWidth=runParameters['fiberSheetletWidth'],radialFiberAngle=runParameters['radialFiberAngle'],fiberLength=runParameters['lr'],saveSubFolder=saveFolder,loadSubFolder=saveFolder)
        return refHeart
    def getUnloadedGeometry(self,editParameters=None,casename=None,meshname=None,targetPressure=None,targetVolume=None,tryPressureRatio=0.2,targetMeshCoords=None,savename=None,iterationNumber=0,toRunCountFolder=False,inverseHeart=True):
        #when targetMeshCoords is not None, set target pressure as the maximum pressure to try
        #= mesh.coordinates()[:] at target volume
        runParameters=self.defaultParameters.copy()
        if editParameters is not None:
            for key in runParameters:
                if key in editParameters:
                    runParameters[key]=editParameters[key]
        else:
            editParameters={}
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
        if targetPressure is None:
            targetPressure=runParameters['EDP_LV']
        if targetVolume is None:
            targetVolume=runParameters['EDV_LV']
        if targetMeshCoords is not None:
            heart=lcleeHeart.heart(casename,meshname,runParameters)
            self.solveVolume(heart,targetVolume,voladj=0.05)
            fenics.ALE.move(heart.mesh,heart.w.sub(0))
            meshRMS_error=np.array([pointsMeansqDiff_kabsch(heart.mesh.coordinates()[:],targetMeshCoords),float('inf')])
            #meshRMS_error=np.array([pointsMeansqDiff_kabsch(heart.mesh.coordinates()[:],targetMeshCoords)*0.5,float('inf')])
            meshRMS_pressurebound=np.array([0,targetPressure])
            tryPressure=targetPressure
            adjtryPressure=targetPressure
        else:
            tryPressure=tryPressureRatio*targetPressure
            adjtryPressure=min(tryPressureRatio*0.33,(1-tryPressureRatio)*0.33,0.1)
        if inverseHeart:
            heart=self.adjustMeshToPressureWithInverse(casename,meshname,tryPressure,savename=savename,iterationNumber=iterationNumber,saveFolder=toRunCountFolder)
        else:
            heart=self.adjustMeshToPressure(casename,meshname,tryPressure,savename=savename,usePrevious=0,iterationNumber=iterationNumber,saveFolder=toRunCountFolder)
        self.savehdf5(heart,casename,meshname,'tempadjmesh_'+meshname,saveFolder=toRunCountFolder)
        lcleeHeart.changeFiberAngles(casename,'tempadjmesh_'+meshname,runParameters['endo_angle'],runParameters['epi_angle'],fiberSheetletAngle=runParameters['fiberSheetletAngle'],fiberSheetletWidth=runParameters['fiberSheetletWidth'],radialFiberAngle=runParameters['radialFiberAngle'],fiberLength=runParameters['lr'],saveSubFolder=toRunCountFolder,loadSubFolder=toRunCountFolder)
        
        heart=lcleeHeart.heart(casename+toRunCountFolder,'tempadjmesh_'+meshname,runParameters)
        self.solveVolume(heart,targetVolume,voladj=0.05)
        
        
            
        if heart.uflforms.cavitypressure()*0.0075<targetPressure:
            tunePressure=1
        else:
            tunePressure=-1
        while abs(heart.uflforms.cavitypressure()*0.0075/targetPressure-1)>10**-4.:
            if targetMeshCoords is not None:
                fenics.ALE.move(heart.mesh,heart.w.sub(0))
                temperror=pointsMeansqDiff_kabsch(heart.mesh.coordinates()[:],targetMeshCoords)
                #temperror=pointsMeansqDiff_kabsch(heart.mesh.coordinates()[:],targetMeshCoords)*0.5
                #heart=lcleeHeart.heart(casename,'tempadjmesh_'+meshname,runParameters)
                #self.solveVolume(heart,initialVolume,voladj=0.05)
                #fenics.ALE.move(heart.mesh,heart.w.sub(0))
                #temperror+=pointsMeansqDiff_kabsch(heart.mesh.coordinates()[:],initialMeshCoordinates)*0.5
                replaceInd=np.argmax(meshRMS_error)
                meshRMS_error[replaceInd]=temperror
                meshRMS_pressurebound[replaceInd]=tryPressure
                replaceInd=np.argmax(meshRMS_error)
                tryPressure=meshRMS_pressurebound[replaceInd]*2./3.+meshRMS_pressurebound[1-replaceInd]/3.
                adjtryPressure=meshRMS_pressurebound[1]-meshRMS_pressurebound[0]
            else:
                tryPressureRatio=tryPressureRatio+tunePressure*adjtryPressure
                tryPressure=tryPressureRatio*targetPressure
                logging.info('Trying Pressure '+repr(tryPressure))
            if inverseHeart:
                heart=self.adjustMeshToPressureWithInverse(casename,meshname,tryPressure,savename=savename,iterationNumber=iterationNumber,saveFolder=toRunCountFolder)
            else:
                heart=self.adjustMeshToPressure(casename,meshname,tryPressure,savename=savename,usePrevious=0,iterationNumber=iterationNumber,saveFolder=toRunCountFolder,tolerance=10.**-3)
            self.savehdf5(heart,casename,meshname,'tempadjmesh_'+meshname,saveFolder=toRunCountFolder)
            lcleeHeart.changeFiberAngles(casename,'tempadjmesh_'+meshname,runParameters['endo_angle'],runParameters['epi_angle'],fiberSheetletAngle=runParameters['fiberSheetletAngle'],fiberSheetletWidth=runParameters['fiberSheetletWidth'],radialFiberAngle=runParameters['radialFiberAngle'],fiberLength=runParameters['lr'],saveSubFolder=toRunCountFolder,loadSubFolder=toRunCountFolder)
            heart=lcleeHeart.heart(casename+toRunCountFolder,'tempadjmesh_'+meshname,runParameters)
            self.solveVolume(heart,targetVolume,voladj=0.05)
            if targetMeshCoords is not None:
                logging.info('Trying Pressure '+repr(tryPressure)+' between '+repr(meshRMS_pressurebound))
                logging.info('Errors '+repr(meshRMS_error))
            elif (tunePressure*heart.uflforms.cavitypressure()*0.0075)>(tunePressure*targetPressure):
                tunePressure*=-1
                adjtryPressure*=0.33
            elif tunePressure==-1 and (tryPressureRatio<1.5*adjtryPressure):
                adjtryPressure*=0.33
            if tryPressure<10.**-4. and tunePressure==-1:
                break
            if adjtryPressure<10**-3:
                break
        self.savehdf5(heart,casename,meshname,savename,saveFolder=toRunCountFolder)
        lcleeHeart.changeFiberAngles(casename,savename,runParameters['endo_angle'],runParameters['epi_angle'],fiberSheetletAngle=runParameters['fiberSheetletAngle'],fiberSheetletWidth=runParameters['fiberSheetletWidth'],radialFiberAngle=runParameters['radialFiberAngle'],fiberLength=runParameters['lr'],saveSubFolder=toRunCountFolder,loadSubFolder=toRunCountFolder)
        return heart
    def getUnloadedGeometry(self,editParameters=None,casename=None,meshname=None,targetPressure=None,targetVolume=None,tryPressureRatio=0.2,targetMeshCoords=None,savename=None,iterationNumber=0,toRunCountFolder=False,inverseHeart=True):
        #when targetMeshCoords is not None, set target pressure as the maximum pressure to try
        #= mesh.coordinates()[:] at target volume
        runParameters=self.defaultParameters.copy()
        if editParameters is not None:
            for key in runParameters:
                if key in editParameters:
                    runParameters[key]=editParameters[key]
        else:
            editParameters={}
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
        if targetPressure is None:
            targetPressure=runParameters['EDP_LV']
        if targetVolume is None:
            targetVolume=runParameters['EDV_LV']
        if targetMeshCoords is not None:
            heart=lcleeHeart.heart(casename,meshname,runParameters)
            self.solveVolume(heart,targetVolume,voladj=0.05)
            fenics.ALE.move(heart.mesh,heart.w.sub(0))
            meshRMS_error=np.array([pointsMeansqDiff_kabsch(heart.mesh.coordinates()[:],targetMeshCoords),float('inf')])
            #meshRMS_error=np.array([pointsMeansqDiff_kabsch(heart.mesh.coordinates()[:],targetMeshCoords)*0.5,float('inf')])
            meshRMS_pressurebound=np.array([0,targetPressure])
            tryPressure=targetPressure
            adjtryPressure=targetPressure
        else:
            tryPressure=tryPressureRatio*targetPressure
            adjtryPressure=min(tryPressureRatio*0.33,(1-tryPressureRatio)*0.33,0.1)
        if inverseHeart:
            heart=self.adjustMeshToPressureWithInverse(casename,meshname,tryPressure,savename=savename,iterationNumber=iterationNumber,saveFolder=toRunCountFolder)
        else:
            heart=self.adjustMeshToPressure(casename,meshname,tryPressure,savename=savename,usePrevious=0,iterationNumber=iterationNumber,saveFolder=toRunCountFolder)
        self.savehdf5(heart,casename,meshname,'tempadjmesh_'+meshname,saveFolder=toRunCountFolder)
        lcleeHeart.changeFiberAngles(casename,'tempadjmesh_'+meshname,runParameters['endo_angle'],runParameters['epi_angle'],fiberSheetletAngle=runParameters['fiberSheetletAngle'],fiberSheetletWidth=runParameters['fiberSheetletWidth'],radialFiberAngle=runParameters['radialFiberAngle'],fiberLength=runParameters['lr'],saveSubFolder=toRunCountFolder,loadSubFolder=toRunCountFolder)
        
        heart=lcleeHeart.heart(casename+toRunCountFolder,'tempadjmesh_'+meshname,runParameters)
        self.solveVolume(heart,targetVolume,voladj=0.05)
        
        
            
        if heart.uflforms.cavitypressure()*0.0075<targetPressure:
            tunePressure=1
        else:
            tunePressure=-1
        while abs(heart.uflforms.cavitypressure()*0.0075/targetPressure-1)>10**-4.:
            if targetMeshCoords is not None:
                fenics.ALE.move(heart.mesh,heart.w.sub(0))
                temperror=pointsMeansqDiff_kabsch(heart.mesh.coordinates()[:],targetMeshCoords)
                #temperror=pointsMeansqDiff_kabsch(heart.mesh.coordinates()[:],targetMeshCoords)*0.5
                #heart=lcleeHeart.heart(casename,'tempadjmesh_'+meshname,runParameters)
                #self.solveVolume(heart,initialVolume,voladj=0.05)
                #fenics.ALE.move(heart.mesh,heart.w.sub(0))
                #temperror+=pointsMeansqDiff_kabsch(heart.mesh.coordinates()[:],initialMeshCoordinates)*0.5
                replaceInd=np.argmax(meshRMS_error)
                meshRMS_error[replaceInd]=temperror
                meshRMS_pressurebound[replaceInd]=tryPressure
                replaceInd=np.argmax(meshRMS_error)
                tryPressure=meshRMS_pressurebound[replaceInd]*2./3.+meshRMS_pressurebound[1-replaceInd]/3.
                adjtryPressure=meshRMS_pressurebound[1]-meshRMS_pressurebound[0]
            else:
                tryPressureRatio=tryPressureRatio+tunePressure*adjtryPressure
                tryPressure=tryPressureRatio*targetPressure
                logging.info('Trying Pressure '+repr(tryPressure))
            if inverseHeart:
                heart=self.adjustMeshToPressureWithInverse(casename,meshname,tryPressure,savename=savename,iterationNumber=iterationNumber,saveFolder=toRunCountFolder)
            else:
                heart=self.adjustMeshToPressure(casename,meshname,tryPressure,savename=savename,usePrevious=0,iterationNumber=iterationNumber,saveFolder=toRunCountFolder,tolerance=10.**-3)
            self.savehdf5(heart,casename,meshname,'tempadjmesh_'+meshname,saveFolder=toRunCountFolder)
            lcleeHeart.changeFiberAngles(casename,'tempadjmesh_'+meshname,runParameters['endo_angle'],runParameters['epi_angle'],fiberSheetletAngle=runParameters['fiberSheetletAngle'],fiberSheetletWidth=runParameters['fiberSheetletWidth'],radialFiberAngle=runParameters['radialFiberAngle'],fiberLength=runParameters['lr'],saveSubFolder=toRunCountFolder,loadSubFolder=toRunCountFolder)
            heart=lcleeHeart.heart(casename+toRunCountFolder,'tempadjmesh_'+meshname,runParameters)
            self.solveVolume(heart,targetVolume,voladj=0.05)
            if targetMeshCoords is not None:
                logging.info('Trying Pressure '+repr(tryPressure)+' between '+repr(meshRMS_pressurebound))
                logging.info('Errors '+repr(meshRMS_error))
            elif (tunePressure*heart.uflforms.cavitypressure()*0.0075)>(tunePressure*targetPressure):
                tunePressure*=-1
                adjtryPressure*=0.33
            elif tunePressure==-1 and (tryPressureRatio<1.5*adjtryPressure):
                adjtryPressure*=0.33
            if tryPressure<10.**-4. and tunePressure==-1:
                break
            if adjtryPressure<10**-3:
                break
        self.savehdf5(heart,casename,meshname,savename,saveFolder=toRunCountFolder)
        lcleeHeart.changeFiberAngles(casename,savename,runParameters['endo_angle'],runParameters['epi_angle'],fiberSheetletAngle=runParameters['fiberSheetletAngle'],fiberSheetletWidth=runParameters['fiberSheetletWidth'],radialFiberAngle=runParameters['radialFiberAngle'],fiberLength=runParameters['lr'],saveSubFolder=toRunCountFolder,loadSubFolder=toRunCountFolder)
        return heart
    def writeRunParameters(self,folder,runParameters):
        os.makedirs(folder,exist_ok=True)
        with open(folder+"/runParameters.txt", "w") as f:
            for key, val in runParameters.items():
                f.write("{0:s} , {1:s}\n".format(key,str(val)))
    def decodeStringValue(self,strVal):
        if strVal[-1]=='\n':
            strVal=strVal[:-1]
        while strVal[-1]==' ':
            strVal=strVal[:-1]
        while strVal[0]==' ':
            strVal=strVal[1:]
        if strVal == 'None':
            return None
        elif strVal[0]=='[' and strVal[0]==']':
            temp=strVal.split(', ')
            for n in range(len(temp)):
                temp[n]=self.decodeStringValue(temp[n])
            return temp
        else:
            try:
                temp=int(strVal)
            except:
                try:
                    temp=float(strVal)
                except:
                    temp=strVal
            return temp
    def readRunParameters(self,folder):
        runParameters={}
        with open(folder+"/runParameters.txt", "r") as f:
            lines=f.readlines()
        for line in lines:
            if line[0]=='#':
                continue
            temp=line.split(sep=' , ')
            runParameters[temp[0]]=self.decodeStringValue(temp[1])
        return runParameters
    def getLastcycleCircuitResults(self,BCL):
        cir_results=np.loadtxt(self.casename+"/"+str(self.runCount)+'/'+'circuit_results.txt',skiprows=1)
        with open(self.casename+"/"+str(self.runCount)+'/'+'circuit_results.txt','r') as f:
            cir_results_header=f.readline()
        if cir_results_header[-1]=='\n':
            cir_results_header=cir_results_header[:-1]
        cir_results_header=cir_results_header[1:]
        cir_results[:,0]*=1000. #set to ms
        while cir_results[-1,0]>=(2.*BCL):
            cir_results[:,0]-=BCL
        cir_results=cir_results[cir_results[:,0]>=0]
        cir_results=cir_results[cir_results[:,0]<BCL]
        np.savetxt(self.casename+"/"+str(self.runCount)+'/'+'circuit_results_lastcycle.txt',cir_results,header=cir_results_header,comments='#')
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
        elif self.defaultRunMode=='unloadedGeometryRun':
            runParameters=self.unloadedGeometryRun(editParameters=editParameters,**run_kwargs)
        elif self.defaultRunMode=='fullWindkesselRun':
            runParameters=self.fullWindkesselRun(editParameters=editParameters,**run_kwargs)
        elif self.defaultRunMode=='flowrateWindkesselRun':
            runParameters=self.flowrateWindkesselRun(editParameters=editParameters,**run_kwargs)
        return runParameters
    def getEDPVR(self,editParameters=None,meshname=None,meshFolder=None,volstep=50,extendVol=0.1,minESV=None,maxEDV=None,saveaddstr='',voladj=0.1,toRunCountFolder=False):
        runParameters=self.defaultParameters.copy()
        if editParameters is not None:
            for key in runParameters:
                if key in editParameters:
                    runParameters[key]=editParameters[key]
        else:
            editParameters={}
        self.getLVbehavior(runParameters=runParameters,meshname=meshname,meshFolder=meshFolder,volstep=volstep,extendVol=extendVol,minESV=minESV,maxEDV=maxEDV,saveaddstr=saveaddstr,voladj=voladj,starttime=0,endtime=0,toRunCountFolder=toRunCountFolder)
    def getESPVR(self,editParameters=None,meshname=None,meshFolder=None,volstep=50,extendVol=0.1,minESV=None,maxEDV=None,saveaddstr='',voladj=0.1,toRunCountFolder=False):
        runParameters=self.defaultParameters.copy()
        if editParameters is not None:
            for key in runParameters:
                if key in editParameters:
                    runParameters[key]=editParameters[key]
        else:
            editParameters={}
        self.getLVbehavior(runParameters=runParameters,meshname=meshname,meshFolder=meshFolder,volstep=volstep,extendVol=extendVol,minESV=minESV,maxEDV=maxEDV,saveaddstr=saveaddstr,voladj=voladj,starttime=float(int(runParameters['t0'])),endtime=float(int(runParameters['t0'])),toRunCountFolder=toRunCountFolder)
    def unloadedGeometryRun(self,editParameters=None,inverseHeart=True,tryPressureRatio=0.2):
        runParameters=self.defaultParameters.copy()
        if editParameters is not None:
            for key in runParameters:
                if key in editParameters:
                    runParameters[key]=editParameters[key]
        else:
            editParameters={}
        self.writeRunParameters(self.casename+"/"+str(self.runCount),runParameters)
        #if self.runCount<=4:
        #    return runParameters
        self.generateMesh(Laxis=np.array([runParameters['Laxis_X'],runParameters['Laxis_Y'],runParameters['Laxis_Z']]),endo_angle=runParameters['endo_angle'],epi_angle=runParameters['epi_angle'],fiberSheetletAngle=runParameters['fiberSheetletAngle'],fiberSheetletWidth=runParameters['fiberSheetletWidth'],radialFiberAngle=runParameters['radialFiberAngle'],fiberLength=runParameters['lr'],clipratio=runParameters['clip_ratio'],toRunCountFolder=True)
        heartRef=lcleeHeart.heart(self.casename+"/"+str(self.runCount),self.meshname,runParameters)
        volRef=heartRef.uflforms.cavityvol()
        heart=self.getUnloadedGeometry(editParameters=runParameters,casename=self.casename+"/"+str(self.runCount),savename=self.meshname+'_unloadedmesh',toRunCountFolder=False,inverseHeart=inverseHeart,tryPressureRatio=tryPressureRatio)
        self.solveVolume(heart,runParameters['EDV_LV'],voladj=0.05)
        fenics.ALE.move(heart.mesh,heart.w.sub(0))
        np.savetxt(self.casename+"/"+str(self.runCount)+"/coordsDiaMesh.txt",heart.mesh.coordinates()[:])
        self.savehdf5(heart,self.casename+"/"+str(self.runCount),self.meshname+'_unloadedmesh',self.meshname+'_unloadedmesh_atDiastole')
        lcleeHeart.changeFiberAngles(self.casename+"/"+str(self.runCount),self.meshname+'_unloadedmesh_atDiastole',runParameters['endo_angle'],runParameters['epi_angle'],fiberSheetletAngle=runParameters['fiberSheetletAngle'],fiberSheetletWidth=runParameters['fiberSheetletWidth'],radialFiberAngle=runParameters['radialFiberAngle'],fiberLength=runParameters['lr'])
        heart=lcleeHeart.heart(self.casename+"/"+str(self.runCount),self.meshname+'_unloadedmesh',runParameters)
        self.solveVolume(heart,volRef,voladj=0.05)
        fenics.ALE.move(heart.mesh,heart.w.sub(0))
        self.savehdf5(heart,self.casename+"/"+str(self.runCount),self.meshname+'_unloadedmesh',self.meshname+'_unloadedmesh_at'+self.meshname)
        lcleeHeart.changeFiberAngles(self.casename+"/"+str(self.runCount),self.meshname+'_unloadedmesh_at'+self.meshname,runParameters['endo_angle'],runParameters['epi_angle'],fiberSheetletAngle=runParameters['fiberSheetletAngle'],fiberSheetletWidth=runParameters['fiberSheetletWidth'],radialFiberAngle=runParameters['radialFiberAngle'],fiberLength=runParameters['lr'])
        return runParameters
    def LVbehaviorRun(self,editParameters=None,unloadGeo=True,folderToLVbehavior=None,runCycles=10):
        #if using displacement to get unloaded geometry, set unloadGeo='displacement'
        if isinstance(folderToLVbehavior,str):
            if folderToLVbehavior[0]!='/':
                folderToLVbehavior='/'+folderToLVbehavior
        runParameters=self.defaultParameters.copy()
        if editParameters is not None:
            for key in runParameters:
                if key in editParameters:
                    runParameters[key]=editParameters[key]
        else:
            editParameters={}
        self.writeRunParameters(self.casename+"/"+str(self.runCount),runParameters)
        if min(self.defaultParameters.getParameterRelation(editParameters.keys())+[2])<1:
            self.generateMesh(Laxis=np.array([runParameters['Laxis_X'],runParameters['Laxis_Y'],runParameters['Laxis_Z']]),endo_angle=runParameters['endo_angle'],epi_angle=runParameters['epi_angle'],fiberSheetletAngle=runParameters['fiberSheetletAngle'],fiberSheetletWidth=runParameters['fiberSheetletWidth'],radialFiberAngle=runParameters['radialFiberAngle'],fiberLength=runParameters['lr'],clipratio=runParameters['clip_ratio'],toRunCountFolder=True)
        elif not(os.path.isfile(self.casename+'/'+self.meshname+'.hdf5')):
            if not((folderToLVbehavior is not None and min(self.defaultParameters.getParameterRelation(editParameters.keys())+[2])>0) and os.path.isfile(self.casename+folderToLVbehavior+'/'+self.meshname+'_unloadedmesh.hdf5')):
                self.generateMesh(Laxis=np.array([runParameters['Laxis_X'],runParameters['Laxis_Y'],runParameters['Laxis_Z']]),endo_angle=runParameters['endo_angle'],epi_angle=runParameters['epi_angle'],fiberSheetletAngle=runParameters['fiberSheetletAngle'],fiberSheetletWidth=runParameters['fiberSheetletWidth'],radialFiberAngle=runParameters['radialFiberAngle'],fiberLength=runParameters['lr'],clipratio=runParameters['clip_ratio'],toRunCountFolder=False)
        if unloadGeo:
            meshname=self.meshname+'_unloadedmesh'
        else:
            meshname=self.meshname
        
        if unloadGeo:
            if unloadGeo=='displacement':
                inverseHeart=False
            else:
                inverseHeart=True
            if folderToLVbehavior is not None and min(self.defaultParameters.getParameterRelation(editParameters.keys())+[2])>0:
                if not(os.path.isfile(self.casename+folderToLVbehavior+'/'+self.meshname+'_unloadedmesh.hdf5')):
                    self.getUnloadedGeometry(editParameters=runParameters,savename=self.meshname+'_unloadedmesh',targetPressure=runParameters['EDP_LV'],targetVolume=runParameters['EDV_LV'],toRunCountFolder=folderToLVbehavior,inverseHeart=inverseHeart)
            elif min(self.defaultParameters.getParameterRelation(editParameters.keys())+[2])<1:
                self.getUnloadedGeometry(editParameters=runParameters,casename=self.casename+'/'+str(self.runCount),targetPressure=runParameters['EDP_LV'],targetVolume=runParameters['EDV_LV'],savename=self.meshname+'_unloadedmesh',inverseHeart=inverseHeart)
            else:
                self.getUnloadedGeometry(editParameters=runParameters,savename=self.meshname+'_unloadedmesh',targetPressure=runParameters['EDP_LV'],targetVolume=runParameters['EDV_LV'],toRunCountFolder=True,inverseHeart=inverseHeart)
        volstep_ej=(runParameters['EDV_LV']-runParameters['ESV_LV'])/runParameters['EDV_LV']
        if volstep_ej<0.6:
            volstep_esv=0.4*runParameters['EDV_LV']
        else:
            volstep_esv=runParameters['ESV_LV']
        volstep=volstep_esv*0.9*(runParameters['EDV_LV']/volstep_esv*1.1/0.9)**(np.linspace(0,1,num=50))[::-1]
        volstep=np.concatenate((volstep[:1]*1.3,volstep[:1]*1.1,volstep,volstep[-1:]*0.9,volstep[-1:]*0.7),axis=0)
        if folderToLVbehavior is None or min(self.defaultParameters.getParameterRelation(editParameters.keys())+[2])<=1:
            if min(self.defaultParameters.getParameterRelation(editParameters.keys())+[2])<1:
                self.getLVbehavior(runParameters=runParameters,meshname=meshname,meshFolder=str(self.runCount),volstep=volstep,toRunCountFolder=True)
            elif unloadGeo:
                self.getLVbehavior(runParameters=runParameters,meshname=meshname,meshFolder=str(self.runCount),volstep=volstep,toRunCountFolder=True)
            else:
                self.getLVbehavior(runParameters=runParameters,meshname=meshname,volstep=volstep,toRunCountFolder=True)
        elif not(os.path.isfile(self.casename+folderToLVbehavior+'/'+meshname+"_Press_VolTime.txt")):
            if unloadGeo:
                self.getLVbehavior(runParameters=runParameters,meshname=meshname,meshFolder=folderToLVbehavior,volstep=volstep,toRunCountFolder=folderToLVbehavior)
            else:
                self.getLVbehavior(runParameters=runParameters,meshname=meshname,volstep=volstep,toRunCountFolder=folderToLVbehavior)
        if folderToLVbehavior is None:
            ngspice_py.generateLVtable(self.casename+"/"+str(self.runCount)+'/'+meshname,runParameters['BCL'],timetopeak=runParameters['t0'],loading_casename=self.casename+"/"+str(self.runCount)+'/'+meshname,scale_T0_LV=runParameters['Windkessel_scale_T0_LV'])
        else:
            ngspice_py.generateLVtable(self.casename+"/"+str(self.runCount)+'/'+meshname,runParameters['BCL'],timetopeak=runParameters['t0'],loading_casename=self.casename+folderToLVbehavior+'/'+meshname,scale_T0_LV=runParameters['Windkessel_scale_T0_LV'])
        cmd = "cp "+self.casename+'/'+self.meshname+'_rvflowrate.txt'+" " + self.casename+"/"+str(self.runCount)+'/'+meshname+'_rvflowrate.txt'
        os.system(cmd)
        
        if runParameters['ES_time'] is None:
            ngspice_py.createLVcircuit(self.casename+"/"+str(self.runCount)+'/'+meshname,runParameters)
        else:
            ngspice_py.createLVcircuit(self.casename+"/"+str(self.runCount)+'/'+meshname,runParameters,skipVariableList=["timetopeaktension"])
        ngspice_py.simLVcircuit_align_EDvol_and_EStime(self.casename+"/"+str(self.runCount)+'/'+meshname,runParameters['BCL']*runCycles,self.casename+"/"+str(self.runCount)+'/'+meshname+'_lvcirtable.txt',runParameters['BCL'],runParameters['ES_time'],runParameters['t0'],lvinputvar='table2dtrackrelaxphase',initLVvol=runParameters['EDV_LV'],vla0=runParameters['vla0'],vra0=runParameters['vra0'])
        self.getLastcycleCircuitResults(runParameters['BCL'])
        return runParameters
    def flowrateWindkesselRun(self,editParameters=None,lvufile=None,lvinputvar=None,stepTime=10):
        #if using displacement to get unloaded geometry, set unloadGeo='displacement'
        runParameters=self.defaultParameters.copy()
        if editParameters is not None:
            for key in runParameters:
                if key in editParameters:
                    runParameters[key]=editParameters[key]
        else:
            editParameters={}
        self.writeRunParameters(self.casename+"/"+str(self.runCount),runParameters)
        if lvufile is None:
            lvufile=self.casename+"/"+str(self.runCount)+'/'+self.meshname+'_lvu_outflowrate.txt'
            lvinputvar='I'
        if not(os.path.isfile(lvufile)):
            lvinputvar='I'
            ref_time=np.arange(-1,self.defaultParameters['BCL']*10.1+1)
            ref_volume=self.manualVol(ref_time)
            ref_outflowrate=(ref_volume[2:]-ref_volume[:-2])/2.*1000.
            tempdata=np.array([ref_time[1:-1]/1000.,ref_outflowrate]).T
            np.savetxt(lvufile,tempdata)
        ngspice_py.createLVcircuit(self.casename+"/"+str(self.runCount)+'/'+self.meshname,runParameters,stepTime=stepTime)
        ngspice_py.simLVcircuit(self.casename+"/"+str(self.runCount)+'/'+self.meshname,runParameters['BCL']*10,lvufile,lvinputvar=lvinputvar,initLVvol=runParameters['EDV_LV'],vla0=runParameters['vla0'],vra0=runParameters['vra0'])
        self.getLastcycleCircuitResults(runParameters['BCL'])
        return runParameters
    def fullWindkesselRun(self,editParameters=None,unloadGeo=True,folderToLVbehavior=None,runCycles=10):
        #if using displacement to get unloaded geometry, set unloadGeo='displacement'
        if isinstance(folderToLVbehavior,str):
            if folderToLVbehavior[0]!='/':
                folderToLVbehavior='/'+folderToLVbehavior
        runParameters=self.defaultParameters.copy()
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
            ngspice_py.generateLVtable(self.casename+"/"+str(self.runCount)+'/'+meshname,runParameters['BCL'],timetopeak=runParameters['t0'],loading_casename=self.casename+"/"+str(self.runCount)+'/'+meshname,scale_T0_LV=runParameters['Windkessel_scale_T0_LV'])
        else:
            ngspice_py.generateLVtable(self.casename+"/"+str(self.runCount)+'/'+meshname,runParameters['BCL'],timetopeak=runParameters['t0'],loading_casename=self.casename+folderToLVbehavior+'/'+meshname,scale_T0_LV=runParameters['Windkessel_scale_T0_LV'])
        cmd = "cp "+self.casename+'/'+self.meshname+'_rvflowrate.txt'+" " + self.casename+"/"+str(self.runCount)+'/'+meshname+'_rvflowrate.txt'
        os.system(cmd)
        
        if runParameters['ES_time'] is None:
            ngspice_py.createLVcircuit(self.casename+"/"+str(self.runCount)+'/'+meshname,runParameters)
        else:
            ngspice_py.createLVcircuit(self.casename+"/"+str(self.runCount)+'/'+meshname,runParameters,skipVariableList=["timetopeaktension"])
        ngspice_py.simLVcircuit_align_EDvol_and_EStime(self.casename+"/"+str(self.runCount)+'/'+meshname,runParameters['BCL']*runCycles,self.casename+"/"+str(self.runCount)+'/'+meshname+'_lvcirtable.txt',runParameters['BCL'],runParameters['ES_time'],runParameters['t0'],lvinputvar='table2dtrackrelaxphase',initLVvol=runParameters['EDV_LV'],vla0=runParameters['vla0'],vra0=runParameters['vra0'])
        self.getLastcycleCircuitResults(runParameters['BCL'])
        return runParameters
    def iterativeRun(self,editParameters=None,runTimeList=None,endTime=None,manualPhaseTimeInput=False,unloadGeo=False,setHeart=None,outputResultList=None,tryPressureRatio=0.2):
        '''
        setHeart: list or lcleeHeart object: list [folder name, meshname]
        '''
        runParameters=self.defaultParameters.copy()
        if editParameters is not None:
            for key in runParameters:
                if key in editParameters:
                    runParameters[key]=editParameters[key]
        else:
            editParameters={}
        self.writeRunParameters(self.casename+"/"+str(self.runCount),runParameters)
        
        if setHeart is None:
            if min(self.defaultParameters.getParameterRelation(editParameters.keys())+[2])<1:
                self.generateMesh(Laxis=np.array([runParameters['Laxis_X'],runParameters['Laxis_Y'],runParameters['Laxis_Z']]),endo_angle=runParameters['endo_angle'],epi_angle=runParameters['epi_angle'],fiberSheetletAngle=runParameters['fiberSheetletAngle'],fiberSheetletWidth=runParameters['fiberSheetletWidth'],radialFiberAngle=runParameters['radialFiberAngle'],fiberLength=runParameters['lr'],clipratio=runParameters['clip_ratio'],toRunCountFolder=True)
            elif not(os.path.isfile(self.casename+'/'+self.meshname+'.hdf5')):
                self.generateMesh(Laxis=np.array([runParameters['Laxis_X'],runParameters['Laxis_Y'],runParameters['Laxis_Z']]),endo_angle=runParameters['endo_angle'],epi_angle=runParameters['epi_angle'],fiberSheetletAngle=runParameters['fiberSheetletAngle'],fiberSheetletWidth=runParameters['fiberSheetletWidth'],radialFiberAngle=runParameters['radialFiberAngle'],fiberLength=runParameters['lr'],clipratio=runParameters['clip_ratio'],toRunCountFolder=False)
        if outputResultList is None:
            resultWriter=fenicsResultWriter(self.casename+"/"+str(self.runCount),['deformation','stress'])
        elif isinstance(outputResultList,str):
            resultWriter=fenicsResultWriter(self.casename+"/"+str(self.runCount),[outputResultList])
        else:
            resultWriter=fenicsResultWriter(self.casename+"/"+str(self.runCount),outputResultList)
        #displacementfile = fenics.File(self.casename+"/"+str(self.runCount)+"/deformation/u_disp.pvd")
        #stress_File = fenics.File(self.casename+"/"+str(self.runCount)+"/stress/_stress.pvd") #Joy Changed here

        
        ########################### Fenics's Newton  #########################################################

        # Closed loop cycle
        BCL = runParameters['BCL']
        if self.LVage=='adult':
            V_ven = runParameters['V_ven']
            V_art = runParameters['V_art']
            V_LA = runParameters['V_LA']
  
        cycle = 0.
        t = 0
        tstep = 0
        dt=1.
        t_array=[0]
        
        

        ####### Loading phase for LV ####################################################
        if self.manualVol is None:
            startVolume=runParameters['EDV_LV']
        elif runTimeList is None:
            startVolume=self.manualVol(0.)
        else:
            startVolume=self.manualVol(runTimeList[0])
        #if Cavityvol.vol>(1.1*runParameters['EDV_LV']):
        # 	raise Exception('Starting cavity volume,'+str(Cavityvol.vol)+', is more than target EDV_LV,'+str(runParameters['EDV_LV']))
        
        if setHeart is not None:
            if isinstance(setHeart,list):
                heart=lcleeHeart.heart(setHeart[0],setHeart[1],runParameters)
            else:
                heart=setHeart
            
        elif runParameters['EDP_LV'] is None:
            if min(self.defaultParameters.getParameterRelation(editParameters.keys())+[2])<1:
                heart=lcleeHeart.heart(self.casename+'/'+str(self.runCount),self.meshname,runParameters)
            else:
                heart=lcleeHeart.heart(self.casename,self.meshname,runParameters)
        else:
            if min(self.defaultParameters.getParameterRelation(editParameters.keys())+[2])<1:
        	    self.getUnloadedGeometry(editParameters=runParameters,casename=self.casename+'/'+str(self.runCount),targetPressure=runParameters['EDP_LV'],targetVolume=runParameters['EDV_LV'],savename=self.meshname+'_unloadedmesh',tryPressureRatio=tryPressureRatio)
        	    heart=lcleeHeart.heart(self.casename+'/'+str(self.runCount),self.meshname+'_unloadedmesh',runParameters)
            else:
        	    self.getUnloadedGeometry(editParameters=runParameters,savename=self.meshname+'_unloadedmesh',targetPressure=runParameters['EDP_LV'],targetVolume=runParameters['EDV_LV'],toRunCountFolder=True,tryPressureRatio=tryPressureRatio)
        	    heart=lcleeHeart.heart(self.casename,self.meshname+'_unloadedmesh',runParameters)
        heart.dt.dt = 1.
        if manualPhaseTimeInput and self.manualPhaseTime is not None:
            if runTimeList is None:
                heart.t_a.t_a=0.
            else:
                heart.t_a.t_a=self.manualPhaseTime(runTimeList[0])
        self.solveVolume(heart,startVolume)
        LVcav_array = [heart.uflforms.cavityvol()]
        Pcav_array = [heart.uflforms.cavitypressure()*0.0075]
        #Joy Changed from here
        #Stress calculation 
        
        ###stress_File << cauchy ## wei xuan remove
        if(fenics.MPI.rank(heart.comm) == 0): #added to output the volume, stroke volume data before the while loop
            fdataPV = open(self.casename+"/"+str(self.runCount)+"/PV_.txt", "w")
            fdatals = open(self.casename+"/"+str(self.runCount)+"/ls.txt", "w")
            fdataSV = open(self.casename+"/"+str(self.runCount)+"/SVp.txt", "w")
            fdataSVt = open(self.casename+"/"+str(self.runCount)+"/SVpt.txt", "w")
            fdata_stress = open( self.casename+"/"+str(self.runCount)+"/_stress_.txt", "w")
            fdataWork = open(self.casename+"/"+str(self.runCount)+"/work.txt", "w")
        p_cav = heart.uflforms.cavitypressure()
        V_cav = heart.uflforms.cavityvol()
        if(fenics.MPI.rank(heart.comm) == 0):
        	logger.info("Cycle number = "+repr(cycle)+ " cell time = "+repr(t)+ " tstep = "+repr(tstep)+" dt = "+repr(heart.dt.dt))
        	print(tstep, p_cav*0.0075 , V_cav, file=fdataPV)
        
        #if(fenics.MPI.rank(heart.comm) == 0):
        #    displacementfile << heart.w.sub(0)
        #Joy changed  from here
        #cauchy1 =  heart.uflforms.Cauchy1() + heart.activeforms.cauchy()
        #cauchy = fenics.project(cauchy1,fenics.TensorFunctionSpace(heart.mesh, "DG", 1), form_compiler_parameters={"representation":"uflacs"})
        #cauchy.rename("cauchy_stress","cauchy_stress")
        #stress_File << cauchy
        
        resultWriter(heart)
        cauchy1 =  heart.uflforms.Cauchy1() + heart.activeforms.cauchy()
        cauchy = fenics.project(cauchy1,fenics.TensorFunctionSpace(heart.mesh, "DG", 1), form_compiler_parameters={"representation":"uflacs"})
        sigma_fiber_LV = heart.activeforms.CalculateFiberStress(sigma = cauchy, e_fiber = heart.f0, Vol = heart.mesh_volume, Mesh = heart.mesh)
       
       	#Joy chnaged till here
        #save files
       
        # Closed-loop phase
        if self.LVage!='adult' and self.manualVol is None:
            ngspice_py.createLVcircuit(self.casename+'/'+self.meshname,runParameters)
        if runTimeList is not None:
            runTimeList=list(runTimeList)
            runTimeList.pop(0)
        if endTime is None:
            endTime = BCL*self.numofCycle
        while(cycle < self.numofCycle):
        	
        	#################################################################
        	p_cav = heart.uflforms.cavitypressure()
        	V_cav = heart.uflforms.cavityvol()
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
        	
        	if(t >= 0.0 and t < 4.0):
        		dt = 0.50
        	elif (t >= (runParameters['t0']-0.5) and t < (runParameters['t0']+0.5)):
        		dt = 3
        	else :
        		dt = 1.0
        
        	logger.info("t_a="+repr(t))
        
        	
        	if self.manualVol is not None:
        		V_cav=self.manualVol(tstep)
        	elif self.LVage=='adult':
        		V_cav,V_art,V_ven,V_LA=self.runWindkessel(heart.comm,tstep,dt,p_cav,V_cav,V_art,V_ven,V_LA)
        	else:
        		V_cav=self.runWindkessel(heart.comm,tstep,dt)

            
        	t_array.append(t)
        	LVcav_array.append(V_cav)
        	Pcav_array.append(p_cav*0.0075)
            
            #these lines added as between the two timesteps if there is a large volume it will crash, therefore below increases volume gradually
        	if manualPhaseTimeInput and self.manualPhaseTime is not None:
        		t2=self.manualPhaseTime(t)
        		heart.t_a.t_a=t2
        		if(t2 >= 0.0 and t2 < 4.0):
        			heart.dt.dt = 0.50
        		elif (t2 >= (runParameters['t0']-0.5) and t2 < (runParameters['t0']+0.5)):
        			heart.dt.dt = 3
        		else :
        			heart.dt.dt = 1.0
        	else:
        		heart.dt.dt=dt#runTimeList[0]-t
        		heart.t_a.t_a=t   
        	self.solveVolume(heart,V_cav)
        	#self.solveTa(heart,t,peaktime=runParameters['t0'])
            
        	
        	p_cav = heart.uflforms.cavitypressure()
        	V_cav = heart.uflforms.cavityvol()

        	if(fenics.MPI.rank(heart.comm) == 0):
        		logger.info("Cycle number = "+repr(cycle)+ " cell time = "+repr(t)+ " tstep = "+repr(tstep)+" dt = "+repr(heart.dt.dt))
        		print(tstep, p_cav*0.0075 , V_cav, file=fdataPV)
        	ls0 = runParameters['l0']
        	ls = fenics.sqrt(fenics.dot(heart.f0, heart.Cmat*heart.f0))*ls0
        	ls1 = fenics.project(ls,heart.Q).vector().get_local()[:]
        	eca = fenics.project(heart.activeforms.ECa(), heart.Q).vector().get_local()[:]
        	t_r = fenics.project(heart.activeforms.tr(), heart.Q).vector().get_local()[:]
        
        	if(fenics.MPI.rank(heart.comm) == 0):
        		print(heart.t_a.t_a, min(ls1), max(ls1), min(eca), max(eca), min(t_r), max(t_r), file=fdatals)

        	#if(fenics.MPI.rank(heart.comm) == 0):
        	#	displacementfile << heart.w.sub(0)
        	#Joy Chnaged from here
        	#cauchy = fenics.project(cauchy1,fenics.TensorFunctionSpace(heart.mesh, "DG", 1), form_compiler_parameters={"representation":"uflacs"})
        	#cauchy.rename("cauchy_stress","cauchy_stress")
        	#print (cauchy)
        	#print (cauchy.vector().array()[:])
        	#stress_File << cauchy
  
        	resultWriter(heart)
        	cauchy1 =  heart.uflforms.Cauchy1() + heart.activeforms.cauchy()
        	cauchy = fenics.project(cauchy1,fenics.TensorFunctionSpace(heart.mesh, "DG", 1), form_compiler_parameters={"representation":"uflacs"})
        	sigma_fiber_LV = heart.activeforms.CalculateFiberStress(sigma = cauchy, e_fiber = heart.f0, Vol = heart.mesh_volume, Mesh = heart.mesh)
        
        	if(fenics.MPI.rank(heart.comm) == 0):
        		print (fdata_stress, tstep, sigma_fiber_LV, file=fdata_stress )
        		print(tstep, heart.uflforms.strainEnergy(integrate=True) , file=fdataWork)
        	#Joy changed Till here          
        
        if(fenics.MPI.rank(heart.comm) == 0):
        	fdataPV.close()
        
        if (tstep >= runParameters['BCL']):
        	t_array=np.array(t_array)
        	temp_Ind=np.nonzero((t_array[2:]-t_array[:-2])>0)[0][-1]
        	temp_ind=np.nonzero(np.logical_and(t_array[:temp_Ind]>t_array[temp_Ind],t_array[:temp_Ind]<t_array[temp_Ind+2]))[0][-1]

        	LVcav_arrayT=LVcav_array[(temp_ind-1):]
        	Pcav_arrayT=Pcav_array[(temp_ind-1):]
        	plt.plot(LVcav_arrayT, Pcav_arrayT)
        	plt.xlabel("Volume")
        	plt.ylabel("Pressure")
        	plt.legend()
        	plt.savefig(self.casename+"/PV.png")
        	
        	print(max(LVcav_arrayT),min(LVcav_arrayT),max(LVcav_arrayT)-min(LVcav_arrayT), max(Pcav_arrayT),min(Pcav_arrayT),max(Pcav_arrayT)-min(Pcav_arrayT), file=fdataSV) 
        	print(LVcav_arrayT,Pcav_arrayT, file=fdataSVt)
        
        if(fenics.MPI.rank(heart.comm) == 0):
        	fdataSV.close()
        	fdataSVt.close()
        	fdata_stress.close() #Joy Changed here	
        	fdataWork.close()
        ######################################################################################################
        return runParameters
class fenicsResultWriter:
    def __new__(cls, *args, **kwargs):
        return super().__new__(cls)
    def __init__(self,savePath,outputResultList):
        self.savePath=savePath
        self.displacementFile=None
        self.stressFile=None
        self.strainFile=None
        self.workFile=None
        self.fiberFile=None
        if "displacement" in outputResultList:
            os.makedirs(self.savePath+"/deformation",exist_ok=True)
            self.displacementFile = fenics.File(self.savePath+"/deformation/u_disp.pvd")
        if "stress" in outputResultList:
            os.makedirs(self.savePath+"/stress",exist_ok=True)
            self.stressFile = fenics.File(self.savePath+"/stress/_stress.pvd") #Joy Changed here
        if "strain" in outputResultList:
            os.makedirs(self.savePath+"/strain",exist_ok=True)
            self.strainFile = []
            for n in ['ff','ss','nn','fs','fn','ns']:
                self.strainFile.append(fenics.File(self.savePath+"/strain/_strain_"+n+".pvd"))
        if "strain_energy_density" in outputResultList:
            os.makedirs(self.savePath+"/strain_energy_density",exist_ok=True)
            self.workFile = fenics.File(self.savePath+"/strain_energy_density/_strain_energy_density.pvd")
        if "fiber" in outputResultList:
            os.makedirs(self.savePath+"/fiber",exist_ok=True)
            self.fiberFile=[]
            for n in ['ff','ss','nn']:
                self.fiberFile.append(fenics.File(self.savePath+"/fiber/_fiber_"+n+".pvd"))
    def __call__(self,heart):
        if(fenics.MPI.rank(heart.comm) == 0):
            if self.displacementFile is not None:
                self.displacementFile << heart.w.sub(0)
            if self.stressFile is not None:
                cauchy1 =  heart.uflforms.Cauchy1() + heart.activeforms.cauchy()
                cauchy = fenics.project(cauchy1,fenics.TensorFunctionSpace(heart.mesh, "DG", 1), form_compiler_parameters={"representation":"uflacs"})
                cauchy.rename("cauchy_stress","cauchy_stress")
                self.stressFile << cauchy
            if self.strainFile is not None:
                Ea=heart.uflforms.Emat()
                f0 = heart.uflforms.parameters["fiber"]
                s0 = heart.uflforms.parameters["sheet"]
                n0 = heart.uflforms.parameters["sheet-normal"]
                Eff = fenics.project(fenics.inner(f0, Ea*f0),fenics.FunctionSpace(heart.mesh, "CG", 1), form_compiler_parameters={"representation":"uflacs"})
                Eff.rename("Eff_strain","Eff_strain")
                self.strainFile[0] << Eff
                Ess = fenics.project(fenics.inner(s0, Ea*s0),fenics.FunctionSpace(heart.mesh, "CG", 1), form_compiler_parameters={"representation":"uflacs"})
                Ess.rename("Ess_strain","Ess_strain")
                self.strainFile[1] << Ess
                Enn = fenics.project(fenics.inner(n0, Ea*n0),fenics.FunctionSpace(heart.mesh, "CG", 1), form_compiler_parameters={"representation":"uflacs"})
                Enn.rename("Enn_strain","Enn_strain")
                self.strainFile[2] << Enn
                Efs = fenics.project(fenics.inner(f0, Ea*s0),fenics.FunctionSpace(heart.mesh, "CG", 1), form_compiler_parameters={"representation":"uflacs"})
                Efs.rename("Efs_strain","Efs_strain")
                self.strainFile[3] << Efs
                Efn = fenics.project(fenics.inner(f0, Ea*n0),fenics.FunctionSpace(heart.mesh, "CG", 1), form_compiler_parameters={"representation":"uflacs"})
                Efn.rename("Efn_strain","Efn_strain")
                self.strainFile[4] << Efn
                Ens = fenics.project(fenics.inner(n0, Ea*s0),fenics.FunctionSpace(heart.mesh, "CG", 1), form_compiler_parameters={"representation":"uflacs"})
                Ens.rename("Ens_strain","Ens_strain")
                self.strainFile[5] << Ens
            if self.workFile is not None:
                work1=heart.uflforms.strainEnergy()
                work=fenics.project(work1,fenics.FunctionSpace(heart.mesh, "CG", 1), form_compiler_parameters={"representation":"uflacs"})
                work.rename("strain_energy_density","strain_energy_density")
                self.workFile << work
            if self.fiberFile is not None:
                Fmat1=heart.uflforms.Fmat()
                f0 = heart.uflforms.parameters["fiber"]
                s0 = heart.uflforms.parameters["sheet"]
                n0 = heart.uflforms.parameters["sheet-normal"]
                fiber_ff = fenics.project(Fmat1*f0,fenics.VectorFunctionSpace(heart.mesh, "CG", 1), form_compiler_parameters={"representation":"uflacs"})
                fiber_ff.rename("fiber_ff","fiber_ff")
                self.fiberFile[0] << fiber_ff
                fiber_ss = fenics.project(Fmat1*s0,fenics.VectorFunctionSpace(heart.mesh, "CG", 1), form_compiler_parameters={"representation":"uflacs"})
                fiber_ss.rename("fiber_ss","fiber_ss")
                self.fiberFile[1] << fiber_ss
                fiber_nn = fenics.project(Fmat1*n0,fenics.VectorFunctionSpace(heart.mesh, "CG", 1), form_compiler_parameters={"representation":"uflacs"})
                fiber_nn.rename("fiber_nn","fiber_nn")
                self.fiberFile[2] << fiber_nn
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
def plotLVBehavior(mainPath,labelList=None,timeSlice=None,volumeSlice=None,fmt='png',show=False):
    if not(isinstance(mainPath,list)):
        mainPath=[mainPath]
    removeVolumeratio=0.1
    VolTime_total=[]
    timeSpace=[]
    volumeSpace=[]
    for n in range(len(mainPath)):
        timeSpace.append(np.loadtxt(mainPath[n]+'_Press_timeSpace.txt'))
        volumeSpace.append(np.loadtxt(mainPath[n]+'_Press_volumeSpace.txt'))
        VolTime_base=np.loadtxt(mainPath[n]+'_Press_VolTime_base.txt')
        VolTime=np.loadtxt(mainPath[n]+'_Press_VolTime.txt')
        
        remove_volumeSpace=int(len(volumeSpace)*removeVolumeratio)
        remove_timeSpace=int(len(timeSpace)*0.5)
        volumeSpace[-1]=volumeSpace[-1][remove_volumeSpace:-remove_volumeSpace]
        VolTime_base=VolTime_base[remove_volumeSpace:-remove_volumeSpace]
        VolTime=VolTime[remove_volumeSpace:-remove_volumeSpace]
        timeSpace[-1]=timeSpace[-1][:-remove_timeSpace]
        VolTime=VolTime[:,:-remove_timeSpace]
        
        VolTime_total.append(VolTime_base.reshape((-1,1))+VolTime)
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
                label=labelList[m]
            plt.plot(volumeSpace[m],data,color=cList[m],label=label)
        ax.set_title('LV Behavior at {0:.1f}ms'.format(timeSlice))
        ax.set_xlabel('Volume (mL)')
        ax.set_ylabel('Pressure (mmHg)')
        if len(mainPath)>1:
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
                label=labelList[m]
            plt.plot(volumeSpace[m][::-1],data,color=cList[m],label=label)
        ax.set_title('LV Behavior at {0:.1f}mL'.format(volumeSlice))
        ax.set_xlabel('Time (ms)')
        ax.set_ylabel('Pressure (mmHg)')
        if len(mainPath)>1:
            ax.legend(loc=0)
        plt.savefig(mainPath[0]+'_LVBehaviorPlot_volume_{0:.1f}mL.png'.format(volumeSlice),dpi=fig.dpi,format=fmt,bbox_inches='tight')
        fig.clf()
    
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    X, Y = np.meshgrid(timeSpace, volumeSpace)
    
    # Customize the z axis.
    vmax=VolTime_total.max()
    vmin=VolTime_total.min()
    surf = ax.plot_surface(X, Y, VolTime_total,cmap=cm.coolwarm,linewidth=0, antialiased=False)

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
            return returnFunc(fdatacl[:,0]*1000., fdatacl[:,name.index(varName)])
        elif self.runMode=='iterativeRun':
            try:
                fdatacl = np.loadtxt(folderName+"/circuit_results.txt",skiprows=1)
                name=[' ', 'p_cav',' ', 'V_cav',' ', 'PLA',' ','PRV',' ','PRA', ' ','Ppa1',' ','Paa',' ','Plvvalv' ]
                return returnFunc(fdatacl[:,0]*1000., fdatacl[:,name.index(varName)])
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
            valueArray=valueArray[valueArray[:,0]>=(valueArray[:,0].max()-runParameters['BCL'])]
            cycletime=valueArray[:,0]%runParameters['BCL']
            diastolictimeInd=np.argmin(np.minimum(cycletime,runParameters['BCL']-cycletime))
            minval=np.min(valueArray[:,1])
            maxval=valueArray[diastolictimeInd,1]#np.max(valueArray[:,1])
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
            valueArray=valueArray[valueArray_p_cav[:,0]>=(valueArray_p_cav[:,0].max()-runParameters['BCL'])]
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
            valueArray=valueArray[valueArray_p_cav[:,0]>=(valueArray_p_cav[:,0].max()-runParameters['BCL'])]
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
            valueArray=valueArray[valueArray_p_cav[:,0]>=(valueArray_p_cav[:,0].max()-runParameters['BCL'])]
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
            valueArray=valueArray[valueArray_p_cav[:,0]>=(valueArray_p_cav[:,0].max()-runParameters['BCL'])]
            return func(valueArray[:,1].max(),ref)
        elif varName =='t_isorelax':
            try:
                ref=np.loadtxt(self.matchingFolder+"/"+varName+".txt")[0]
            except:
                ref=None
            valueArray=self.getspline(folderName,'V_cav')
            valueArray=valueArray[valueArray[:,0]>=(valueArray[:,0].max()-runParameters['BCL'])]
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
            valueArray=valueArray[valueArray[:,0]>=(valueArray[:,0].max()-runParameters['BCL'])]
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
    LVclosed_object.generateMesh(Laxis=np.array([LVclosed_object.defaultParameters['Laxis_X'],LVclosed_object.defaultParameters['Laxis_Y'],LVclosed_object.defaultParameters['Laxis_Z']]),endo_angle=LVclosed_object.defaultParameters['endo_angle'],epi_angle=LVclosed_object.defaultParameters['epi_angle'],fiberSheetletAngle=LVclosed_object.defaultParameters['fiberSheetletAngle'],fiberSheetletWidth=LVclosed_object.defaultParameters['fiberSheetletWidth'],radialFiberAngle=LVclosed_object.defaultParameters['radialFiberAngle'],fiberLength=LVclosed_object.defaultParameters['lr'],clipratio=LVclosed_object.defaultParameters['clip_ratio'],toRunCountFolder=True)
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
        heart=lcleeHeart.heart(LVclosed_object.casename+'/'+LVclosed_object.runCount,LVclosed_object.meshname+'_unloadedmesh',LVclosed_object.defaultParameters)
        LVclosed_object.solveVolume(heart,LVclosed_object.defaultParameters['ESV_LV'])
        temp_ESP=heart.uflforms.cavitypressure()*0.0075
        LVclosed_object.solveVolume(heart,LVclosed_object.defaultParameters['EDV_LV'])
        temp_EDP=heart.uflforms.cavitypressure()*0.0075
        end_start_Diastolic_pressure=temp_EDP-temp_ESP
        LVclosed_object.defaultParameters['EDP_LV']=temp_EDP
        LVclosed_object.writeRunParameters(LVclosed_object.casename+'/'+LVclosed_object.runCount,LVclosed_object.defaultParameters)
        if end_start_Diastolic_pressure*tune>maxLALV_pressure*tune:
            tune*=-1.
            stlPressure_adj*=0.3
    
    cmd = "cp -r " + LVclosed_object.casename+'/'+folder+'/'+LVclosed_object.meshname+'* '+ LVclosed_object.casename+'/'+folder+'/best_fit'
    os.system(cmd)
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
    np.savetxt(LVclosed_object.casename+'/'+folder+'/ref/p_cav.txt',np.array([[0.,LVclosed_object.defaultParameters['EDP_LV']],[LVclosed_object.defaultParameters['t0'],500.]]))
    cost=cost_function(LVclosed_object.defaultRunMode,matchingFolder=LVclosed_object.casename+'/'+folder+'/ref')
    cost.addCost('LeftVentriclePressureMeanSquare')
    linker=optimiser_linker(LVclosed_object,cost)
    linker.addVariable(optVarList)
    linker.set_kwargs(endTime=LVclosed_object.defaultParameters['t0']+3,outputResultList=[],tryPressureRatio=tryPressureRatio)
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
    np.savetxt(LVclosed_object.casename+'/'+folder+'/ref/p_cav.txt',np.array([[0.,LVclosed_object.defaultParameters['EDP_LV']],[time_in_ms,pressure_in_mmHg]]))
    cost=cost_function(LVclosed_object.defaultRunMode,matchingFolder=LVclosed_object.casename+'/'+folder+'/ref')
    cost.addCost('LeftVentriclePressureMeanSquare')
    linker=optimiser_linker(LVclosed_object,cost)
    linker.addVariable(['T0_LV'])
    linker.set_kwargs(endTime=time_in_ms+3,outputResultList=[],setHeart=setHeart)
    iniPara=np.array([LVclosed_object.defaultParameters['T0_LV']])
    optimize.minimize(linker,iniPara,bounds=[(1e4,LVclosed_object.defaultParameters['T0_LV']*1.1)],tol=0.001)

def optimiseWinkesselParameters(LVclosed_object,strokeVolume,maxLVLApressureDiff,maxLVAApressureDiff,step=0,stopstep=3,Q_regurge_flowratio=0.1,optVarList=None,costList=None,baseFolder=None,folder=None,folderToLVbehavior=None,runCycles=20):
    if folder is None:
        folder='optimiseWinkesselParameters'
    if baseFolder is not None:
        folder=baseFolder+'/'+folder
    if optVarList is None:
        optVarList=['lvlvvalvk','lvregurgevalveratio','Windkessel_scale_T0_LV']
        logs=[3.,3.,0]
    else:
        optVarList=['lvlvvalvk','lvregurgevalveratio','Windkessel_scale_T0_LV']+optVarList
        logs=[3.,3.,0]+list([0]*len(optVarList))
    if costList is None:
        costList=['StrokeVolumeMeanSquare','MaxLAVALVPressureMeanSquare','MaxLVVALVPressureMeanSquare']
    os.makedirs(LVclosed_object.casename+'/'+folder+'/ref',exist_ok=True)
    weights=[]
    for n in range(len(costList)):
        if 'StrokeVolume' in costList[n]:
            weights.append(10./LVclosed_object.defaultParameters['ESV_LV'])
            np.savetxt(LVclosed_object.casename+'/'+folder+'/ref/V_sv.txt',np.array([strokeVolume]))
        if 'MaxLAVALVPressure' in costList[n]:
            weights.append(1./maxLVLApressureDiff)
            np.savetxt(LVclosed_object.casename+'/'+folder+'/ref/pmax_LAVALVregurge.txt',np.array([maxLVLApressureDiff]))
        if 'MaxLVVALVPressure' in costList[n]:
            weights.append(2./(maxLVLApressureDiff+maxLVAApressureDiff))
            np.savetxt(LVclosed_object.casename+'/'+folder+'/ref/pmax_LVVALV.txt',np.array([maxLVAApressureDiff]))
    LVclosed_object.defaultRunMode='flowrateWindkesselRun'
    LVLA_Q=LVclosed_object.manualVol(np.arange(LVclosed_object.defaultParameters['BCL']*10.)/10.)
    LVLA_Q=(LVLA_Q[1:]-LVLA_Q[:-1])/0.1
    maxLVLA_Q=-LVLA_Q.min()*1000.
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
        LVclosed_object.defaultParameters['lvregurgevalveratio']=results['lvregurgevalveratio']
        if 'lvlvvalvr' in optVarList:
            LVclosed_object.defaultParameters['lvlvvalvr']=results['lvlvvalvr']
        os.makedirs(LVclosed_object.casename+'/'+folder+'/init_lvlvvalvk/ref',exist_ok=True)
        np.savetxt(LVclosed_object.casename+'/'+folder+'/init_lvlvvalvk/ref/pmax_LAVALVregurge.txt',np.array([maxLVLApressureDiff]))
        np.savetxt(LVclosed_object.casename+'/'+folder+'/init_lvlvvalvk/ref/pmax_LVVALV.txt',np.array([maxLVAApressureDiff]))
        
        for n in [1]:
            LVclosed_object.runCount=folder+'/init_lvlvvalvk/0'
            cost=cost_function(LVclosed_object.defaultRunMode,matchingFolder=LVclosed_object.casename+'/'+folder+'/init_lvlvvalvk/ref')
            cost.addCost(['MaxLAVALVPressureMeanSquare','MaxLVVALVPressureMeanSquare'],weights=[1./maxLVLApressureDiff,2./(maxLVLApressureDiff+maxLVAApressureDiff)])
            linker=optimiser_linker(LVclosed_object,cost)
            if 'lvlvvalvr' in optVarList:
                linker.addVariable(['lvlvvalvk','lvregurgevalveratio','lvlvvalvr'],[(LVclosed_object.defaultParameters['lvlvvalvk']/3.**n,LVclosed_object.defaultParameters['lvlvvalvk']*3.**n),(LVclosed_object.defaultParameters['lvregurgevalveratio']/10.**n,LVclosed_object.defaultParameters['lvregurgevalveratio']*10.**n),(LVclosed_object.defaultParameters['lvlvvalvr']/2.**n,LVclosed_object.defaultParameters['lvlvvalvr']*3.**n)],log=[3.,10.,3.])
                optimize.minimize(linker,linker.invOperate_para(np.array([LVclosed_object.defaultParameters['lvlvvalvk'],LVclosed_object.defaultParameters['lvregurgevalveratio'],LVclosed_object.defaultParameters['lvlvvalvr']])),jac='3-point',options={'eps':0.00001,'maxiter':8},tol=0.0001)
            else:
                linker.addVariable(['lvlvvalvk','lvregurgevalveratio'],[(LVclosed_object.defaultParameters['lvlvvalvk']/3.**n,LVclosed_object.defaultParameters['lvlvvalvk']*3.**n),(LVclosed_object.defaultParameters['lvregurgevalveratio']/10.**n,LVclosed_object.defaultParameters['lvregurgevalveratio']*10.**n)],log=[3.,10.])
                optimize.minimize(linker,linker.invOperate_para(np.array([LVclosed_object.defaultParameters['lvlvvalvk'],LVclosed_object.defaultParameters['lvregurgevalveratio']])),jac='3-point',options={'eps':0.00001,'maxiter':5},tol=0.0001)
            results=LVclosed_object.readRunParameters(LVclosed_object.casename+'/'+folder+'/init_lvlvvalvk/best_fit')
            LVclosed_object.defaultParameters['lvlvvalvk']=results['lvlvvalvk']
            LVclosed_object.defaultParameters['lvregurgevalveratio']=results['lvregurgevalveratio']
            if 'lvlvvalvr' in optVarList:
                LVclosed_object.defaultParameters['lvlvvalvr']=results['lvlvvalvr']
    if os.path.isfile(LVclosed_object.casename+'/'+folder+'/init_lvlvvalvk/best_fit/runParameters.txt'):
        results=LVclosed_object.readRunParameters(LVclosed_object.casename+'/'+folder+'/init_lvlvvalvk/best_fit')
    else:
        results=LVclosed_object.readRunParameters(LVclosed_object.casename+'/'+folder+'/init_Qratio/best_fit')
    LVclosed_object.defaultParameters['lvlvvalvk']=results['lvlvvalvk']
    LVclosed_object.defaultParameters['lvregurgevalveratio']=results['lvregurgevalveratio']
    if 'lvlvvalvr' in optVarList:
        LVclosed_object.defaultParameters['lvlvvalvr']=results['lvlvvalvr']
    LVclosed_object.defaultRunMode='fullWindkesselRun'
    if step<=2 and stopstep>=2:
    
        #get appropriate Windkessel_scale_T0_LV
        if step==2 and os.path.isfile(LVclosed_object.casename+'/'+folder+'/init_WindScale/best_fit/runParameters.txt'):
            init_Windkessel_scale_T0_LV=LVclosed_object.readRunParameters(LVclosed_object.casename+'/'+folder+'/init_WindScale/best_fit')['Windkessel_scale_T0_LV']
        else:
            LVclosed_object.runCount=folder+'/init_WindScale/0'
            LVclosed_object.fullWindkesselRun(unloadGeo=True,folderToLVbehavior=folderToLVbehavior,runCycles=runCycles)
            results_LVBehavior=getCircuitResult(LVclosed_object.casename+'/'+LVclosed_object.runCount,'v(lvvol)',last_cycle=True)
            temp_SV=results_LVBehavior.max()-results_LVBehavior.min()
            init_Windkessel_scale_T0_LV=min(0.9999999,max(0.0000001,LVclosed_object.defaultParameters['Windkessel_scale_T0_LV']*strokeVolume/temp_SV))
        #initialise 'Windkessel_scale_T0_LV' with Stroke Volume
        os.makedirs(LVclosed_object.casename+'/'+folder+'/init_WindScale/ref',exist_ok=True)
        np.savetxt(LVclosed_object.casename+'/'+folder+'/init_WindScale/ref/V_sv.txt',np.array([strokeVolume]))
        LVclosed_object.runCount=folder+'/init_WindScale/0'
        cost=cost_function(LVclosed_object.defaultRunMode,matchingFolder=LVclosed_object.casename+'/'+folder+'/init_WindScale/ref')
        cost.addCost(['StrokeVolumeMeanSquare'])
        linker=optimiser_linker(LVclosed_object,cost)
        linker.addVariable(['Windkessel_scale_T0_LV'],[(init_Windkessel_scale_T0_LV*0.5,min(1,init_Windkessel_scale_T0_LV*1.5))])
        linker.set_kwargs(folderToLVbehavior=folderToLVbehavior,runCycles=runCycles)
        bounds=[(float('-inf'),float('inf'))]
        optimize.minimize(linker,linker.invOperate_para(np.array([init_Windkessel_scale_T0_LV])),jac='3-point',options={'eps':0.000001,'maxiter':1},tol=0.0001)
        results=LVclosed_object.readRunParameters(LVclosed_object.casename+'/'+folder+'/init_WindScale/best_fit')['Windkessel_scale_T0_LV']
        LVclosed_object.defaultParameters['Windkessel_scale_T0_LV']=results
        for n in range(1):
            if abs(results-strokeVolume)/strokeVolume<0.01:
                break
            LVclosed_object.runCount=folder+'/init_WindScale/0'
            cost=cost_function(LVclosed_object.defaultRunMode,matchingFolder=LVclosed_object.casename+'/'+folder+'/init_WindScale/ref')
            cost.addCost(['StrokeVolumeMeanSquare'])
            linker=optimiser_linker(LVclosed_object,cost)
            linker.addVariable(['Windkessel_scale_T0_LV'],[(init_Windkessel_scale_T0_LV*0.5,min(1,init_Windkessel_scale_T0_LV*1.5))])
            linker.set_kwargs(folderToLVbehavior=folderToLVbehavior,runCycles=runCycles)
            bounds=[(float('-inf'),float('inf'))]
            optimize.minimize(linker,linker.invOperate_para(np.array([init_Windkessel_scale_T0_LV])),jac='3-point',options={'eps':0.000001,'maxiter':2},tol=0.0001)
            results=LVclosed_object.readRunParameters(LVclosed_object.casename+'/'+folder+'/init_WindScale/best_fit')['Windkessel_scale_T0_LV']
            LVclosed_object.defaultParameters['Windkessel_scale_T0_LV']=results
    results=LVclosed_object.readRunParameters(LVclosed_object.casename+'/'+folder+'/init_WindScale/best_fit')['Windkessel_scale_T0_LV']
        
    LVclosed_object.defaultParameters['Windkessel_scale_T0_LV']=results
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
                if optVarList[n]=='Windkessel_scale_T0_LV':
                    all_limits.append((LVclosed_object.defaultParameters['Windkessel_scale_T0_LV']*0.5,min(1.,LVclosed_object.defaultParameters['Windkessel_scale_T0_LV']*1.5)))
                    iniPara.append(LVclosed_object.defaultParameters['Windkessel_scale_T0_LV'])
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
            

def optimiseAllWinkesselParameters(LVclosed_object,strokeVolume,maxLVLApressureDiff,maxLVLVVALVpressureDiff,Q_regurge_flowratio=0.1,folderToLVbehavior=None,runCycles=20):
    baseFolder='optimiseWinkesselParameters'
    folder='0'
    EDVol=LVclosed_object.defaultParameters['EDV_LV']
    ESVol=LVclosed_object.defaultParameters['ESV_LV']
    os.makedirs(LVclosed_object.casename+'/'+baseFolder+'/'+folder,exist_ok=True)
    if not(os.path.isfile(LVclosed_object.casename+'/'+baseFolder+'/'+folder+'/init_WindScale/best_fit/circuit_results_lastcycle.txt')):
        if not(os.path.isfile(LVclosed_object.casename+'/'+baseFolder+'/'+folder+'/init_lvlvvalvk/best_fit/circuit_results_lastcycle.txt')):
            optimiseWinkesselParameters(LVclosed_object,strokeVolume,maxLVLApressureDiff,maxLVLVVALVpressureDiff,step=0,stopstep=2,Q_regurge_flowratio=Q_regurge_flowratio,baseFolder=baseFolder,folder=folder,folderToLVbehavior=folderToLVbehavior,runCycles=runCycles)
        else:
            optimiseWinkesselParameters(LVclosed_object,strokeVolume,maxLVLApressureDiff,maxLVLVVALVpressureDiff,step=2,stopstep=2,Q_regurge_flowratio=Q_regurge_flowratio,baseFolder=baseFolder,folder=folder,folderToLVbehavior=folderToLVbehavior,runCycles=runCycles)
    for startiter in range(1,10):
        if not(os.path.isfile(LVclosed_object.casename+'/'+baseFolder+'/'+str(startiter)+'/init_WindScale/best_fit/circuit_results_lastcycle.txt')):
            break
    for n in range(startiter,max(10,startiter+2)):
        folder=str(n)
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
        while time[-1]>LVclosed_object.defaultParameters['BCL']:
            time=time[:-1]
            Vol_LV=Vol_LV[:-1]
        if time[0]>0:
            time=np.concatenate(([0.],time))
            Vol_LV=np.concatenate(([Vol_LV[0]],Vol_LV))
        argmin=np.argmin(Vol_LV)
        Vol_LV[:argmin]=(Vol_LV[:argmin]-Vol_LV[argmin])/(Vol_LV[0]-Vol_LV[argmin])*(EDVol-ESVol)+ESVol
        Vol_LV[argmin:]=(Vol_LV[argmin:]-Vol_LV[argmin])/(Vol_LV[-1]-Vol_LV[argmin])*(EDVol-ESVol)+ESVol
        if time[-1]<LVclosed_object.defaultParameters['BCL']:
            time=np.concatenate((time,[LVclosed_object.defaultParameters['BCL']]))
            Vol_LV=np.concatenate((Vol_LV,[EDVol]))
        #Vol_LV=(Vol_LV-Vol_LV.min())/(Vol_LV[0]-Vol_LV.min())*(EDVol-ESVol)+ESVol
        np.savetxt(LVclosed_object.casename+'/'+baseFolder+'/'+str(n)+'/volume.txt',np.array([time,Vol_LV]).T)
        Q_regurge_flowratio=real_Qratio*0.7+Q_regurge_flowratio*0.3
        LVclosed_object.setManualVolumeFile(LVclosed_object.casename+'/'+baseFolder+'/'+str(n)+'/volume.txt')
        optimiseWinkesselParameters(LVclosed_object,strokeVolume,maxLVLApressureDiff,maxLVLVVALVpressureDiff,step=0,stopstep=2,Q_regurge_flowratio=Q_regurge_flowratio,optVarList=None,baseFolder=baseFolder,folder=folder,folderToLVbehavior=folderToLVbehavior,runCycles=runCycles)
