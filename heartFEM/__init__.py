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
'''
_version='2.3.1'
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
from matplotlib import pylab as plt
from petsc4py import PETSc

import math
import csv

import vtk
from heartFEM import ngspice_py
from heartFEM import heartParameters
from mpi4py import MPI as pyMPI
try:
    import lcleeHeart
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
def pointsMeansqDiff_kabsch(from_a,to_b):
    R,u=kabsch(from_a,to_b)
    new_a=R.dot(from_a.T).T+u.reshape((1,-1))
    return np.sum((new_a-to_b)**2.,axis=1).mean()
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
        while time>timedata[-1]:
            time-=self.defaultParameters['BCL']
        if time<timedata[0]:
            raise Exception('time '+str(time)+' , is not inside file given:'+self.manualVolFile)
        spl=scipyinterpolate.splrep(timedata, voldata)
        return scipyinterpolate.splev(np.array([time]), spl)[0]
    def setManualVolumeFile(self,volfile,timecol=0,timescale=1.,volcol=1,volscale=1.):
        self.manualVolFile=volfile
        self.manualVolFiletimecol=timecol
        self.manualVolFilevolcol=volcol
        self.manualVolFiletimescale=timescale
        self.manualVolFilevolscale=volscale
    def setDefaultWindkessel(self,modelString):
        self.defaultParameters.setDefaultWindkessel(modelString)
        if modelString[:5]=='fetal':
            try:
                self.setRVcurrentFourier(self.casename+'/'+self.meshname+'_rvflowrate.txt')
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
    def setRVcurrentFourier(self,filename,fourierTerms=1):
        fourierTerms=int(fourierTerms)
        data=np.loadtxt(filename)
        period=2.*data[-1,0]-data[-2,0]
        def fourierFit(x, *a):
            ret = 0.
            for deg in range(fourierTerms):
                ret += a[deg*2] * np.sin((deg+1) *2.* np.pi / period * x+a[deg*2+1])
            return ret
        popt, pcov = curve_fit(fourierFit, data[:,0], data[:,1],p0=np.ones(2*fourierTerms))

        self.defaultParameters['rvfunc']='fourier'+str(fourierTerms)
        self.defaultParameters['rvfuncarg']=[]
        for n in range(len(popt)):
            self.defaultParameters['rvfuncarg'].append('rvfuncarg'+str(n))
            popt_temp=popt[n]
            if n%2==1:
                while popt_temp<-np.pi:
                    popt_temp+=2.*np.pi
                while popt_temp>np.pi:
                    popt_temp-=2.*np.pi
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
    def generateMesh(self,Laxis=None,endo_angle='',epi_angle='',clipratio=0.95,meshsize=0.6,saveaddstr='',toRunCountFolder=False):
        if Laxis is None:
            Laxis=np.array([np.array([self.defaultParameters['Laxis_X'],self.defaultParameters['Laxis_Y'],self.defaultParameters['Laxis_Z']])])
        if toRunCountFolder:
            if isinstance(toRunCountFolder,str):
                if toRunCountFolder[0]!='/':
                    toRunCountFolder='/'+toRunCountFolder
            else:
                toRunCountFolder='/'+str(self.runCount)
            return lcleeHeart.generateMesh(self.casename,self.meshname,Laxis,endo_angle=endo_angle,epi_angle=epi_angle,clipratio=clipratio,meshsize=meshsize,meshname=self.meshname+saveaddstr,saveSubFolder=toRunCountFolder)
        else:
            return lcleeHeart.generateMesh(self.casename,self.meshname,Laxis,endo_angle=endo_angle,epi_angle=epi_angle,clipratio=clipratio,meshsize=meshsize,meshname=self.meshname+saveaddstr)
    def runWinkessel(self,comm,tstep,dt_dt,*args):
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
        	heart.dt.dt=dt
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
        maxtimeind=0
        for curr_volN in range(len(volumeSpace)):
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
            if len(press_volTime_temp)>maxtimeind:
                maxtimeind=len(press_volTime_temp)
            press_volTime_base[curr_volN]=min(press_volTime_temp[-1],press_volTime_temp[0])
            press_volTime[curr_volN,:len(press_volTime_temp)]=np.array(press_volTime_temp)-press_volTime_base[curr_volN]
            if len(timeSpace)==1:
                volPress[curr_volN,1]=press_volTime_base[curr_volN]+press_volTime[curr_volN,0]
            if(fenics.MPI.rank(heart.comm) == 0):
                if len(timeSpace)>1:
                    np.savetxt(self.casename+toRunCountFolder+"/"+meshname+"_Press_VolTime"+saveaddstr+".txt",press_volTime,header=str(curr_volN)+'solved rowVolm: '+str(volumeSpace)+'\ncolTime: '+str(timeSpace))
                    np.savetxt(self.casename+toRunCountFolder+"/"+meshname+"_Press_VolTime_base"+saveaddstr+".txt",press_volTime_base)
                else:
                    np.savetxt(self.casename+toRunCountFolder+"/"+meshname+"_"+volPressStr+saveaddstr+".txt",volPress,header="Volume, Pressure")
        # press_volTime=press_volTime[:,:maxtimeind]
        #timeSpace=timeSpace[:maxtimeind]
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
    def adjustMeshToPressure(self,casename,meshname,targetPressure,savename=None,softAdjust=float('inf'),usePrevious=3,iterationNumber=float('inf'),initHeart=None,prevDeform=None,saveFolder=''):
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
        while abs(heart.uflforms.cavitypressure()*0.0075/targetPressure-1)>10**-4. or abs(heart.uflforms.cavityvol()/targetVolume-1)>10**-4.:##alternate convergence criteria##(np.mean((heart.mesh.coordinates()[:]-du.vector()[:].reshape((-1,3))-refHeart.mesh.coordinates()[:])**2.).max()**1.5/refHeart.uflforms.cavityvol())<10.**-6.:
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
            else:
                refHeart=lcleeHeart.heart(casename,meshname,runParameters)
                newDeform=refHeart.w.sub(0).copy(deepcopy=True)
                newDeform.vector()[:]*=0.
                newDeform.vector()[:]-=heart.w.sub(0).vector()[:]
                fenics.ALE.move(refHeart.mesh,newDeform)
                self.savehdf5(refHeart,casename,meshname,'tempadjmesh_'+meshname,saveFolder=saveFolder)
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
        if not(usePrevious) and (usePreviousiteration>0) and (abs(heart.uflforms.cavitypressure()*0.0075/targetPressure-1)>10**-4. or abs(heart.uflforms.cavityvol()/targetVolume-1)>10**-4.):
            newDeform=heart.w.sub(0).copy(deepcopy=True)
            newDeform.vector()[:]*=-1
            newDeform.vector()[:]+=prevDeformVectors
            fenics.ALE.move(heart.mesh,newDeform)
            self.savehdf5(heart,casename,meshname,'tempadjmesh_'+meshname,saveFolder=saveFolder)
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
    def getUnloadedGeometry(self,editParameters=None,casename=None,meshname=None,targetPressure=None,targetVolume=None,tryPressureRatio=0.5,targetMeshCoords=None,savename=None,iterationNumber=0,toRunCountFolder=False):
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

        heart=self.adjustMeshToPressure(casename,meshname,tryPressure,savename=savename,usePrevious=0,iterationNumber=iterationNumber,saveFolder=toRunCountFolder)
        self.savehdf5(heart,casename,meshname,'tempadjmesh_'+meshname,saveFolder=toRunCountFolder)
            
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
            heart=self.adjustMeshToPressure(casename,meshname,tryPressure,savename=savename,usePrevious=0,iterationNumber=iterationNumber,saveFolder=toRunCountFolder)
            self.savehdf5(heart,casename,meshname,'tempadjmesh_'+meshname,saveFolder=toRunCountFolder)
            heart=lcleeHeart.heart(casename+toRunCountFolder,'tempadjmesh_'+meshname,runParameters)
            self.solveVolume(heart,targetVolume,voladj=0.05)
            if targetMeshCoords is not None:
                logging.info('Trying Pressure '+repr(tryPressure)+' between '+repr(meshRMS_pressurebound))
                logging.info('Errors '+repr(meshRMS_error))
            elif (tunePressure*heart.uflforms.cavitypressure()*0.0075)>(tunePressure*targetPressure):
                tunePressure*=-1
                adjtryPressure*=0.33
            if adjtryPressure<10**-6:
                break
        self.savehdf5(heart,casename,meshname,savename,saveFolder=toRunCountFolder)
        return heart
    def __call__(self,editParameters=None,**kwargs):
        if self.defaultRun_kwargs is None:
            run_kwargs={}
        else:
            run_kwargs=self.defaultRun_kwargs
        for key in kwargs:
            run_kwargs[key]=kwargs[key]
        if self.defaultRunMode=='iterativeRun':
            runParameters=self.iterativeRun(editParameters=editParameters,**run_kwargs)
        if self.defaultRunMode=='LVbehaviorRun':
            runParameters=self.LVbehaviorRun(editParameters=editParameters,**run_kwargs)
        if self.defaultRunMode=='unloadedGeometry':
            runParameters=self.unloadedGeometryRun(editParameters=editParameters,**run_kwargs)
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
    def unloadedGeometryRun(self,editParameters=None):
        runParameters=self.defaultParameters.copy()
        if editParameters is not None:
            for key in runParameters:
                if key in editParameters:
                    runParameters[key]=editParameters[key]
        else:
            editParameters={}
        os.makedirs(self.casename+"/"+str(self.runCount),exist_ok=True)
        with open(self.casename+'/'+str(self.runCount)+"/runParameters.txt", "w") as f:
                for key, val in runParameters.items():
                    f.write("{0:s} , {1:s}\n".format(key,str(val)))
        #if self.runCount<=4:
        #    return runParameters
        print(self.runCount,"start generate mesh")
        self.generateMesh(Laxis=np.array([runParameters['Laxis_X'],runParameters['Laxis_Y'],runParameters['Laxis_Z']]),endo_angle=runParameters['endo_angle'],epi_angle=runParameters['epi_angle'],clipratio=runParameters['clip_ratio'],toRunCountFolder=True)
        print(self.runCount,"end generate mesh")
        heart=self.getUnloadedGeometry(editParameters=runParameters,casename=self.casename+"/"+str(self.runCount),savename=self.meshname+'_unloadedmesh',toRunCountFolder=False)
        self.solveVolume(heart,runParameters['EDV_LV'],voladj=0.05)
        fenics.ALE.move(heart.mesh,heart.w.sub(0))
        np.savetxt(self.casename+"/"+str(self.runCount)+"/coordsDiaMesh.txt",heart.mesh.coordinates()[:])
        return runParameters
    def LVbehaviorRun(self,editParameters=None,unloadGeo=False,folderToLVbehavior=None):
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
        os.makedirs(self.casename+"/"+str(self.runCount),exist_ok=True)
        with open(self.casename+'/'+str(self.runCount)+"/runParameters.txt", "w") as f:
            for key, val in runParameters.items():
                f.write("{0:s} , {1:s}\n".format(key,str(val)))
        if min(self.defaultParameters.getParameterRelation(editParameters.keys())+[2])<1:
            self.generateMesh(Laxis=np.array([runParameters['Laxis_X'],runParameters['Laxis_Y'],runParameters['Laxis_Z']]),endo_angle=runParameters['endo_angle'],epi_angle=runParameters['epi_angle'],clipratio=runParameters['clip_ratio'],toRunCountFolder=True)
        elif not(os.path.isfile(self.casename+'/'+self.meshname+'.hdf5')):
            self.generateMesh(Laxis=np.array([runParameters['Laxis_X'],runParameters['Laxis_Y'],runParameters['Laxis_Z']]),endo_angle=runParameters['endo_angle'],epi_angle=runParameters['epi_angle'],clipratio=runParameters['clip_ratio'],toRunCountFolder=False)
        if unloadGeo:
            meshname=self.meshname+'_unloadedmesh'
        else:
            meshname=self.meshname
        
        if unloadGeo:
            if folderToLVbehavior is not None and min(self.defaultParameters.getParameterRelation(editParameters.keys())+[2])>0:
                if not(os.path.isfile(self.casename+folderToLVbehavior+'/'+self.meshname+'_unloadedmesh.hdf5')):
                    self.getUnloadedGeometry(editParameters=runParameters,savename=self.meshname+'_unloadedmesh',toRunCountFolder=folderToLVbehavior)
            elif min(self.defaultParameters.getParameterRelation(editParameters.keys())+[2])<1:
                self.getUnloadedGeometry(editParameters=runParameters,casename=self.casename+'/'+str(self.runCount),savename=self.meshname+'_unloadedmesh')
            else:
                self.getUnloadedGeometry(editParameters=runParameters,savename=self.meshname+'_unloadedmesh',toRunCountFolder=True)
        volstep=runParameters['ESV_LV']*0.9*(runParameters['EDV_LV']/runParameters['ESV_LV']*1.1/0.9)**(np.linspace(0,1,num=50))[::-1]
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
            ngspice_py.generateLVtable(self.casename+"/"+str(self.runCount)+'/'+meshname,runParameters['BCL'],loading_casename=self.casename+'/'+meshname)
        else:
            ngspice_py.generateLVtable(self.casename+"/"+str(self.runCount)+'/'+meshname,runParameters['BCL'],loading_casename=self.casename+folderToLVbehavior+'/'+meshname)
        cmd = "cp "+self.casename+'/'+self.meshname+'_rvflowrate.txt'+" " + self.casename+"/"+str(self.runCount)+'/'+meshname+'_rvflowrate.txt'
        os.system(cmd)
        if 'ES_time' not in runParameters:
            ngspice_py.createLVcircuit(self.casename+"/"+str(self.runCount)+'/'+meshname,runParameters)
            ngspice_py.simLVcircuit(self.casename+"/"+str(self.runCount)+'/'+meshname,runParameters['BCL']*10,self.casename+"/"+str(self.runCount)+'/'+meshname+'_lvcirtable.txt',lvinputvar='table2d',initLVvol=runParameters['EDV_LV'])
        else:
            if runParameters['ES_time'] is None:
                ngspice_py.createLVcircuit(self.casename+"/"+str(self.runCount)+'/'+meshname,runParameters)
            else:
                ngspice_py.createLVcircuit(self.casename+"/"+str(self.runCount)+'/'+meshname,runParameters,skipVariableList=["timetopeaktension"])
            ngspice_py.simLVcircuit_align_EDvol_and_EStime(self.casename+"/"+str(self.runCount)+'/'+meshname,runParameters['BCL']*10,self.casename+"/"+str(self.runCount)+'/'+meshname+'_lvcirtable.txt',runParameters['BCL'],runParameters['ES_time'],runParameters['t0'],lvinputvar='table2d',initLVvol=runParameters['EDV_LV'])
        return runParameters
    def iterativeRun(self,editParameters=None,always_remesh=False,runTimeList=None):
        runParameters=self.defaultParameters.copy()
        if editParameters is not None:
            for key in runParameters:
                if key in editParameters:
                    runParameters[key]=editParameters[key]
        else:
            editParameters={}
        os.makedirs(self.casename+"/"+str(self.runCount),exist_ok=True)
        os.makedirs(self.casename+"/"+str(self.runCount)+"/stress",exist_ok=True)
        with open(self.casename+'/'+str(self.runCount)+"/runParameters.txt", "w") as f:
                for key, val in runParameters.items():
                    f.write("{0:s} , {1:s}\n".format(key,str(val)))
        pressureIteration_factor=0.1
        
        displacementfile = fenics.File(self.casename+"/"+str(self.runCount)+"/deformation/u_disp.pvd")
        activestressfile = fenics.File(self.casename+"/"+str(self.runCount)+"/activestress/stress.pvd") #CW
        stress_File = fenics.File(self.casename+"/"+str(self.runCount)+"/stress/_stress.pvd") #Joy Changed here

        while True: #setup
            if always_remesh or str(self.runCount)=='1' or 'Laxis_X' in editParameters or 'Laxis_Y' in editParameters or 'Laxis_Z' in editParameters or 'endo_angle' in editParameters or 'epi_angle' in editParameters:
                if runParameters['Laxis_X'] is None or runParameters['Laxis_Y'] is None or runParameters['Laxis_Z'] is None:
                    self.getLongAxis()
                    runParameters['Laxis_X']=self.defaultParameters['Laxis_X']
                    runParameters['Laxis_Y']=self.defaultParameters['Laxis_Y']
                    runParameters['Laxis_Z']=self.defaultParameters['Laxis_Z']
                self.generateMesh(Laxis=np.array([runParameters['Laxis_X'],runParameters['Laxis_Y'],runParameters['Laxis_Z']]),endo_angle=runParameters['endo_angle'] ,epi_angle=runParameters['epi_angle'])
            heart=lcleeHeart.heart(self.casename,self.meshname,runParameters)
            #solver,uflforms,activeforms,t_a,dt,Cavityvol,f0,mesh_volume,mesh,comm=lcleeHeart.setupHeart(self.casename,self.meshname,runParameters)
            
            ########################### Fenics's Newton  #########################################################
            heart.Cavityvol.vol = heart.uflforms.cavityvol()
            
            # Closed loop cycle
            BCL = runParameters['BCL']
            if self.LVage=='adult':
                V_ven = runParameters['V_ven']
                V_art = runParameters['V_art']
                V_LA = runParameters['V_LA']
  
            cycle = 0.
            t = 0
            tstep = 0
            heart.dt.dt = 1.
            dt=1.
            
            t_array=[0]
            LVcav_array = [heart.uflforms.cavityvol()]
            Pcav_array = [heart.uflforms.cavitypressure()*0.0075]
            
            #Joy Changed from here
            #Stress calculation 
            cauchy1 =  heart.uflforms.Cauchy1() + heart.activeforms.cauchy()
            cauchy = fenics.project(cauchy1,fenics.TensorFunctionSpace(heart.mesh, "DG", 1), form_compiler_parameters={"representation":"uflacs"})
            cauchy.rename("cauchy_stress","cauchy_stress")
            ###stress_File << cauchy ## wei xuan remove

            ####### Loading phase for LV ####################################################
            if self.manualVol is None:
                startVolume=runParameters['EDV_LV']
            elif runTimeList is None:
                startVolume=self.manualVol(0.)
            else:
                startVolume=self.manualVol(runTimeList[0])
            #if Cavityvol.vol>(1.1*runParameters['EDV_LV']):
            # 	raise Exception('Starting cavity volume,'+str(Cavityvol.vol)+', is more than target EDV_LV,'+str(runParameters['EDV_LV']))
            self.solveVolume(heart,startVolume)
            
            if runParameters['EDP_LV'] is None or self.manualVol is not None:
                break
            else:
                ###need to edit
            	self.getUnloadedGeometry(casename=self.casename,meshname=self.meshname,targetPressure=runParameters['EDP_LV'],targetVolume=runParameters['EDV_LV'],savename=self.meshname+'_unloadedmesh',iterationNumber=0)
            	heart=lcleeHeart.heart(self.casename,self.meshname+'_unloadedmesh',runParameters)
        if(fenics.MPI.rank(heart.comm) == 0): #added to output the volume, stroke volume data before the while loop
            fdataPV = open(self.casename+"/"+str(self.runCount)+"/PV_.txt", "w")
            fdatals = open(self.casename+"/"+str(self.runCount)+"/ls.txt", "w")
            fdataSV = open(self.casename+"/"+str(self.runCount)+"/SVp.txt", "w")
            fdataSVt = open(self.casename+"/"+str(self.runCount)+"/SVpt.txt", "w")
            fdata_stress = open( self.casename+"/"+str(self.runCount)+"/stress/_stress_.txt", "w")

        p_cav = heart.uflforms.cavitypressure()
        V_cav = heart.uflforms.cavityvol()
        if(fenics.MPI.rank(heart.comm) == 0):
        	logger.info("Cycle number = "+repr(cycle)+ " cell time = "+repr(t)+ " tstep = "+repr(tstep)+" dt = "+repr(heart.dt.dt))
        	print(tstep, p_cav*0.0075 , V_cav, file=fdataPV)


        if(fenics.MPI.rank(heart.comm) == 0):
            displacementfile << heart.w.sub(0)
        #Joy changed  from here
        cauchy = fenics.project(cauchy1,fenics.TensorFunctionSpace(heart.mesh, "DG", 1), form_compiler_parameters={"representation":"uflacs"})
        cauchy.rename("cauchy_stress","cauchy_stress")
        stress_File << cauchy

        sigma_fiber_LV = heart.activeforms.CalculateFiberStress(sigma = cauchy, e_fiber = heart.f0, Vol = heart.mesh_volume, Mesh = heart.mesh)
       
       	#Joy chnaged till here
        #save files
       
        # Closed-loop phase
        if self.LVage!='adult' and self.manualVol is None:
            ngspice_py.createLVcircuit(self.casename+'/'+self.meshname,runParameters)
        if runTimeList is not None:
            runTimeList=list(runTimeList)
            runTimeList.pop(0)
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
        		V_cav,V_art,V_ven,V_LA=self.runWinkessel(heart.comm,tstep,dt,p_cav,V_cav,V_art,V_ven,V_LA)
        	else:
        		V_cav=self.runWinkessel(heart.comm,tstep,dt)

            
        	t_array.append(t)
        	LVcav_array.append(V_cav)
        	Pcav_array.append(p_cav*0.0075)
            
            #these lines added as between the two timesteps if there is a large volume it will crash, therefore below increases volume gradually
           
        	self.solveVolume(heart,V_cav)
        	self.solveTa(heart,t,peaktime=runParameters['t0'])
            
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

        	if(fenics.MPI.rank(heart.comm) == 0):
        		displacementfile << heart.w.sub(0)
        	#Joy Chnaged from here
        	cauchy = fenics.project(cauchy1,fenics.TensorFunctionSpace(heart.mesh, "DG", 1), form_compiler_parameters={"representation":"uflacs"})
        	cauchy.rename("cauchy_stress","cauchy_stress")
        	#print (cauchy)
        	#print (cauchy.vector().array()[:])
        	stress_File << cauchy
  
        	sigma_fiber_LV = heart.activeforms.CalculateFiberStress(sigma = cauchy, e_fiber = heart.f0, Vol = heart.mesh_volume, Mesh = heart.mesh)
        
        	if(fenics.MPI.rank(heart.comm) == 0):
        		print (fdata_stress, tstep, sigma_fiber_LV, file=fdata_stress ) 
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
        ######################################################################################################
        return runParameters
        
class optimiser_linker:
    def __new__(cls, *args, **kwargs):
        return super().__new__(cls)
    def __init__(self,runSimulationClass,calculateCostClass):
        self.runSimulation=runSimulationClass
        self.calculateCost=calculateCostClass
        self.variableList=[]
        self.minCost=float('inf')
    def addVariable(self,variables):
        if isinstance(variables,str):
            variables=[variables]
        self.variableList+=variables
    def __call__(self,para):
        if len(self.variableList)==0:
            raise Exception('Select variable to optimise with optimiser.addVariable(VARIABLE_STRING[s])')
        os.makedirs(self.runSimulation.casename+'/best_fit',exist_ok=True)
        editParameters={}
        for n in range(len(self.variableList)):
            editParameters[self.variableList[n]]=para[n]
        try:
            self.runSimulation.runCount+=1
            runParameters=self.runSimulation(editParameters=editParameters)  
        except Exception as e_inst:
            with open(self.runSimulation.casename+'/'+str(self.runSimulation.runCount)+"/runParameters.txt", "r+") as f:
                 old = f.read() # read everything in the file
                 f.seek(0) # rewind
                 f.write("Status , FAILED\n cost , "+str(type(e_inst))+'  '+str(type(e_inst))+"\n" + old) # write the new line before
            return float('inf')
        else:
            cost_current=self.calculateCost(self.runSimulation.casename+'/'+str(self.runSimulation.runCount),runParameters)
            with open(self.runSimulation.casename+'/'+str(self.runSimulation.runCount)+"/runParameters.txt", "r+") as f:
                 old = f.read() # read everything in the file
                 f.seek(0) # rewind
                 f.write("Status , SUCCESS\n cost , "+str(cost_current)+"\n" + old) # write the new line before
            if cost_current<self.minCost:
                self.minCost=cost_current
                cmd = "rm -r " + self.runSimulation.casename+'/best_fit'
                os.system(cmd)
                os.makedirs(self.runSimulation.casename+'/best_fit',exist_ok=True)
                cmd = "cp -r " + self.runSimulation.casename+'/'+str(self.runSimulation.runCount)+'/* '+ self.runSimulation.casename+'/best_fit'
                os.system(cmd)
            return cost_current
        
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
        self.collateNameSwitch={'MeanSquare':self.MeanSquare
                            }
        self.varNameSwitch={'LeftVentriclePressure':'p_cav',
                         'LeftVentricleVolume':'V_cav',
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
        f=0.
        for n in range(len(self.costs)):
            func,var=self.decodeCost(self.costs[n])
            f+=func(folderName,runParameters,var)*self.weights[n]
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
            weights=[weights]
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
    def getspline(self,folderName,varName):
        if self.runMode=='LVbehaviorRun':
            fdatacl = np.loadtxt(folderName+"/circuit_results.txt",skiprows=1)
            name=[' ', 'p_cav',' ', 'V_cav',' ', 'PLA',' ','PRV',' ','PRA', ' ','Ppa1',' ','Paa' ]
            return scipyinterpolate.splrep(fdatacl[:,0]*1000., fdatacl[:,name.index(varName)])
        elif self.runMode=='iterativeRun':
            fdatacl = np.loadtxt(folderName+"/cl.txt")
            name=['tstep', 'p_cav', 'V_cav', 't_la', 'Pven', 'PLV', 'Part', 'PLA', 'Qla', 'Qao', 'Qper', 'Qmv', 'V_ven', 'V_cav', 'V_art', 'V_LA']
            return scipyinterpolate.splrep(fdatacl[:,0], fdatacl[:,name.index(varName)])
        elif self.runMode=='unloadedGeometry':
            fdatacl = np.loadtxt(folderName+"/"+varName+".txt")
            return fdatacl
    def MeanSquare(self,folderName,runParameters,varName):
        if varName[:6]=="coords":
            ref=np.loadtxt(self.matchingFolder+"/"+varName+".txt")
            simulated=self.getspline(folderName,varName)
            return np.mean((simulated-ref)**2.)*3.
        else:
            ref=np.loadtxt(self.matchingFolder+"/"+varName+".txt")
            spl=self.getspline(folderName,varName)
            ref[:,0]=ref[:,0]+runParameters['BCL']*int((spl[0][-6]-ref[0,0])/runParameters['BCL'])
            ref[:,0][ref[:,0]>spl[0][-7]]-=runParameters['BCL']
            simulated = scipyinterpolate.splev(ref[:,0], spl)
            return np.mean((simulated-ref[:,1])**2.)
    def getFunc(self,varName,collateFunction=None):
        if collateFunction is None:
            collateFunction=self.MeanSquare
        def func(folderName,runParameters):
            return collateFunction(folderName,runParameters,varName)
        return func
        
