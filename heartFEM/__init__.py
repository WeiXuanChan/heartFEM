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
from mpi4py import MPI as pyMPI

import lcleeHeart.
 
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
class LVclosed:
    def __new__(cls, *args, **kwargs):
        return super().__new__(cls)
    def __init__(self,defaultParameters=None,defaultAge='fetal'):
        
        #define defaults
        self.meshname = "t0"
        self.casename = os.getcwd()
        self.numofCycle=12
        self.LVage=defaultAge
        
        self.defaultParameters={}
        self.manualVol=None
        self.defaultParameters['topid'] = 4 #id for the base it is on top of geometry
        self.defaultParameters['endoid'] = 2 #id for the inner surface it is on top of geometry
        self.defaultParameters['epiid'] = 1 #id for the outer surface it is on top of geometry
        
        self.defaultParameters['Laxis_X']=None # longitudinal axis X
        self.defaultParameters['Laxis_Y']=None # longitudinal axis Y
        self.defaultParameters['Laxis_Z']=None # longitudinal axis Z
        #self.defaultParameters['endo_angle']=30#83 # endo fiber angle 
        #self.defaultParameters['epi_angle']=-50#-56 # epi fiber angle 
        
        self.defaultParameters['Kspring_constant']=90 #Kspring constant which model as pericardial cavity , unit in Pa
        
        self.defaultParameters['Tact_constant'] = 1e5 #Tact_constant we dont use it, unit is in Pa 
        self.defaultParameters['T0_LV'] = 60e3 #active tension forces unit in Pa
        
        self.defaultParameters['EDV_LV'] = 3.004778703264237
        self.defaultParameters['EDP_LV'] = None
        self.defaultParameters['ESV_LV'] = 1.4
        
        self.defaultParameters['BCL'] = 400.0 #set for the duration of cardiac cycle value is in ms, for fetal is 400ms. for adult is 800ms
        self.defaultParameters['lr'] = 1.85
        
        
        #Active Material
        if self.LVage=='adult':
            self.defaultParameters['Ca0'] = 4.35 #peak intracellular calcium concentration, µM
            self.defaultParameters['Ca0max'] = 4.35 #maximum peak intracellular calcium concentration, µM 
            self.defaultParameters['B'] = 4.75 #governs shape of peak isometric tension-sarcomere length relation, µm−1
            self.defaultParameters['t0'] = 222#238 #200.5#170 #132.5 #time to peak tension, ms
            self.defaultParameters['l0'] = 1.58#1.58 #sarcomere length at which no active tension develops,µm
            self.defaultParameters['m'] = 1048*0.5#1049 #slope of linear relaxation duration-sarcomere length relation, ms µm−1
            self.defaultParameters['b'] = -1600*0.5#-1429 #time-intercept of linear relaxation duration-sarcomere length relation, ms
        else:
            self.defaultParameters['Ca0'] = 4.35 #peak intracellular calcium concentration, µM
            self.defaultParameters['Ca0max'] = 4.35 #maximum peak intracellular calcium concentration, µM 
            self.defaultParameters['B'] = 4.75 #governs shape of peak isometric tension-sarcomere length relation, µm−1
            self.defaultParameters['t0'] = 222#238 #150.5#170 #132.5 #time to peak tension, ms
            self.defaultParameters['l0'] = 1.58#1.58 #sarcomere length at which no active tension develops,µm
            self.defaultParameters['m'] = 1048*0.5#1049 #slope of linear relaxation duration-sarcomere length relation, ms µm−1
            self.defaultParameters['b'] = -1600*0.5#0.5*(m_adult*l0+b_adult)-m_fetal #time-intercept of linear relaxation duration-sarcomere length relation, ms
        self.changeDefaultParameters(defaultParameters)
         
        self.constrains=['Kspring','translation','rotation']
        self.runCount=0
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
        if modelString=='adult':
            self.defaultParameters['AV'] = 160.0 #this is value set for left atrium which will be used in below 
        
            #### For Calculating P_LA ######################################## 
            self.defaultParameters['Ees_la'] = 60*1e0; #LCL #this is value for End-systolic elastance for left atrium Pa/m 
            self.defaultParameters['A_la'] = 58.05*10; #LCL #Scaling factor for EDPVR Pa
            self.defaultParameters['B_la'] = 0.049; #Exponent for EDPVR mL-1
            self.defaultParameters['V0_la'] = 0.25#1.0*5; #Volume axis intercept mL 
            self.defaultParameters['Tmax_la'] = 125; #Time to end-systole msec
            self.defaultParameters['tau_la'] = 25; #Time constant of relaxation msec
            
            self.defaultParameters['Cao'] = 0.0065*0.1 #this is value for compliance/capaciator of aortic valve  mLPa
            self.defaultParameters['Cven'] = 0.01*2 #this is value for compliance/capaciator of vein,mLPa
            self.defaultParameters['Vart0'] = 7.5#66*2;#450; #this is value for intial volume for artery  mL
            self.defaultParameters['Vven0'] = 96.5#34*10*2; #this is value for intial volume for vein mL
            self.defaultParameters['Rao'] = 10*5500*0.8*16*0.25*0.25*8*2*2*0.5;#10000.0; #this is value for resistance of aortic valve  ,unit is PamsecmL-1
            self.defaultParameters['Rven'] = 1000*0.5;#20000.0; #this is value for resistance of vein unit is PamsecmL-1
            self.defaultParameters['Rper'] = 100*140000*0.5*0.5*0.5*0.5*4*0.5*0.6*0.65*0.65*2*2#170000.0; #this is value for peripheral arteries resistance unit is PamsecmL-1
            self.defaultParameters['Rmv'] = 2500*2*2;#10000.0; #this is value for mitral valve unit is PamsecmL-1
            self.defaultParameters['V_ven'] = 92.5#37*10*2*12; #this is value for volume of vein mL
            self.defaultParameters['V_art'] = 18.5#74*2.5;#440 #this is value for volume of arteries mL
            self.defaultParameters['V_LA'] = 0.3#1.2*20; #this is value for volume of left atrium mL
        elif modelString[:5]=='fetal':
            
            self.defaultParameters['aac'] =0.05
            self.defaultParameters['ao1c'] =0.08
            self.defaultParameters['ao2c'] =0.07
            self.defaultParameters['ao3c'] =0.04
            self.defaultParameters['ao4c'] =0.05
            self.defaultParameters['pa1c'] =0.08
            self.defaultParameters['pa2c'] =0.08
            self.defaultParameters['lungc'] =0.4
            self.defaultParameters['cac'] =0.01
            self.defaultParameters['brc'] =0.3
            self.defaultParameters['ubc'] =0.85
            self.defaultParameters['intec'] =0.25
            self.defaultParameters['kidc'] =0.02
            self.defaultParameters['hec'] =3.
            self.defaultParameters['placc'] =1.5
            self.defaultParameters['legc'] =4.
            self.defaultParameters['uvc'] = 0.3
            self.defaultParameters['svcc'] =1.
            self.defaultParameters['ivcc'] =0.6
            self.defaultParameters['rac'] =1.
            self.defaultParameters['lac'] =2.

            self.defaultParameters['dal'] =0.006
            self.defaultParameters['ao1cal'] =0.08
            self.defaultParameters['aaao1l'] =0.002
            self.defaultParameters['pa1pa2l'] =0.02
            self.defaultParameters['raravalvl'] =0.0016
            self.defaultParameters['lalavalvl'] =0.0
            
            self.defaultParameters['aaao1r'] =0.12
            self.defaultParameters['ao1ao2r'] =0.4
            self.defaultParameters['ao2ao3r'] =0.04
            self.defaultParameters['ao3ao4r'] =0.06
            self.defaultParameters['pa1pa2r'] =0.07
            self.defaultParameters['pa2lungr'] =13.5
            self.defaultParameters['dar'] =0.01
            self.defaultParameters['ao1car'] =0.3
            self.defaultParameters['cabrr'] =3.
            self.defaultParameters['brsvcr'] =8.5
            self.defaultParameters['ao1ubr'] =8.
            self.defaultParameters['ubsvcr'] =4.9
            self.defaultParameters['ao3her'] =81.
            self.defaultParameters['ao3inter'] =34
            self.defaultParameters['inteher'] =7
            self.defaultParameters['ao3kidr'] =3.5
            self.defaultParameters['kidivcr'] =14
            self.defaultParameters['ao4placr'] =3.9
            self.defaultParameters['placuvr'] =3.4
            self.defaultParameters['ao4legr'] =3.5
            self.defaultParameters['legivcr'] =0.6
            self.defaultParameters['uvher'] =0.5
            self.defaultParameters['heivcr'] =0.16
            self.defaultParameters['dvr'] =1.3
            self.defaultParameters['svcrar'] =0.2
            self.defaultParameters['ivcrar'] =0.12
            self.defaultParameters['lunglar'] =2
            self.defaultParameters['for'] =0.
            self.defaultParameters['raravalvr'] =0.
            self.defaultParameters['rvrvvalvr'] =0.08
            self.defaultParameters['lvlvvalvr'] =0.
            self.defaultParameters['lalavalvr'] =0.
            
            self.defaultParameters['fok'] =0.4
            self.defaultParameters['dak'] =0.009
            self.defaultParameters['dvk'] =0.26
            self.defaultParameters['raravalvk'] =0.002
            self.defaultParameters['rvrvvalvk'] =0.001
            self.defaultParameters['lalavalvk'] =0.002
            self.defaultParameters['lvlvvalvk'] =0.001
            
            self.defaultParameters['fob'] =0.625
            self.defaultParameters['dab'] =2.
            self.defaultParameters['dvb'] =2.
            self.defaultParameters['raravalvb'] =2.
            self.defaultParameters['rvrvvalvb'] =2.
            self.defaultParameters['lalavalvb'] =2.
            self.defaultParameters['lvlvvalvb'] =2.
            try:
                self.setRVcurrentFourier(self.casename+'/'+self.meshname+'_rvflowrate.txt')
            except Exception as e:
                print(e)
                print('RV fourier not loaded',self.casename+'/'+self.meshname+'_rvflowrate.txt')
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
        #scaleDict={'default':,'aac':,...}
        if 'default' not in scaleDict:
            scaleDict['default']=1.
        if compstr!='':
            for comp in WindkesselComponents + WindkessellinkComponents:
                if comp+compstr in self.defaultParameters:
                    if comp+compstr in scaleDict:
                        self.defaultParameters[comp+compstr]*=scaleDict[comp+compstr]
                    elif compstr in scaleDict:
                        self.defaultParameters[comp+compstr]*=scaleDict[compstr]
                    else:
                        self.defaultParameters[comp+compstr]*=scaleDict['default']
        else:
            self.scaleWinkessel(scaleDict,compstr='r')
            self.scaleWinkessel(scaleDict,compstr='l')
            self.scaleWinkessel(scaleDict,compstr='c')
            self.scaleWinkessel(scaleDict,compstr='k')
            self.scaleWinkessel(scaleDict,compstr='b')
    def scaleWinkesselwithAge(self,ageInWeeks,poweradjustDict=None,compstr=''):
        if poweradjustDict is None:
            poweradjustDict={}
        if compstr!='':
            for comp in WindkesselComponents + WindkessellinkComponents:
                if comp+compstr in self.defaultParameters:
                    if comp+compstr in defaultAgeScalePower:
                        agepower=defaultAgeScalePower[comp+compstr]
                    else:
                        agepower=defaultAgeScalePower['default'+compstr]
                    if comp+compstr in poweradjustDict:
                        agepower-=defaultAgeScalePower['default'+compstr]+poweradjustDict[comp+compstr]*defaultAgeScalePower['default'+compstr]
                    elif compstr in poweradjustDict:
                        agepower-=defaultAgeScalePower['default'+compstr]+poweradjustDict[compstr]*defaultAgeScalePower['default'+compstr]
                    self.defaultParameters[comp+compstr]*=((10.**(0.2508+0.1458*ageInWeeks-0.0016*ageInWeeks**2.))/(10.**(0.2508 + 0.1458*38.-0.0016*38.**2.)))**agepower
        else:
            self.scaleWinkesselwithAge(ageInWeeks,poweradjustDict=poweradjustDict,compstr='r')
            self.scaleWinkesselwithAge(ageInWeeks,poweradjustDict=poweradjustDict,compstr='l')
            self.scaleWinkesselwithAge(ageInWeeks,poweradjustDict=poweradjustDict,compstr='c')
            self.scaleWinkesselwithAge(ageInWeeks,poweradjustDict=poweradjustDict,compstr='k')
            self.scaleWinkesselwithAge(ageInWeeks,poweradjustDict=poweradjustDict,compstr='b')
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
        pdata = lcleeHeart..readSTL(meshfilename)
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
        print('Estimated Laxis as ',Laxis)
    def generateMesh(self,Laxis,endo_angle='',epi_angle='',clipratio=0.95,saveaddstr=''):
        return lcleeHeart..generateMesh(self.casePath,self.stlname,Laxis,endo_angle=endo_angle,epi_angle=epi_angle,clipratio=clipratio,meshname=self.stlname+saveaddstr)
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
        		print ("out=", out)
        	Part = 1.0/Cao*(V_art - Vart0);
        	Pven = 1.0/Cven*(V_ven - Vven0);
        	PLV = p_cav;
            
        	if (tstep < (BCL - AV)):
        		t_la = tstep + AV;
        	else: 
        		t_la = tstep - BCL + AV;
        
        	if(fenics.MPI.rank(comm) == 0):
        		print ("t_LA = ", t_la)
        
        		PLA = et(t_la,Tmax_la,tau_la)*Ees_la*(V_LA - V0_la) + (1 - et(t_la,Tmax_la,tau_la))*A_la*(math.exp(B_la*(V_LA - V0_la)) - 1);
        	##################################################################################################################################
        
        	if(fenics.MPI.rank(comm) == 0):
        		print ("P_ven = ",Pven);
        		print ("P_LV = ", PLV);
        		print ("P_art = ", Part);		
        		print ("P_LA = ", PLA);
        
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
        		print ("Q_LA = ", Qla) ;
        		print ("Q_ao = ", Qao) ;
        		print ("Q_per = ", Qper);
        		print ("Q_mv = ", Qmv) ;
        
        	V_cav_prev = V_cav
        	V_art_prev = V_art
        	V_ven_prev = V_ven
        
        	V_cav = V_cav + dt_dt*(Qla - Qao);
        	V_art = V_art + dt_dt*(Qao - Qper);
        	V_ven = V_ven + dt_dt*(Qper - Qmv);
        	V_LA = V_LA + dt_dt*(Qmv - Qla);
                    
        	if(fenics.MPI.rank(comm) == 0):
        		print ("V_ven = ", V_ven);
        		print ("V_LV = ", V_cav);
        		print ("V_art = ", V_art);
        		print ("V_LA = ", V_LA);
                
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
    def getLVbehavior(self,runParameters=None,volstep=50,extendVol=0.1,minESV=None,maxEDV=None,saveaddstr='',voladj=0.1,starttime=0,endtime=None):
        if runParameters is None:
            runParameters=self.defaultParameters
        solver,uflforms,activeforms,t_a,dt,Cavityvol,f0,mesh_volume,mesh,comm=lcleeHeart..setupHeart(self.casename,self.meshname,runParameters)
        
        ########################### Fenics's Newton  #########################################################
        Cavityvol.vol = uflforms.cavityvol()
        
        # Closed loop cycle
        BCL = runParameters['BCL']
        
        
        #######set volume solve space#############
        EDV=runParameters['EDV_LV']
        if maxEDV is None:
            maxEDV=EDV*(1+extendVol)
        if 'ESV_LV' in runParameters:
            ESV=runParameters['ESV_LV']
        else:
            ESV=Cavityvol.vol
        if minESV is None:
            minESV=ESV*(1-extendVol)
        volumeSpace=minESV*(maxEDV/minESV)**(np.linspace(0,1,num=volstep))[::-1]
        ####### Loading phase to MAX volume ####################################################

        if Cavityvol.vol>(volumeSpace[0]):
        	tuneVol=-1
        else:
            tuneVol=1
        while (Cavityvol.vol*tuneVol)<(volumeSpace[0]*tuneVol):
        	if (((1+tuneVol*voladj)*Cavityvol.vol)*tuneVol)<(volumeSpace[0]*tuneVol):
        		Cavityvol.vol *= (1+tuneVol*voladj)
        	else:
        		Cavityvol.vol = volumeSpace[0]
        	solver.solvenonlinear()
        	p_cav = uflforms.cavitypressure()
        	V_cav = uflforms.cavityvol()
        	print ("PLV = ", p_cav*.0075, " VLV = ", V_cav)	

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
            timeSpace=timeSpace[timeSpace<=endtime]
        if(fenics.MPI.rank(comm) == 0):
            np.savetxt(self.casename+"/"+self.meshname+"_Press_volumeSpace"+saveaddstr+".txt",np.array(volumeSpace))
            np.savetxt(self.casename+"/"+self.meshname+"_Press_timeSpace"+saveaddstr+".txt",np.array(timeSpace))
        print('========START============')
        press_volTime_base=np.zeros(len(volumeSpace))
        press_volTime=np.zeros((len(volumeSpace),len(timeSpace)))
        maxtimeind=0
        for curr_volN in range(len(volumeSpace)):
            t_a.t_a=timeSpace[0]
            dt.dt=timeSpace[1]-timeSpace[0]
            press_volTime_temp=[]
            samePressure=0
            while (Cavityvol.vol*(1-voladj))>volumeSpace[curr_volN]:
            	Cavityvol.vol *= (1-voladj)
            	solver.solvenonlinear()
            	p_cav = uflforms.cavitypressure()
            	V_cav = uflforms.cavityvol()
            	print ("PLV = ", p_cav*.0075, " VLV = ", V_cav)	
            for curr_timeN in range(len(timeSpace)):
            	t = timeSpace[curr_timeN]
            	if curr_timeN<(len(timeSpace)-1):
            		dt.dt = timeSpace[curr_timeN+1]-timeSpace[curr_timeN]
            	else :
            		dt.dt = 1.0
            	t_a.t_a = t
            	print (curr_volN,'/',len(volumeSpace),curr_timeN,'/',len(timeSpace)," t_a=",t_a)

            	Cavityvol.vol = volumeSpace[curr_volN]
            
            	solver.solvenonlinear()
            	print ("PLV = ", uflforms.cavitypressure()*.0075, " VLV = ", uflforms.cavityvol())
            	press_volTime_temp.append(uflforms.cavitypressure()*.0075)
            	if len(press_volTime_temp)>1:
                    if abs(press_volTime_temp[-1]-press_volTime_temp[-2])<10**-12:
                        samePressure+=1
                    else:
                        samePressure=0
                    if samePressure>1:
                        break
            if len(press_volTime_temp)>maxtimeind:
                maxtimeind=len(press_volTime_temp)
            press_volTime_base[curr_volN]=press_volTime_temp[-1]
            press_volTime[curr_volN,:len(press_volTime_temp)]=np.array(press_volTime_temp)-press_volTime_base[curr_volN]
            if(fenics.MPI.rank(comm) == 0):
                np.savetxt(self.casename+"/"+self.meshname+"_Press_VolTime"+saveaddstr+".txt",press_volTime,header=str(curr_volN)+'solved rowVolm: '+str(volumeSpace)+'\ncolTime: '+str(timeSpace))
                np.savetxt(self.casename+"/"+self.meshname+"_Press_VolTime_base"+saveaddstr+".txt",press_volTime_base)
        # press_volTime=press_volTime[:,:maxtimeind]
        #timeSpace=timeSpace[:maxtimeind]
        if(fenics.MPI.rank(comm) == 0):
            #np.savetxt(self.casename+"/"+self.meshname+"_Press_timeSpace"+saveaddstr+".txt",np.array(timeSpace))
            np.savetxt(self.casename+"/"+self.meshname+"_Press_VolTime"+saveaddstr+".txt",np.array(press_volTime),header='rowVolm: '+str(volumeSpace)+'\ncolTime: '+str(timeSpace))
            np.savetxt(self.casename+"/"+self.meshname+"_Press_VolTime_base"+saveaddstr+".txt",press_volTime_base)
    def adjustMeshToPressure():
        return
        
    def __call__(self,editParameters=None,always_remesh=False,runTimeList=None,getEmat=False):
        runParameters=self.defaultParameters.copy()
        if editParameters is not None:
            for key in runParameters:
                if key in editParameters:
                    runParameters[key]=editParameters[key]
        else:
            editParameters={}
        if isinstance(self.runCount,int):
            self.runCount+=1
        os.makedirs(self.casename+"/"+str(self.runCount),exist_ok=True)
        os.makedirs(self.casename+"/"+str(self.runCount)+"/stress",exist_ok=True)
        with open(self.casename+'/'+str(self.runCount)+"/runParameters.txt", "w") as f:
                for key, val in runParameters.items():
                    f.write("{0:s} , {1:s}\n".format(key,str(val)))
        pressureIteration_factor=0.1
        
        displacementfile = fenics.File(self.casename+"/"+str(self.runCount)+"/deformation/u_disp.pvd")
        activestressfile = fenics.File(self.casename+"/"+str(self.runCount)+"/activestress/stress.pvd") #CW
        stress_File = fenics.File(self.casename+"/"+str(self.runCount)+"/stress/_stress.pvd") #Joy Changed here
        strain_File = fenics.File(self.casename+"/"+str(self.runCount)+"/strain/_strain.pvd") #Laura Changed here
        ECCstrain_File = fenics.File(self.casename+"/"+str(self.runCount)+"/strain/_ECCstrain.pvd") #Laura Changed here
        ELLstrain_File = fenics.File(self.casename+"/"+str(self.runCount)+"/strain/_ELLstrain.pvd") #Laura Changed here

        while True: #setup
            if always_remesh or str(self.runCount)=='1' or 'Laxis_X' in editParameters or 'Laxis_Y' in editParameters or 'Laxis_Z' in editParameters or 'endo_angle' in editParameters or 'epi_angle' in editParameters:
                if runParameters['Laxis_X'] is None or runParameters['Laxis_Y'] is None or runParameters['Laxis_Z'] is None:
                    self.getLongAxis()
                    runParameters['Laxis_X']=self.defaultParameters['Laxis_X']
                    runParameters['Laxis_Y']=self.defaultParameters['Laxis_Y']
                    runParameters['Laxis_Z']=self.defaultParameters['Laxis_Z']
                self.generateMesh(np.array([runParameters['Laxis_X'],runParameters['Laxis_Y'],runParameters['Laxis_Z']]),runParameters['endo_angle'] ,runParameters['epi_angle'])
            solver,uflforms,activeforms,t_a,dt,Cavityvol,f0,mesh_volume,mesh,comm=lcleeHeart..setupHeart(self.casename,self.meshname,runParameters)
            
            ########################### Fenics's Newton  #########################################################
            Cavityvol.vol = uflforms.cavityvol()
            
            # Closed loop cycle
            BCL = runParameters['BCL']
            if self.LVage=='adult':
                V_ven = runParameters['V_ven']
                V_art = runParameters['V_art']
                V_LA = runParameters['V_LA']
  
            cycle = 0.
            t = 0
            tstep = 0
            dt.dt = 1.
            
            t_array=[0]
            LVcav_array = [uflforms.cavityvol()]
            Pcav_array = [uflforms.cavitypressure()*0.0075]
            
            #Joy Changed from here
            #Stress calculation 
            cauchy1 =  uflforms.Cauchy1() + activeforms.cauchy()
            cauchy = fenics.project(cauchy1,fenics.TensorFunctionSpace(mesh, "DG", 1), form_compiler_parameters={"representation":"uflacs"})
            cauchy.rename("cauchy_stress","cauchy_stress")
            ###stress_File << cauchy ## wei xuan remove

            ####### Loading phase for LV ####################################################
            current_starting_cavity_volume=Cavityvol.vol
            if self.manualVol is None:
                startVolume=runParameters['EDV_LV']
            elif runTimeList is None:
                startVolume=self.manualVol(0.)
            else:
                startVolume=self.manualVol(runTimeList[0])
            #if Cavityvol.vol>(1.1*runParameters['EDV_LV']):
            # 	raise Exception('Starting cavity volume,'+str(Cavityvol.vol)+', is more than target EDV_LV,'+str(runParameters['EDV_LV']))
            while Cavityvol.vol<startVolume:
            	if (1.1*Cavityvol.vol)<startVolume:
            		Cavityvol.vol *= 1.1
            	else:
            		Cavityvol.vol = startVolume
            	solver.solvenonlinear()
            	p_cav = uflforms.cavitypressure()
            	V_cav = uflforms.cavityvol()
            	print ("PLV = ", p_cav*.0075, " VLV = ", V_cav)
                
            while Cavityvol.vol>startVolume:
            	if (0.9*Cavityvol.vol)<startVolume:
            		Cavityvol.vol *= 0.9
            	else:
            		Cavityvol.vol = startVolume
            	solver.solvenonlinear()
            	p_cav = uflforms.cavitypressure()
            	V_cav = uflforms.cavityvol()
            	print ("PLV = ", p_cav*.0075, " VLV = ", V_cav)
            
                   
            if runParameters['EDP_LV'] is None or self.manualVol is not None:
                break
            else:
                ###need to edit
            	p_cav = uflforms.cavitypressure()
            	if np.abs((p_cav-runParameters['EDP_LV'])/runParameters['EDP_LV'])>10**-3 and abs(pressureIteration_factor)>10**-3:
            		meshfile=self.casename+ '/'+self.meshname+ '.stl'
            		if not(os.path.isfile(self.casename+ '/'+self.meshname+ '_originalBackup.stl')):
            			cmd = "cp " + meshfile+' '+ self.casename+ '/'+self.meshname+ '_originalBackup.stl'
            			os.system(cmd)
            		if p_cav<runParameters['EDP_LV']:
            		    if pressureIteration_factor>0:
            		    	pressureIteration_factor*=-0.3
            		else:
            		    if pressureIteration_factor<0:
            		    	pressureIteration_factor*=-0.3
            		stl_scale_factor=1+pressureIteration_factor
            		if stl_scale_factor<=0:
            		    raise Exception('Unable to reach target pressure')
            		pdata = lcleeHeart.readSTL(meshfilename)
            		pdata = lcleeHeart.scale_pdata(stl_scale_factor)
            		lcleeHeart.writeSTL(pdata, meshfile)
            	else:
            		break

        if(fenics.MPI.rank(comm) == 0): #added to output the volume, stroke volume data before the while loop
            fdataPV = open(self.casename+"/"+str(self.runCount)+"/PV_.txt", "w")
            fdatals = open(self.casename+"/"+str(self.runCount)+"/ls.txt", "w")
            fdataSV = open(self.casename+"/"+str(self.runCount)+"/SVp.txt", "w")
            fdataSVt = open(self.casename+"/"+str(self.runCount)+"/SVpt.txt", "w")
            fdata_stress = open( self.casename+"/"+str(self.runCount)+"/stress/_stress_.txt", "w")

        p_cav = uflforms.cavitypressure()
        V_cav = uflforms.cavityvol()
        if(fenics.MPI.rank(comm) == 0):
        	print("Cycle number = ", cycle, " cell time = ", t, " tstep = ", tstep, " dt = ", dt.dt)
        	print(tstep, p_cav*0.0075 , V_cav, file=fdataPV)


        if(fenics.MPI.rank(comm) == 0):
            displacementfile << w.sub(0)
        #Joy changed  from here
        cauchy = project(cauchy1,TensorFunctionSpace(mesh, "DG", 1), form_compiler_parameters={"representation":"uflacs"})
        cauchy.rename("cauchy_stress","cauchy_stress")
       	#print (cauchy)
       	#print (cauchy.vector().array()[:])
        stress_File << cauchy

        if getEmat is True:
            Emat =  uflforms.Emat()
            Emat = project(Emat,TensorFunctionSpace(mesh, "DG", 1), form_compiler_parameters={"representation":"uflacs"})
            Emat.rename("green_strain","green_strain")
            strain_File << Emat

            EmatECC =  uflforms.EmatECC()
            EmatECC = project(EmatECC,FunctionSpace(mesh, "DG", 1), form_compiler_parameters={"representation":"uflacs"})
            EmatECC.rename("ECCgreen_strain","ECCgreen_strain")
            ECCstrain_File << EmatECC

            EmatELL =  uflforms.EmatELL()
            EmatELL = project(EmatELL,FunctionSpace(mesh, "DG", 1), form_compiler_parameters={"representation":"uflacs"})
            EmatELL.rename("ELLgreen_strain","ELLgreen_strain")
            ELLstrain_File << EmatELL
            
            #reader = vtk.vtkXMLUnstructuredGridReader()
            #reader.SetFileName(strain_File)
            #reader.Update()  # Needed because of GetScalarRange
            #output = reader.GetOutput()
            #potential = output.GetPointData().GetArray("potential")

        sigma_fiber_LV = activeforms.CalculateFiberStress(sigma = cauchy, e_fiber = f0, Vol = mesh_volume, Mesh = mesh)
       
        #if(fenics.MPI.rank(comm) == 0):
       	#	print (fdata_stress, 0.0, sigma_fiber_LV)
       	#Joy chnaged till here
        #save files
       
        # Closed-loop phase
        if self.LVage!='adult' and self.manualVol is None:
            ngspice_py.createLVcircuit(self.casename+'/'+self.meshname,runParameters)
        if runTimeList is not None:
            runTimeList.pop(0)
        while(cycle < self.numofCycle):
        	
        	#################################################################
        	p_cav = uflforms.cavitypressure()
        	V_cav = uflforms.cavityvol()
        	if runTimeList is None:
        		tstep = tstep + dt.dt
        	elif len(runTimeList)<=0:
        		break
        	else:
        		tstep =runTimeList.pop(0)
        	cycle = math.floor(tstep/BCL)
        	print ("cycle=",cycle)
        	t = tstep - cycle*BCL
        
        	if(t >= 0.0 and t < 4.0):
        		dt.dt = 0.50
        	elif (t >= 200.0 and t < 201.0):
        		dt.dt = 3
        	else :
        		dt.dt = 1.0
        
        	t_a.t_a = t
        	print ("t_a=",t_a)
        
        	
        	if self.manualVol is not None:
        		V_cav=self.manualVol(tstep)
        	elif self.LVage=='adult':
        		V_cav,V_art,V_ven,V_LA=self.runWinkessel(comm,tstep,dt.dt,p_cav,V_cav,V_art,V_ven,V_LA)
        	else:
        		V_cav=self.runWinkessel(comm,tstep,dt.dt)

            
        	t_array.append(t)
        	LVcav_array.append(V_cav)
        	Pcav_array.append(p_cav*0.0075)
            
            #these lines added as between the two timesteps if there is a large volume it will crash, therefore below increases volume gradually
           

        	while Cavityvol.vol<V_cav:
        		if (1.1*Cavityvol.vol)<V_cav:
        			Cavityvol.vol *= 1.1
        		else:
        			Cavityvol.vol = V_cav
        		solver.solvenonlinear()
        		print ("PLV = ", p_cav*.0075, " VLV = ", V_cav)
 
        	while Cavityvol.vol>V_cav:
        		if (0.9*Cavityvol.vol)<V_cav:
        			Cavityvol.vol *= 0.9
        		else:
        			Cavityvol.vol = V_cav
        		solver.solvenonlinear()
        		print ("PLV = ", p_cav*.0075, " VLV = ", V_cav)

            
        	p_cav = uflforms.cavitypressure()
        	V_cav = uflforms.cavityvol()

        	if(fenics.MPI.rank(comm) == 0):
        		print("Cycle number = ", cycle, " cell time = ", t, " tstep = ", tstep, " dt = ", dt.dt)
        		print(tstep, p_cav*0.0075 , V_cav, file=fdataPV)
        	

        	ls = sqrt(dot(f0, Cmat*f0))*ls0
        	ls1 = project(ls,Q).vector().get_local()[:]
        	eca = project(activeforms.ECa(), Q).vector().get_local()[:]
        	t_r = project(activeforms.tr(), Q).vector().get_local()[:]
        
        	if(fenics.MPI.rank(comm) == 0):
        		print(t_a.t_a, min(ls1), max(ls1), min(eca), max(eca), min(t_r), max(t_r), file=fdatals)

        	if(fenics.MPI.rank(comm) == 0):
        		displacementfile << w.sub(0)
        	#Joy Chnaged from here
        	cauchy = project(cauchy1,TensorFunctionSpace(mesh, "DG", 1), form_compiler_parameters={"representation":"uflacs"})
        	cauchy.rename("cauchy_stress","cauchy_stress")
        	print (cauchy)
        #	print (cauchy.vector().array()[:])
        	stress_File << cauchy
  
        	if getEmat is True:
        		Emat = uflforms.Emat()
        		Emat = project(Emat,TensorFunctionSpace(mesh, "DG", 1), form_compiler_parameters={"representation":"uflacs"})
        		Emat.rename("green_strain","green_strain")
        		strain_File << Emat

        		EmatECC =  uflforms.EmatECC()
        		EmatECC = project(EmatECC,FunctionSpace(mesh, "DG", 1), form_compiler_parameters={"representation":"uflacs"})
        		EmatECC.rename("ECCgreen_strain","ECCgreen_strain")
        		ECCstrain_File << EmatECC

        		EmatELL =  uflforms.EmatELL()
        		EmatELL = project(EmatELL,FunctionSpace(mesh, "DG", 1), form_compiler_parameters={"representation":"uflacs"})
        		EmatELL.rename("ELLgreen_strain","ELLgreen_strain")
        		ELLstrain_File << EmatELL


        	sigma_fiber_LV = activeforms.CalculateFiberStress(sigma = cauchy, e_fiber = f0, Vol = mesh_volume, Mesh = mesh)
        
        	if(fenics.MPI.rank(comm) == 0):
        		print (fdata_stress, tstep, sigma_fiber_LV, file=fdata_stress ) 
        	#Joy changed Till here          
        
        if(fenics.MPI.rank(comm) == 0):
        	fdataPV.close()
        
        if (tstep >= runParameters['BCL']):
        	t_array=np.array(t_array)
        	temp_Ind=np.nonzero((t_array[2:]-t_array[:-2])>0)[0][-1]
        	temp_ind=np.nonzero(np.logical_and(t_array[:temp_Ind]>t_array[temp_Ind],t_array[:temp_Ind]<t_array[temp_Ind+2]))[0][-1]

        	#print >>fdataSV,tstep, p_cav*0.0075, V_cav
        	LVcav_arrayT=LVcav_array[(temp_ind-1):]
        	Pcav_arrayT=Pcav_array[(temp_ind-1):]
        	plt.plot(LVcav_arrayT, Pcav_arrayT)
        	plt.xlabel("Volume")
        	plt.ylabel("Pressure")
        	plt.legend()
        	plt.savefig(self.casename+"/PV.png")
        	
        	print(max(LVcav_arrayT),min(LVcav_arrayT),max(LVcav_arrayT)-min(LVcav_arrayT), max(Pcav_arrayT),min(Pcav_arrayT),max(Pcav_arrayT)-min(Pcav_arrayT), file=fdataSV) 
        	print(LVcav_arrayT,Pcav_arrayT, file=fdataSVt)
        
        if(fenics.MPI.rank(comm) == 0):
        	fdataSV.close()
        	fdataSVt.close()
        	fdata_stress.close() #Joy Changed here	
        ######################################################################################################
        return runParameters

    # def getaverageLVstrain(runTimeList=None):
    #     uflforms = Forms(params)
    #     Emat = uflforms.Emat()
    #     eCC=numpy.dot(numpy.array(Emat)[:,0], eC)
    #     print(eCC)
        
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
        os.makedirs(self.casename+'/best_fit',exist_ok=True)
        editParameters={}
        for n in self.variableList:
            editParameters[self.variableList[n]]=para[n]
        try:
            runParameters=self.runSimulation(editParameters=editParameters)           
        except Exception as e_inst:
            with open(self.casename+'/'+str(self.runCount)+"/runParameters.txt", "r+") as f:
                 old = f.read() # read everything in the file
                 f.seek(0) # rewind
                 f.write("Status , FAILED\n cost , "+str(type(e_inst))+'  '+str(type(e_inst))+"\n" + old) # write the new line before
            return float('inf')
        else:
            cost_current=self.cost_function(self.runSimulation.casename+'/'+str(self.runSimulation.runCount),runParameters)
            with open(self.casename+'/'+str(self.runCount)+"/runParameters.txt", "r+") as f:
                 old = f.read() # read everything in the file
                 f.seek(0) # rewind
                 f.write("Status , FAILED\n cost , "+str(cost_current)+"\n" + old) # write the new line before
            if cost_current<self.minCost:
                self.minCost=cost_current
                cmd = "cp -r " + self.casename+'/'+str(self.runCount)+'/* '+ self.casename+'/best_fit/summary'
                os.system(cmd)
                cmd = "cp -r " + self.casename+'/activestress '+ self.casename+'/best_fit'
                os.system(cmd)
                cmd = "cp -r " + self.casename+'/deformation '+ self.casename+'/best_fit'
                os.system(cmd)
                cmd = "cp -r " + self.casename+'/stress '+ self.casename+'/best_fit'
                os.system(cmd)
                cmd = "cp " + self.casename+"/PV.png "+ self.casename+'/best_fit'
                os.system(cmd)
            return cost_current
class cost_function:
    def __new__(cls, *args, **kwargs):
        return super().__new__(cls)
    def __init__(self,matchingFolder='ref'):
        self.costs=[]
        self.weights=[]
        if os.path.isdir(matchingFolder):
            self.matchingFolder=matchingFolder
        else:
            self.matchingFolder=os.getcwd()+'/'+matchingFolder
        self.costSwitch={'CavityPressureMeanSquare':self.getFunc('p_cav'),
                         'CavityVolumeMeanSquare':self.getFunc('V_cav'),
                         't_laMeanSquare':self.getFunc('t_la'),
                         'VentriclePressureMeanSquare':self.getFunc('Pven'),
                         'LVPressureMeanSquare':self.getFunc('PLV'),
                         'ArterialPressureMeanSquare':self.getFunc('Part'), 
                         'LAPressureMeanSquare':self.getFunc('PLA'), 
                         'QlaMeanSquare':self.getFunc('Qla'), 
                         'AortaQMeanSquare':self.getFunc('Qao'), 
                         'perQMeanSquare':self.getFunc('Qper'), 
                         'mvQMeanSquare':self.getFunc('Qmv'), 
                         'VentricleVolumeMeanSquare':self.getFunc('V_ven'), 
                         'CatvityVolumeMeanSquare':self.getFunc('V_cav'), 
                         'ArterialVolumeMeanSquare':self.getFunc('V_art'), 
                         'LAVolumeMeanSquare':self.getFunc('V_LA')
                         }
    def __call__(self,folderName,runParameters):
        f=0.
        for n in range(len(self.costs)):
            f+=self.costs[n](folderName)*self.weights[n]
        return f
    def addCost(self,costString,weights=None):
        #put weight to negative to cost which is better if higher
        if weights is None:
            weights=1.
        if isinstance(costString,str):
            costString=[costString]
        if isinstance(weights,(float,int)):
            weights=[weights]
        for n in range(len(costString)):
            self.costs+=self.costSwitch[costString[n]]
            self.weights+=weights[n]
    def removeCost(self,costString):
        if isinstance(costString,str):
            costString=[costString]
        for string in costString:
            if string in self.costs:
                temp_ind=self.costs.index(string)
                self.costs.pop(temp_ind)
                self.weights.pop(temp_ind)
    def getspline(self,folderName,varName):
        fdatacl = np.loadtxt(folderName+"/cl.txt")
        name=['tstep', 'p_cav', 'V_cav', 't_la', 'Pven', 'PLV', 'Part', 'PLA', 'Qla', 'Qao', 'Qper', 'Qmv', 'V_ven', 'V_cav', 'V_art', 'V_LA']
        return scipyinterpolate.splrep(fdatacl[:,0], fdatacl[:,name.index(varName)])
    def MeanSquare(self,folderName,runParameters,varName):
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
        
