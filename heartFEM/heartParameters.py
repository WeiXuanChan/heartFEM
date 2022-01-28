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
  Author: w.x.chan@gmail.com         05APR2021           - Created
  Author: w.x.chan@gmail.com         26APR2021           - v1.2.0
                                                            -added endo_angle and epi_angle
  Author: w.x.chan@gmail.com         19MAY2021           - v1.2.1
                                                            -debug func changeDefaultParameters
                                                            -debug when defaultParameters and defaultAge is swapped during initialization
  Author: w.x.chan@gmail.com         19MAY2021           - v3.0.0
                                                            -added reguritation for heart chambers
                                                            - added vra0 and vla0
                                                            -debug these values: 'rac', 'lac', 'pa1pa2l', 'lalavalvl', 'lvlvvalvr', 
                                                            - debug scale with age
  Author: w.x.chan@gmail.com         08JUL2021           - v3.1.0
                                                            -added 'outOfplaneAngle'
  Author: w.x.chan@gmail.com         14JUL2021           - v3.3.0
                                                            -removed outOfplaneDeg and added 'fiberSheetletAngle','fiberSheetletWidth','radialFiberAngle''Windkessel_scale_T0_LV'
  Author: w.x.chan@gmail.com         12Aug2021           - v3.5.0
                                                            -added 'Windkessel_scale_T0_LV'
                                                            -added Aortic_stenosis to multiple resistances from LV to AO
  Author: w.x.chan@gmail.com         20Sep2021           - v3.5.1
                                                            -added switchvalve
  Author: w.x.chan@gmail.com         29Sep2021           - v3.6.0
                                                            -added regurgitation at aortic and pulmonary valve
                                                            -added "StrainEnergyDensityFunction_Coef"
  Author: w.x.chan@gmail.com         05Nov2021           - v4.0.0
                                                            -rename all parameters
                                                            - change time units to milliseconds for Windkessel
                                                            - adjustment to default power scaling
'''
_version='3.6.0'
import logging
import numpy as np
import os
try:
    import dolfin as fenics
except:
    import fenics as fenics
logger = logging.getLogger(__name__)

parameters_for_FEniCS=['spring constant for LV pericardial cavity in Pa',
                       'maximum LV fiber tension in Pa',
                       'LV end systolic volume in mL',
                       'LV sarcomere length',
                       'duration of one cardiac cycle in ms',
                       'peak intracellular calcium concentration in uM',
                       'maximum peak intracellular calcium concentration in uM',
                       'exponential coefficient for relation of peak isometric tension and sarcomere length in um-1',
                       'peak intracellular calcium concentration in uM',
                       'LV sarcomere stretch ratio threshold where no tension develops',
                       'slope of linear relation of relaxation duration and sarcomere length in ms um-1',
                       'time intercept of linear relation of relaxation duration and sarcomere length in ms']                 
parameters_for_mesh=['LV base mesh surface ID',
                     'LV endocardial mesh surface ID',
                     'LV epicardial mesh surface ID',
                     'LV longitudinal axis superior X',
                     'LV longitudinal axis superior Y',
                     'LV longitudinal axis superior Z',
                     'ratio to clip mesh',
                     'LV endocardial fiber angle in degree',
                     'LV epicardial fiber angle in degree',
                     'LV fiber sheetlet angle in degrees',
                     'LV fiber sheetlet width',
                     'LV fiber radial angle in degrees',
                     'LV end diastolic volume in mL',
                     'LV end diastolic pressure in mmHg']
WindkesselComponents=['lv','la','rv','ra','aa','ao1','ao2','ao3','ao4','br','ca','ub','he','inte','ivc','kid','leg','lung','pa1','pa2','plac','svc','uv']
WindkessellinkComponents=['aaao1','ao1ao2','ao2ao3','ao3ao4','pa1pa2','pa2lung','da','ao1ca','cabr','brsvc','ao1ub','ubsvc','ao3he','ao3inte','intehe','ao3kid','kidivc','ao4plac','placuv','ao4leg','legivc','uvhe','heivc','dv','svcra','ivcra','lungla','fo','raravalv','lalavalv','lvlvvalv','rvrvvalv']
defaultAgeScalePower={'defaultr':-1.,'pa2lungr':-1.2,'lunglar':-1.2,'cabrr':-1.1,'brsvcr':-1.1,'dvr':-0.55,
                      'defaultl':-0.33,
                      'defaultc':1.,'brc':1.471,'lungc':1.6,'rac':0.5,'lac':0.5,
                      'defaultk':0.,'fok':-0.6,'dak': -2.5,'dvk':-0.88,'raravalvk':-1.33,'rvrvvalvk':-1.33,'lalavalvk':-1.33,'lvlvvalvk':-1.33,
                      'defaultb':0.}
defaultAdjustmentToScaling={'r_scale':1.21,'c_scale':0.27}



class heartParameters(dict):
    def __new__(cls, *args, **kwargs):
        return super().__new__(cls)
    def __init__(self,defaultParameters=None,defaultAge=None):
        '''
        NGSPICE parameters:
            'ES_time'          # target End-systolic time to match while adjusting time to peak tension
            ADULT
              FOR calculating loading phase
                'Ees_la'      #LCL #this is value for End-systolic elastance for left atrium Pa/m 
                'A_la'        #LCL #Scaling factor for EDPVR Pa
                'B_la'        #Exponent for EDPVR mL-1
                'V0_la'       #Volume axis intercept mL 
                'Tmax_la'     #Time to end-systole msec
                'tau_la'      #Time constant of relaxation msec
              'Cao'      #this is value for compliance/capaciator of aortic valve  mLPa
              'Cven'     #this is value for compliance/capaciator of vein,mLPa
              'Vart0'    #this is value for intial volume for artery  mL
              'Vven0'    #this is value for intial volume for vein mL
              'Rao'      #this is value for resistance of aortic valve  ,unit is PamsecmL-1
              'Rven'     #this is value for resistance of vein unit is PamsecmL-1
              'Rper'     #this is value for peripheral arteries resistance unit is PamsecmL-1
              'Rmv'      #this is value for mitral valve unit is PamsecmL-1
              'V_ven'    #this is value for volume of vein mL
              'V_art'    #this is value for volume of arteries mL
              'V_LA'     #this is value for volume of left atrium mL
            FETAL
              refer to:  Figure 1 in Pennati et al.'s (2000) "Scaling Approach to Study the Changes Through the Gestationof Human Fetal Cardiac and Circulatory Behaviors"
              concatenation of 2 strings and a character (valv represents valve)
              First string is from the organ
              Second string is to the organ
              Last character represents the components: (r) resistance, (c) capacitance, (l) inductance
                k and b are the non-linear relation term where pressure = k*flowrate**b
        '''
        tempAge='fetal28'
        if isinstance(defaultParameters,str):
            if defaultParameters[:5]=='adult' or defaultParameters[:5]=='fetal':
                tempAge=defaultParameters
        if not(isinstance(defaultParameters,dict)):
            if isinstance(defaultAge,dict):
                defaultParameters=defaultAge.copy()
            else:
                defaultParameters=None
        if not(isinstance(defaultAge,str)):
            defaultAge=tempAge
        elif defaultAge[:5]!='adult' and defaultAge[:5]!='fetal':
            defaultAge=tempAge
        self['heart to solve for inverse']=False
        self['path of folder with case']='.'
        self['filename of stl']='t0'
        self['name of mesh']='t0'
        self['subfolder to save files']=''
        self['subfolder to load files']=''
        self['common communicator']=fenics.Mesh().mpi_comm()
             
        self['LV base mesh surface ID'] = 4
        self['LV endocardial mesh surface ID'] = 2
        self['LV epicardial mesh surface ID'] = 1
        self['RV endocardial mesh surface ID'] = 3
        
        self['LV longitudinal axis superior X']=None # longitudinal axis X towards base
        self['LV longitudinal axis superior Y']=None # longitudinal axis Y
        self['LV longitudinal axis superior Z']=None # longitudinal axis Z
        self['LV lateral axis X']=None # longitudinal axis X rv towards lv
        self['LV lateral axis Y']=None # longitudinal axis Y
        self['LV lateral axis Z']=None # longitudinal axis Z
        self['ratio to clip mesh']=0.95
        self['mesh element size factor']=0.6
        
        self['spring constant for LV pericardial cavity in Pa']=90
        
        self['LV end diastolic volume in mL'] = 3.004778703264237
        self['LV end diastolic pressure in mmHg'] = None
        self['LV end systolic volume in mL'] = 1.4
        self['duration of LV diastole in ms'] = None
        
        self['RV end diastolic volume in mL'] = 3.004778703264237
        self['RV end diastolic pressure in mmHg'] = None
        self['RV end systolic volume in mL'] = 1.4
        
        self['fiber relaxation based on phase during FEA']=False
        self['LV sarcomere length'] = 1.85
        
        self["strain energy density function coefficient in J mL-1"]=100.
        self["strain energy density function exponential coefficient in fiber direction"]=29.9
        self["strain energy density function exponential coefficient in fiber sheetlet direction"]=13.3
        self["strain energy density function exponential coefficient in fiber normal direction"]=13.3
        self["strain energy density function exponential coefficient in cross fiber normal and fiber sheetlet direction"]=26.6
        self["strain energy density function exponential coefficient in cross fiber and fiber sheetlet direction"]=53.2
        self["strain energy density function exponential coefficient in cross fiber and fiber normal direction"]=53.2
        
        self['LV fiber sheetlet angle in degrees']=0.
        self['LV fiber sheetlet width']=0.
        self['LV fiber radial angle in degrees']=0.
        #Active Material
        if defaultAge[:5]=='adult':
            self['duration of one cardiac cycle in ms'] = 800.0 #set for the duration of cardiac cycle value is in ms, for fetal is 400ms. for adult is 800ms
            self['peak intracellular calcium concentration in uM'] = 4.35 #peak intracellular calcium concentration, uM
            self['maximum peak intracellular calcium concentration in uM'] = 4.35 #maximum peak intracellular calcium concentration, uM 
            self['exponential coefficient for relation of peak isometric tension and sarcomere length in um-1'] = 4.75 #governs shape of peak isometric tension-sarcomere length relation, um-1
            for cavity in ["LA","RA","LV","RV"]:
                self['time to maximum '+cavity+' fiber tension in ms'] = 300.5#238 #200.5#170 #132.5 #time to peak tension, ms
            self['LV sarcomere length threshold where no tension develops']= 1.58
            #self['LV sarcomere stretch ratio threshold where no tension develops'] = 1.58/self['LV sarcomere length']#1.58 #sarcomere length at which no active tension develops,um
            self['slope of linear relation of relaxation duration and sarcomere length in ms um-1'] = 1048#1049 #slope of linear relaxation duration-sarcomere length relation, ms um-1
            self['time intercept of linear relation of relaxation duration and sarcomere length in ms'] = -1600#-1429 #time-intercept of linear relaxation duration-sarcomere length relation, ms
            self['LV endocardial fiber angle in degrees']=80.
            self['LV epicardial fiber angle in degrees']=-70.
            self['maximum LV fiber tension in Pa'] = 200.7e3
        else:
            self['duration of one cardiac cycle in ms'] = 400.0 #set for the duration of cardiac cycle value is in ms, for fetal is 400ms. for adult is 800ms
            self['peak intracellular calcium concentration in uM'] = 4.35 #peak intracellular calcium concentration, uM
            self['maximum peak intracellular calcium concentration in uM'] = 4.35 #maximum peak intracellular calcium concentration, uM 
            self['exponential coefficient for relation of peak isometric tension and sarcomere length in um-1'] = 4.75 #governs shape of peak isometric tension-sarcomere length relation, um-1
            for cavity in ["LA","RA","LV","RV"]:
                self['time to maximum '+cavity+' fiber tension in ms'] = 140.5#238 #150.5#170 #132.5 #time to peak tension, ms
            self['LV sarcomere length threshold where no tension develops']= 1.58
            #self['LV sarcomere stretch ratio threshold where no tension develops'] = 1.58/self['LV sarcomere length']#1.58 #sarcomere length at which no active tension develops,um
            self['slope of linear relation of relaxation duration and sarcomere length in ms um-1'] = 1048*0.5#1049 #slope of linear relaxation duration-sarcomere length relation, ms um-1
            self['time intercept of linear relation of relaxation duration and sarcomere length in ms'] = -1600*0.5#0.5*(m_adult*l0+b_adult)-m_fetal #time-intercept of linear relaxation duration-sarcomere length relation, ms
            self['LV endocardial fiber angle in degrees']=60.
            self['LV epicardial fiber angle in degrees']=-60.
            self['maximum LV fiber tension in Pa'] = 60e3
        self.setDefaultWindkessel(defaultAge)
        self.changeDefaultParameters(defaultParameters)
            
    def duplicate(self):
        return heartParameters(defaultParameters=self.copy())
    def __add__(self,addDict):
        temp=self.duplicate()
        tempDict=addDict.copy()
        for key in addDict:
            temp[key]=tempDict[key]
        return temp
    def setHeartrate(self,beatsperminute):
        self['duration of one cardiac cycle in ms'] =60000./beatsperminute
        self['peak intracellular calcium concentration in uM']=self['duration of one cardiac cycle in ms']*0.375+0.5
        self['slope of linear relation of relaxation duration and sarcomere length in ms um-1']=self['duration of one cardiac cycle in ms']/800.*1048.
        self['time intercept of linear relation of relaxation duration and sarcomere length in ms']=self['duration of one cardiac cycle in ms']/800.*-1600.
    def setDefaultWindkessel(self,modelString):
        if modelString=='adult':
            raise Exception('Adult model not supported in this version.')
            self['AV'] = 160.0 #this is value set for left atrium which will be used in below 
        
            #### For Calculating P_LA ######################################## 
            '''
            self['Ees_la'] = 60*1e0; #LCL #this is value for End-systolic elastance for left atrium Pa/m 
            self['A_la'] = 58.05*10; #LCL #Scaling factor for EDPVR Pa
            self['B_la'] = 0.049; #Exponent for EDPVR mL-1
            self['V0_la'] = 0.25#1.0*5; #Volume axis intercept mL 
            self['Tmax_la'] = 125; #Time to end-systole msec
            self['tau_la'] = 25; #Time constant of relaxation msec
            
            self['Cao'] = 0.0065*0.1 #this is value for compliance/capaciator of aortic valve  mLPa
            self['Cven'] = 0.01*2 #this is value for compliance/capaciator of vein,mLPa
            self['Vart0'] = 7.5#66*2;#450; #this is value for intial volume for artery  mL
            self['Vven0'] = 96.5#34*10*2; #this is value for intial volume for vein mL
            self['Rao'] = 10*5500*0.8*16*0.25*0.25*8*2*2*0.5;#10000.0; #this is value for resistance of aortic valve  ,unit is PamsecmL-1
            self['Rven'] = 1000*0.5;#20000.0; #this is value for resistance of vein unit is PamsecmL-1
            self['Rper'] = 100*140000*0.5*0.5*0.5*0.5*4*0.5*0.6*0.65*0.65*2*2#170000.0; #this is value for peripheral arteries resistance unit is PamsecmL-1
            self['Rmv'] = 2500*2*2;#10000.0; #this is value for mitral valve unit is PamsecmL-1
            self['V_ven'] = 92.5#37*10*2*12; #this is value for volume of vein mL
            self['V_art'] = 18.5#74*2.5;#440 #this is value for volume of arteries mL
            self['V_LA'] = 0.3#1.2*20; #this is value for volume of left atrium mL
            '''
            
            #for ngspice
            self['lac'] =1.
            self['artc'] =0.0032
            self['venc'] =0.28

            self['lvartl'] =0.0 *(1000.**2)
            self['artvenl'] =0.0 *(1000.**2)
            self['venlal'] =0.0 *(1000.**2)
            
            self['lvartr'] =500. *(1000.)
            self['artvenr'] =6000. *(1000.)
            self['venlar'] =100. *(1000.)
            
            self['lalavalvr'] =200.0 *(1000.)
            
            self['lalavalvk'] =0.002 *(1000.**2)
            self['lvlvvalvk'] =0.001 *(1000.**2)
            
            self['lalavalvb'] =2.
            self['lvlvvalvb'] =2.
            
            self['lvregurgevalveratio']=-1
            self['aaregurgevalveratio']=-1
            self['rvregurgevalveratio']=-1
            self['pa1regurgevalveratio']=-1
            
            self['reference ngspice circuit filename']='simpleLV'
        elif modelString[:5]=='fetal':
            #resistance r in mmHg ms mL-1
            #resistance k in mmHg ms2 mL-1
            #compliance c in mmHg mL-1
            #inertance l in mmHg ms2 mL-1
            self['Windkessel circuit to use']='full'
            self['scale Windkessel active pressure']=1.
            
            self['aac'] =0.05
            self['ao1c'] =0.08
            self['ao2c'] =0.07
            self['ao3c'] =0.04
            self['ao4c'] =0.05
            self['pa1c'] =0.08
            self['pa2c'] =0.08
            self['lungc'] =0.4
            self['cac'] =0.01
            self['brc'] =0.3
            self['ubc'] =0.85
            self['intec'] =0.25
            self['kidc'] =0.02
            self['hec'] =3.
            self['placc'] =1.5
            self['legc'] =4.
            self['uvc'] = 0.3
            self['svcc'] =1.
            self['ivcc'] =0.6
            self['rac'] =2.
            self['lac'] =1.
            
            self['rvc'] =2.
            self['lvc'] =2.

            self['dal'] =0.006 *(1000.**2)
            self['ao1cal'] =0.08 *(1000.**2)
            self['aaao1l'] =0.002 *(1000.**2)
            self['pa1pa2l'] =0.002 *(1000.**2)
            self['raravalvl'] =0.0016 *(1000.**2)
            self['lalavalvl'] =0.0016 *(1000.**2)
            self['rvrvvalvl'] =0. *(1000.**2)
            self['lvlvvalvl'] =0. *(1000.**2)
            
            self['aaao1r'] =0.12 *(1000.)
            self['ao1ao2r'] =0.4 *(1000.)
            self['ao2ao3r'] =0.04 *(1000.)
            self['ao3ao4r'] =0.06 *(1000.)
            self['pa1pa2r'] =0.07 *(1000.)
            self['pa2lungr'] =13.5 *(1000.)
            self['dar'] =0.01 *(1000.)
            self['ao1car'] =0.3 *(1000.)
            self['cabrr'] =3. *(1000.)
            self['brsvcr'] =8.5 *(1000.)
            self['ao1ubr'] =8. *(1000.)
            self['ubsvcr'] =4.9 *(1000.)
            self['ao3her'] =81. *(1000.)
            self['ao3inter'] =34 *(1000.)
            self['inteher'] =7 *(1000.)
            self['ao3kidr'] =3.5 *(1000.)
            self['kidivcr'] =14 *(1000.)
            self['ao4placr'] =3.9 *(1000.)
            self['placuvr'] =3.4 *(1000.)
            self['ao4legr'] =3.5 *(1000.)
            self['legivcr'] =0.6 *(1000.)
            self['uvher'] =0.5 *(1000.)
            self['heivcr'] =0.16 *(1000.)
            self['dvr'] =1.3 *(1000.)
            self['svcrar'] =0.2 *(1000.)
            self['ivcrar'] =0.12 *(1000.)
            self['lunglar'] =2 *(1000.)
            self['for'] =0. *(1000.)
            self['raravalvr'] =0. *(1000.)
            self['rvrvvalvr'] =0.08 *(1000.)
            self['lvlvvalvr'] =0.08 *(1000.)
            self['lalavalvr'] =0. *(1000.)
            
            self['fok'] =0.4 *(1000.**2)
            self['dak'] =0.009 *(1000.**2)
            self['dvk'] =0.26 *(1000.**2)
            self['raravalvk'] =0.002 *(1000.**2)
            self['rvrvvalvk'] =0.001 *(1000.**2)
            self['lalavalvk'] =0.002 *(1000.**2)
            self['lvlvvalvk'] =0.001 *(1000.**2)
            
            self['fob'] =0.625
            self['dab'] =2.
            self['dvb'] =2.
            self['raravalvb'] =2.
            self['rvrvvalvb'] =2.
            self['lalavalvb'] =2.
            self['lvlvvalvb'] =2.
            
            self['lvregurgevalveratio']=-1
            self['aaregurgevalveratio']=-1
            self['rvregurgevalveratio']=-1
            self['pa1regurgevalveratio']=-1
            
            self.scaleWinkessel({'default':defaultAdjustmentToScaling['r_scale']},compstr='r')
            self.scaleWinkessel({'default':defaultAdjustmentToScaling['c_scale']},compstr='c')
            if len(modelString)>5:
                self.scaleWinkesselwithAge(float(modelString[5:]))
            
            self['reference ngspice circuit filename']='LV'
        self['duration of LV diastole in ms']=None
        self['vla0']=None
        self['vra0']=None
        
        for cavity in ["LA","RA","LV","RV"]:
            self['Windkessel '+cavity+' source file']='./'+cavity.lower()+'ufile.txt'
            self['Windkessel '+cavity+' source pulse function peak pressure in mmHg']=0
            self['Windkessel '+cavity+' source pulse function time at peak pressure in ms']=1
            self['Windkessel '+cavity+' source pulse function pressure pulse width in ms']=1
            for n in range(4):
                self['Windkessel '+cavity+' source fourier function sine amplitude term '+str(n)]=0.
                self['Windkessel '+cavity+' source fourier function sine degree phase term '+str(n)]=0.
        self['Windkessel LA source function']='pulse'
        self['Windkessel RA source function']='pulse'
        self['Windkessel LV source function']='pressure 2D table with tracking of relaxation phase'
        self['Windkessel RV source function']='fourier current'
    def getParameterRelation(self,parameterStringList,return_integer=True):
        '''
        get a list of intergers which
            1: mesh generation variable
            2: FEniCS (FEM) variables
            3: ngspice (Windkessel) calculation variables
        '''
        relation=[]
        for parameterString in parameterStringList:
            if parameterString in parameters_for_mesh:
                if return_integer:
                    relation.append(0)
                else:
                    relation.append('MESH')
            elif parameterString in parameters_for_FEniCS:
                if return_integer:
                    relation.append(1)
                else:
                    relation.append('FENICS')
            else:
                if return_integer:
                    relation.append(2)
                else:
                    relation.append('NGSPICE')
        return relation
    def scaleWinkessel(self,scaleDict,compstr=''):
        '''
        scaleDict={'default':,'aac':,...}
        '''
        if 'default' not in scaleDict:
            scaleDict['default']=1.
        if compstr!='':
            for comp in WindkesselComponents + WindkessellinkComponents:
                if comp+compstr in self:
                    if comp+compstr in scaleDict:
                        self[comp+compstr]*=scaleDict[comp+compstr]
                    elif compstr in scaleDict:
                        self[comp+compstr]*=scaleDict[compstr]
                    else:
                        self[comp+compstr]*=scaleDict['default']
        else:
            self.scaleWinkessel(scaleDict,compstr='r')
            self.scaleWinkessel(scaleDict,compstr='l')
            self.scaleWinkessel(scaleDict,compstr='c')
            self.scaleWinkessel(scaleDict,compstr='k')
            self.scaleWinkessel(scaleDict,compstr='b')
    def scaleWinkesselwithAge(self,ageInWeeks,poweradjustDict=None,compstr=''):
        '''
        set poweradjustDict to enable default
        '''
        if poweradjustDict is None:
            poweradjustDict={}
        if compstr!='':
            for comp in WindkesselComponents + WindkessellinkComponents:
                if comp+compstr in self:
                    if comp+compstr in defaultAgeScalePower:
                        agepower=defaultAgeScalePower[comp+compstr]
                    else:
                        agepower=defaultAgeScalePower['default'+compstr]#1.33
                    if comp+compstr in poweradjustDict:
                        agepower-=defaultAgeScalePower['default'+compstr]-poweradjustDict[comp+compstr]*defaultAgeScalePower['default'+compstr]
                    elif compstr in poweradjustDict:
                        agepower-=defaultAgeScalePower['default'+compstr]-poweradjustDict[compstr]*defaultAgeScalePower['default'+compstr]
                    self[comp+compstr]*=((10.**(0.2508+0.1458*ageInWeeks-0.0016*ageInWeeks**2.))/(10.**(0.2508 + 0.1458*38.-0.0016*38.**2.)))**agepower
        else:
            self.scaleWinkesselwithAge(ageInWeeks,poweradjustDict=poweradjustDict,compstr='r')
            self.scaleWinkesselwithAge(ageInWeeks,poweradjustDict=poweradjustDict,compstr='l')
            self.scaleWinkesselwithAge(ageInWeeks,poweradjustDict=poweradjustDict,compstr='c')
            self.scaleWinkesselwithAge(ageInWeeks,poweradjustDict=poweradjustDict,compstr='k')
            self.scaleWinkesselwithAge(ageInWeeks,poweradjustDict=poweradjustDict,compstr='b')
    def changeDefaultParameters(self,editParameters=None):
        if editParameters is None:
            return
        else:
            for key in self.keys():
                if key in editParameters:
                    self[key]=editParameters[key]
    def searchParameter(self,parameterString,returnTop=3):
        if parameterString in self.keys():
            if returnTop<1:
                return parameterString
            else:
                return [parameterString]
        wordList=parameterString.split(' ')
        searchScore=[]
        keys=self.keys()
        for key in keys:
            searchScore_temp=[0,0]
            text_all=key.split(' ')
            for word in wordList:
                if word in text_all:
                    temp_best_score=1.
                else:
                    temp_best_score=0.
                    for text in text_all:
                        if (word in text):
                            new_best_score=float(len(word))/len(text)
                            if temp_best_score<new_best_score:
                                temp_best_score=new_best_score
                if temp_best_score>0:
                    searchScore_temp[0]+=temp_best_score
                else:
                    searchScore_temp[1]+=1
            searchScore_temp[1]=(len(wordList)-searchScore_temp[1])/len(wordList)
            searchScore.append(np.round(searchScore_temp[0])+searchScore_temp[1])
        if returnTop<1:
            return list(keys)[np.argmax(searchScore)]
        elif returnTop<=-1:
            return list(keys)[np.argsort(searchScore)[::-1]]
        else:
            return list(keys)[np.argsort(searchScore)[::-1][:int(returnTop)]]            
    def setParameter(self,parameterString,value):
        key=self.searchParameter(parameterString,returnTop=0)
        self[key]=value
        logger.info('set '+key+" = "+repr(value))
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
    def readParameters(self,file):
        parameters={}
        with open(file, "r") as f:
            lines=f.readlines()
        for line in lines:
            if line[0]=='#':
                continue
            temp=line.split(sep=' , ')
            parameters[temp[0]]=self.decodeStringValue(temp[1])
        for key in parameters:
            self.setParameter(key,parameters[key])
    def writeParameters(self,file):
        folder,filename = os.path.split(file)
        os.makedirs(folder,exist_ok=True)
        with open(file, "w") as f:
            for key, val in self.items():
                f.write("{0:s} , {1:s}\n".format(key,str(val)))
    def getFetalPopulation_LVEDP(self,ageInWeeks):
        #from P Johnsona, D J Maxwellb, M J Tynanb, L D Allanc, "Intracardiac pressures in the human fetus"
        return 0.7018*ageInWeeks-7.62928682
    #other possible volume reference: doi=10.3760/cma.j.issn.1004-4477.2016.07.005
    def getFetalPopulation_min_LAwidth(self,ageInWeeks):
        #Yan Kai Mao, Bo Wen Zhao, Feng Hua Zheng, Bei Wang, Xiao Hui Peng, Ran Chen, Mei Pan, "Z-scores for fetal left atrial size and left atrium–descending aorta distance in fetuses with isolated total anomalous pulmonary venous connection"
        return -0.11399+0.04125*ageInWeeks
    def getFetalPopulation_min_LAlength(self,ageInWeeks):
        #Yan Kai Mao, Bo Wen Zhao, Feng Hua Zheng, Bei Wang, Xiao Hui Peng, Ran Chen, Mei Pan, "Z-scores for fetal left atrial size and left atrium–descending aorta distance in fetuses with isolated total anomalous pulmonary venous connection"
        return -0.13371+0.05186*ageInWeeks
    def getFetalPopulation_min_LAarea(self,ageInWeeks):
        #Yan Kai Mao, Bo Wen Zhao, Feng Hua Zheng, Bei Wang, Xiao Hui Peng, Ran Chen, Mei Pan, "Z-scores for fetal left atrial size and left atrium–descending aorta distance in fetuses with isolated total anomalous pulmonary venous connection"
        return -1.4234+0.095489*ageInWeeks
    
    def getFetalPopulation_LAwidth(self,ageInWeeks):
        #Julie  Tan,  MD, Norman H.  Silverman, MD, DSc(Med),  Julien  I. E.  Hoffman, MD, Maria  Villegas,  and  Klaus  G.  Schmidt, MD, "Cardiac Dimensions Determined by Cross-Sectional Echocardiography in  the Normal Human Fetus  from 18  Weeks  to  Term "
        return -1.246+0.1305*ageInWeeks-0.001563*ageInWeeks**2.
    def getFetalPopulation_LAlength(self,ageInWeeks):
        #Julie  Tan,  MD, Norman H.  Silverman, MD, DSc(Med),  Julien  I. E.  Hoffman, MD, Maria  Villegas,  and  Klaus  G.  Schmidt, MD, "Cardiac Dimensions Determined by Cross-Sectional Echocardiography in  the Normal Human Fetus  from 18  Weeks  to  Term "
        return -0.6508+0.0873*ageInWeeks-0.000674*ageInWeeks**2.
    def getFetalPopulation_LAwidth2(self,ageInWeeks):
        #Julie  Tan,  MD, Norman H.  Silverman, MD, DSc(Med),  Julien  I. E.  Hoffman, MD, Maria  Villegas,  and  Klaus  G.  Schmidt, MD, "Cardiac Dimensions Determined by Cross-Sectional Echocardiography in  the Normal Human Fetus  from 18  Weeks  to  Term "
        return -0.3422+0.0498*ageInWeeks
    def getFetalPopulation_LAlength2(self,ageInWeeks):
        #Julie  Tan,  MD, Norman H.  Silverman, MD, DSc(Med),  Julien  I. E.  Hoffman, MD, Maria  Villegas,  and  Klaus  G.  Schmidt, MD, "Cardiac Dimensions Determined by Cross-Sectional Echocardiography in  the Normal Human Fetus  from 18  Weeks  to  Term "
        return -0.6408+0.0707*ageInWeeks
    def getFetalPopulation_LV_SV(self,ageInWeeks):
        #F. S. Molina, C. Faro, A. Sotiriadis, T. Dagklis, K. H. Nicolaides, "Heart stroke volume and cardiac output by four-dimensional ultrasound in normal fetuses"
        #given coefficient seems to be wrong
        #return np.exp(-12.662+0.136*ageInWeeks-(4.715*10.**-4)*ageInWeeks**2.-(5.597*10.**-7)*ageInWeeks**3.
        return 2.27637837-0.360074530*ageInWeeks+(1.68491289*10.**-2)*ageInWeeks**2.-(1.88355120*10.**-4)*ageInWeeks**3.
    def getFetalPopulation_LAEDV(self,ageInWeeks):
        LASV=self.getFetalPopulation_LV_SV(ageInWeeks)
        #Panupong Jiamsripong, Tadaaki Honda, Christina S. Reuss, R. Todd Hurst, Hari P. Chaliki, Diane E. Grill, Stephen L. Schneck, Rochelle Tyler, Bijoy K. Khandheria, Steven J. Lester, "Three methods for evaluation of left atrial volume"
        #method 3
        #ED_LA_volume=-LASV+0.523*((self.getFetalPopulation_LAlength(ageInWeeks)+self.getFetalPopulation_LAlength2(ageInWeeks))/2.)*(self.getFetalPopulation_LAwidth(ageInWeeks))*(self.getFetalPopulation_LAwidth2(ageInWeeks))
        #method 1 with A1=A2
        ED_LA_volume=2.*0.85*(self.getFetalPopulation_min_LAarea(ageInWeeks)**2.)/(self.getFetalPopulation_min_LAlength(ageInWeeks))
        return ED_LA_volume
    def getFetalPopulation_LA_unloadedVolume(self,ageInWeeks):
        LAEDV=self.getFetalPopulation_LAEDV(ageInWeeks)
        LAEDP=self.getFetalPopulation_LVEDP(ageInWeeks)
        LA_Capacitance=heartParameters(defaultAge='fetal'+str(ageInWeeks))['lac']
        return LAEDV-LAEDP*LA_Capacitance