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
'''
_version='3.5.1'
import logging
logger = logging.getLogger(__name__)

parameters_for_FEniCS=['Kspring_constant','Tact_constant','T0_LV','ESV_LV','lr','BCL','Ca0','Ca0max','B','t0','l0','m','b']                 
parameters_for_mesh=['topid','endoid','epiid','Laxis_X','Laxis_Y','Laxis_Z','clip_ratio','endo_angle','epi_angle','fiberSheetletAngle','fiberSheetletWidth','radialFiberAngle','EDV_LV','EDP_LV']
WindkesselComponents=['lv','la','rv','ra','aa','ao1','ao2','ao3','ao4','br','ca','ub','he','inte','ivc','kid','leg','lung','pa1','pa2','plac','svc','uv']
WindkessellinkComponents=['aaao1','ao1ao2','ao2ao3','ao3ao4','pa1pa2','pa2lung','da','ao1ca','cabr','brsvc','ao1ub','ubsvc','ao3he','ao3inte','intehe','ao3kid','kidivc','ao4plac','placuv','ao4leg','legivc','uvhe','heivc','dv','svcra','ivcra','lungla','fo','raravalv','lalavalv','lvlvvalv','rvrvvalv']
defaultAgeScalePower={'defaultr':-1.,'pa2lungr':-1.2,'lunglar':-1.2,'cabrr':-1.1,'brsvcr':-1.1,'dvr':-0.55,
                      'defaultl':-0.33,
                      'defaultc':1.33,'brc':1.471,'lungc':1.6,'ra':0.5,'la':0.5,
                      'defaultk':0.,'fok':-0.6,'dak': -2.5,'dvk':-0.88,'raravalvk':-1.33,'rvrvvalvk':-1.33,'lalavalvk':-1.33,'lvlvvalvk':-1.33,
                      'defaultb':0.}
defaultAdjustmentToScaling={'r_scale':1.21,'c_scale':0.27,'C_adj':1./1.33,'R_adj':-0.58/-1.}



class heartParameters(dict):
    def __new__(cls, *args, **kwargs):
        return super().__new__(cls)
    def __init__(self,defaultParameters=None,defaultAge=None):
        '''
        MESH parameters:
            'topid'           #id for the base it is on top of geometry
            'endoid'          #id for the inner surface it is on top of geometry
            'epiid'           #id for the outer surface it is on top of geometry
            'Laxis_X'         # longitudinal axis X
            'Laxis_Y'         # longitudinal axis Y
            'Laxis_Z'         # longitudinal axis Z
            'clip_ratio'      #ratio to clip the stl 1 is not to clip 0 is to clip all
       FENICS parameters:
            'Kspring_constant'     #Kspring constant which model as pericardial cavity , unit in Pa
            'Tact_constant'        #Tact_constant we dont use it, unit is in Pa 
            'T0_LV'                #active tension forces unit in Pa
            'EDV_LV'               #End-diastolic volume for left ventricle
            'EDP_LV'               #End-diastolic pressure for left ventricle
            'ESV_LV'               #End-systolic volume for left ventricle
            'lr'                   # relaxed sarcomere length
            'BCL'                  # base cycle length (time) in milliseconds 
            'Ca0'                  # peak intracellular calcium concentration
            'Ca0max'               #maximum intracellular calcium concentration for calculating the length-dependent calcium sensitivity variable
            'B'                    # exponential constant for length-dependent calcium sensitivity variable
            't0'                   # time to peak tention
            'l0'                   # length of sarcomere for which below this value, has no tension developed
            'm'                    # slope for linear variation of relaxation time with sarcomere length
            'b'                    # y-intercept for linear variation of relaxation time with sarcomere length
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
        self['topid'] = 4 #id for the base it is on top of geometry
        self['endoid'] = 2 #id for the inner surface it is on top of geometry
        self['epiid'] = 1 #id for the outer surface it is on top of geometry
        
        self['Laxis_X']=None # longitudinal axis X
        self['Laxis_Y']=None # longitudinal axis Y
        self['Laxis_Z']=None # longitudinal axis Z
        self['clip_ratio']=0.95
        
        self['Kspring_constant']=90 #Kspring constant which model as pericardial cavity , unit in Pa
        self['Tact_constant'] = 1e5 #Tact_constant we dont use it, unit is in Pa 
        
        self['EDV_LV'] = 3.004778703264237
        self['EDP_LV'] = None
        self['ESV_LV'] = 1.4
        self['ES_time'] = None
        
        
        self['lr'] = 1.85
        
        self["StrainEnergyDensityFunction_Cff"]=29.9
        self["StrainEnergyDensityFunction_Css"]=13.3
        self["StrainEnergyDensityFunction_Cnn"]=13.3
        self["StrainEnergyDensityFunction_Cns"]=26.6
        self["StrainEnergyDensityFunction_Cfs"]=53.2
        self["StrainEnergyDensityFunction_Cfn"]=53.2
        
        self['fiberSheetletAngle']=0.
        self['fiberSheetletWidth']=0.
        self['radialFiberAngle']=0.
        #Active Material
        if defaultAge[:5]=='adult':
            self['BCL'] = 800.0 #set for the duration of cardiac cycle value is in ms, for fetal is 400ms. for adult is 800ms
            self['Ca0'] = 4.35 #peak intracellular calcium concentration, µM
            self['Ca0max'] = 4.35 #maximum peak intracellular calcium concentration, µM 
            self['B'] = 4.75 #governs shape of peak isometric tension-sarcomere length relation, µm−1
            self['t0'] = 300.5#238 #200.5#170 #132.5 #time to peak tension, ms
            self['l0'] = 1.58#1.58 #sarcomere length at which no active tension develops,µm
            self['m'] = 1048#1049 #slope of linear relaxation duration-sarcomere length relation, ms µm−1
            self['b'] = -1600#-1429 #time-intercept of linear relaxation duration-sarcomere length relation, ms
            self['endo_angle']=80.
            self['epi_angle']=-70.
            self['T0_LV'] = 200.7e3
        else:
            self['BCL'] = 400.0 #set for the duration of cardiac cycle value is in ms, for fetal is 400ms. for adult is 800ms
            self['Ca0'] = 4.35 #peak intracellular calcium concentration, µM
            self['Ca0max'] = 4.35 #maximum peak intracellular calcium concentration, µM 
            self['B'] = 4.75 #governs shape of peak isometric tension-sarcomere length relation, µm−1
            self['t0'] = 150.5#238 #150.5#170 #132.5 #time to peak tension, ms
            self['l0'] = 1.58#1.58 #sarcomere length at which no active tension develops,µm
            self['m'] = 1048*0.5#1049 #slope of linear relaxation duration-sarcomere length relation, ms µm−1
            self['b'] = -1600*0.5#0.5*(m_adult*l0+b_adult)-m_fetal #time-intercept of linear relaxation duration-sarcomere length relation, ms
            self['endo_angle']=60.
            self['epi_angle']=-60.
            self['T0_LV'] = 60e3
        self.setDefaultWindkessel(defaultAge)
        self.changeDefaultParameters(defaultParameters)
    def setHeartrate(self,beatsperminute):
        self['BCL'] =60000./beatsperminute
        self['t0']=self['BCL']*0.375+0.5
        self['m']=self['BCL']/800.*1048.
        self['b']=self['BCL']/800.*-1600.
    def setDefaultWindkessel(self,modelString):
        if modelString=='adult':
            self['AV'] = 160.0 #this is value set for left atrium which will be used in below 
        
            #### For Calculating P_LA ######################################## 
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
        elif modelString[:5]=='fetal':
            self['Windkessel_scale_T0_LV']=1.
            self['Aortic_stenosis']=1.
            
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

            self['dal'] =0.006
            self['ao1cal'] =0.08
            self['aaao1l'] =0.002
            self['pa1pa2l'] =0.002
            self['raravalvl'] =0.0016
            self['lalavalvl'] =0.0016
            self['rvrvvalvl'] =0.
            self['lvlvvalvl'] =0.
            
            self['aaao1r'] =0.12
            self['ao1ao2r'] =0.4
            self['ao2ao3r'] =0.04
            self['ao3ao4r'] =0.06
            self['pa1pa2r'] =0.07
            self['pa2lungr'] =13.5
            self['dar'] =0.01
            self['ao1car'] =0.3
            self['cabrr'] =3.
            self['brsvcr'] =8.5
            self['ao1ubr'] =8.
            self['ubsvcr'] =4.9
            self['ao3her'] =81.
            self['ao3inter'] =34
            self['inteher'] =7
            self['ao3kidr'] =3.5
            self['kidivcr'] =14
            self['ao4placr'] =3.9
            self['placuvr'] =3.4
            self['ao4legr'] =3.5
            self['legivcr'] =0.6
            self['uvher'] =0.5
            self['heivcr'] =0.16
            self['dvr'] =1.3
            self['svcrar'] =0.2
            self['ivcrar'] =0.12
            self['lunglar'] =2
            self['for'] =0.
            self['raravalvr'] =0.
            self['rvrvvalvr'] =0.08
            self['lvlvvalvr'] =0.08
            self['lalavalvr'] =0.
            
            self['fok'] =0.4
            self['dak'] =0.009
            self['dvk'] =0.26
            self['raravalvk'] =0.002
            self['rvrvvalvk'] =0.001
            self['lalavalvk'] =0.002
            self['lvlvvalvk'] =0.001
            
            self['fob'] =0.625
            self['dab'] =2.
            self['dvb'] =2.
            self['raravalvb'] =2.
            self['rvrvvalvb'] =2.
            self['lalavalvb'] =2.
            self['lvlvvalvb'] =2.
            
            self['lvregurger']=-1
            self['lvregurgevalveratio']=-1
            self['rvregurger']=-1
            self['rvregurgevalveratio']=-1
            self['switchvalve']=0
            
            self['ES_time']=None
            self['vla0']=None
            self['vra0']=None
            
            self.scaleWinkessel({'default':defaultAdjustmentToScaling['r_scale']},compstr='r')
            self.scaleWinkessel({'default':defaultAdjustmentToScaling['c_scale']},compstr='c')
            if len(modelString)>5:
                self.scaleWinkesselwithAge(float(modelString[5:]))
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
            poweradjustDict={'r':defaultAdjustmentToScaling['R_adj'],'c':defaultAdjustmentToScaling['C_adj']}
        if compstr!='':
            for comp in WindkesselComponents + WindkessellinkComponents:
                if comp+compstr in self:
                    if comp+compstr in defaultAgeScalePower:
                        agepower=defaultAgeScalePower[comp+compstr]
                    else:
                        agepower=defaultAgeScalePower['default'+compstr]
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
