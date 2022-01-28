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
File: heartModelLinker.py
Description: object to control external heart modules
consolidate units to mL and mmHg
History:
    Date    Programmer SAR# - Description
    ---------- ---------- ----------------------------
  Author: w.x.chan@gmail.com         05NOV2021           - Created
  Author: w.x.chan@gmail.com         05NOV2021           - v4.0.0
'''
_version='4.0.0'
import logging
logger = logging.getLogger('heartModelLinker v'+_version)
import sys
import numpy as np
import copy
from types import ModuleType
try:
    import dolfin as fenics
except:
    import fenics as fenics
from . import lcleeHeart

try:
    import heArt
except Exception as e:
    logger.warning('Unable to load Module "heArt".')
    logger.warning(repr(e))

'''
Set args and kwargs to pass to heart and generatemesh for unknown external modules
'''
convert_Pa_to_mmHg=0.0075
inbuilt_generateMesh_args_search={}
inbuilt_changeFiberAngles_args_search={}
inbuilt_createHeartModel_args_search={}
inbulit_function_names={}
inbuilt_convertParameters={}
def toNumpyArray(*args):
    return np.array(args)
def toList(*args):
    return list(args)
def add(arg1,arg2):
    return arg1+arg2
if 'lcleeHeart' in locals():
    inbulit_function_names['lcleeHeart']={'generateMesh':lcleeHeart.generateMesh,
                                        'changeFiberAngles':lcleeHeart.changeFiberAngles,
                                        'createHeartModel':lcleeHeart.heart}
    inbuilt_generateMesh_args_search['lcleeHeart']={'args':['casePath',
                                                          'stlname',
                                                          'Laxis'],
                                                  'kwargs':['endo_angle',
                                                            'epi_angle',
                                                            'fiberSheetletAngle',
                                                            'fiberSheetletWidth',
                                                            'radialFiberAngle',
                                                            'fiberLength',
                                                            'clipratio',
                                                            'meshsize',
                                                            'meshname',
                                                            'saveSubFolder'],
                                                  'convertParameters':{'casePath':'path of folder with case',
                                                                       'stlname':'filename of stl',
                                                                       'Laxis':[toNumpyArray,'LV longitudinal axis superior X','LV longitudinal axis superior Y','LV longitudinal axis superior Z'],
                                                                       'endo_angle':'LV endocardial fiber angle in degrees',
                                                                       'epi_angle':'LV epicardial fiber angle in degrees',
                                                                       'fiberSheetletAngle':'LV fiber sheetlet angle in degrees',
                                                                       'fiberSheetletWidth':'LV fiber sheetlet width',
                                                                       'radialFiberAngle':'LV fiber radial angle in degrees',
                                                                       'fiberLength':'LV sarcomere length',
                                                                       'clipratio':'ratio to clip mesh',
                                                                       'meshsize':'mesh element size factor',
                                                                       'meshname':'name of mesh',
                                                                       'saveSubFolder':'subfolder to save files'}}
    inbuilt_changeFiberAngles_args_search['lcleeHeart']={'args':['casePath',
                                                                 'meshname',
                                                                 'endo_angle',
                                                                 'epi_angle'],
                                                'kwargs':['fiberSheetletAngle',
                                                          'fiberSheetletWidth',
                                                          'radialFiberAngle',
                                                          'fiberLength',
                                                          'saveSubFolder',
                                                          'loadSubFolder'],
                                                  'convertParameters':{'casePath':'path of folder with case',
                                                                       'meshname':'name of mesh',
                                                                       'endo_angle':'LV endocardial fiber angle in degrees',
                                                                       'epi_angle':'LV epicardial fiber angle in degrees',
                                                                       'fiberSheetletAngle':'LV fiber sheetlet angle in degrees',
                                                                       'fiberSheetletWidth':'LV fiber sheetlet width',
                                                                       'radialFiberAngle':'LV fiber radial angle in degrees',
                                                                       'fiberLength':'LV sarcomere length',
                                                                       'clipratio':'ratio to clip mesh',
                                                                       'meshsize':'mesh element size factor',
                                                                       'saveSubFolder':'subfolder to save files',
                                                                       'loadSubFolder':'subfolder to load files'}}
    inbuilt_createHeartModel_args_search['lcleeHeart']={'args':['casename',
                                                                'meshname',
                                                                'runParameters'],
                                                'kwargs':['inverseHeart',
                                                          'trackphase'],
                                                  'convertParameters':{'casename':'path of folder with case',
                                                                       'meshname':'name of mesh',
                                                                       'trackphase':'fiber relaxation based on phase during FEA',
                                                                       'inverseHeart':'heart to solve for inverse',
                                                                       'topid':'LV base mesh surface ID',
                                                                       'endoid':'LV endocardial mesh surface ID',
                                                                       'epiid':'LV epicardial mesh surface ID',
                                                                       'BCL':'duration of one cardiac cycle in ms',
                                                                       'Kspring_constant':'spring constant for LV pericardial cavity in Pa',
                                                                       'lr':'LV sarcomere length',
                                                                       'l0':'LV sarcomere length threshold where no tension develops',#[np.multiply,'LV sarcomere length','LV sarcomere stretch ratio threshold where no tension develops'],
                                                                       't0':'time to maximum LV fiber tension in ms',
                                                                       'Ca0':'peak intracellular calcium concentration in uM',
                                                                       'Ca0max':'maximum peak intracellular calcium concentration in uM',
                                                                       'B':'exponential coefficient for relation of peak isometric tension and sarcomere length in um-1',
                                                                       'm':'slope of linear relation of relaxation duration and sarcomere length in ms um-1',
                                                                       'b':'time intercept of linear relation of relaxation duration and sarcomere length in ms',
                                                                       'T0_LV':'maximum LV fiber tension in Pa',
                                                                       "StrainEnergyDensityFunction_Coef":"strain energy density function coefficient in J mL-1",
                                                                       "StrainEnergyDensityFunction_Cff":"strain energy density function exponential coefficient in fiber direction",
                                                                       "StrainEnergyDensityFunction_Css":"strain energy density function exponential coefficient in fiber sheetlet direction",
                                                                       "StrainEnergyDensityFunction_Cnn":"strain energy density function exponential coefficient in fiber normal direction",
                                                                       "StrainEnergyDensityFunction_Cns":"strain energy density function exponential coefficient in cross fiber normal and fiber sheetlet direction",
                                                                       "StrainEnergyDensityFunction_Cfs":"strain energy density function exponential coefficient in cross fiber and fiber sheetlet direction",
                                                                       "StrainEnergyDensityFunction_Cfn":"strain energy density function exponential coefficient in cross fiber and fiber normal direction",
                                                                       'EDV_LV':'LV end diastolic volume in mL',
                                                                       'EDP_LV':'LV end diastolic pressure in mmHg',
                                                                       'ESV_LV':'LV end systolic volume in mL',
                                                                       'ES_time':'duration of LV diastole in ms'}}
    inbuilt_createHeartModel_args_search['lcleeHeart']['convertParameters']['runParameters']=inbuilt_createHeartModel_args_search['lcleeHeart']['convertParameters'].copy()
if 'heArt' in locals():
    inbulit_function_names['heArt.BiV']={'generateMesh':lcleeHeart.generateMesh,
                                        'changeFiberAngles':lcleeHeart.changeFiberAngles,
                                        'createHeartModel':heArt.MEmodel}
    inbuilt_generateMesh_args_search['heArt.BiV']={'args':['casePath',
                                                          'stlname',
                                                          'Laxis'],
                                                  'kwargs':['endo_angle',
                                                            'epi_angle',
                                                            'fiberSheetletAngle',
                                                            'fiberSheetletWidth',
                                                            'radialFiberAngle',
                                                            'fiberLength',
                                                            'clipratio',
                                                            'meshsize',
                                                            'meshname',
                                                            'saveSubFolder',
                                                            'RV_vector'],
                                                  'convertParameters':{'casePath':'path of folder with case',
                                                                       'stlname':'filename of stl',
                                                                       'Laxis':[toNumpyArray,'LV longitudinal axis superior X','LV longitudinal axis superior Y','LV longitudinal axis superior Z'],
                                                                       'endo_angle':'LV endocardial fiber angle in degrees',
                                                                       'epi_angle':'LV epicardial fiber angle in degrees',
                                                                       'fiberSheetletAngle':'LV fiber sheetlet angle in degrees',
                                                                       'fiberSheetletWidth':'LV fiber sheetlet width',
                                                                       'radialFiberAngle':'LV fiber radial angle in degrees',
                                                                       'fiberLength':'LV sarcomere length',
                                                                       'clipratio':'ratio to clip mesh',
                                                                       'meshsize':'mesh element size factor',
                                                                       'meshname':'name of mesh',
                                                                       'saveSubFolder':'subfolder to save files',
                                                                       'RV_vector':[toNumpyArray,'LV lateral axis X','LV lateral axis Y','LV lateral axis Z']}}
    inbuilt_changeFiberAngles_args_search['heArt.BiV']={'args':['casePath',
                                                                 'meshname',
                                                                 'endo_angle',
                                                                 'epi_angle'],
                                                'kwargs':['fiberSheetletAngle',
                                                          'fiberSheetletWidth',
                                                          'radialFiberAngle',
                                                          'fiberLength',
                                                          'saveSubFolder',
                                                          'loadSubFolder',
                                                          'RV_vector'],
                                                  'convertParameters':{'casePath':'path of folder with case',
                                                                       'meshname':'name of mesh',
                                                                       'endo_angle':'LV endocardial fiber angle in degrees',
                                                                       'epi_angle':'LV epicardial fiber angle in degrees',
                                                                       'fiberSheetletAngle':'LV fiber sheetlet angle in degrees',
                                                                       'fiberSheetletWidth':'LV fiber sheetlet width',
                                                                       'radialFiberAngle':'LV fiber radial angle in degrees',
                                                                       'fiberLength':'LV sarcomere length',
                                                                       'clipratio':'ratio to clip mesh',
                                                                       'meshsize':'mesh element size factor',
                                                                       'saveSubFolder':'subfolder to save files',
                                                                       'loadSubFolder':'subfolder to load files',
                                                                       'RV_vector':[toNumpyArray,'LV lateral axis X','LV lateral axis Y','LV lateral axis Z']}}
    inbuilt_createHeartModel_args_search['heArt.BiV']={'args':['params',
                                                               'SimDet'],
                                                'kwargs':[],
                                                  'convertParameters':{'SimDet':{"HeartBeatLength": 'duration of one cardiac cycle in ms',
                                                                                  "dt": 1.0,
                                                                                  "writeStep": 10.0,
                                                                                  "GiccioneParams" : {"ParamsSpecified" : True, 
                                                                                             		  "Passive model": {"Name": "Guccione"},
                                                                                             		  "Passive params": {"Cparam": [fenics.Constant,"strain energy density function coefficient in J mL-1"], 
                                                                                        		                     "bff"  : [fenics.Constant,"strain energy density function exponential coefficient in fiber direction"],
                                                                                        		                     "bfx"  : [fenics.Constant,"strain energy density function exponential coefficient in fiber sheetlet direction"],
                                                                                                          		     "bxx"  : [fenics.Constant,"strain energy density function exponential coefficient in cross fiber normal and fiber sheetlet direction"]},
                                                                                             		  "Active model": {"Name": "Guccione"},
                                                                                             		  "Active params": {'m':'slope of linear relation of relaxation duration and sarcomere length in ms um-1',
                                                                                                                        'b':'time intercept of linear relation of relaxation duration and sarcomere length in ms',
                                                                                                                        "B" : 'maximum peak intracellular calcium concentration in uM', 
                                                                                                                        "t0" : 'time to maximum LV fiber tension in ms', 
                                                                                                                        "l0" : 'LV sarcomere length threshold where no tension develops',
                                                                                                                        "Tmax" : [fenics.Constant,'maximum LV fiber tension in Pa'],
                                                                                                                        "Ca0" : 'peak intracellular calcium concentration in uM',
                                                                                                                        "Ca0max" : 'maximum peak intracellular calcium concentration in uM',
                                                                                                                        "lr" : 'LV sarcomere length'},
                                                                                             		  "HomogenousActivation": True,
                                                                                             		  "deg" : 4, 
                                                                                             		  "Kappa": 1e5,
                                                                                             		  "incompressible" : True}, 
                                                                                  "nLoadSteps": 10,
                                                                                  "DTI_EP": False,
                                                                                  "DTI_ME": False,
                                                                                  "d_iso": 1.5*0.01, 
                                                                                  "d_ani_factor": 4.0, 
                                                                                  "ploc": [[-0.083, 5.6,-1.16, 2.0, 1.0]],
                                                                                  "pacing_timing": [[4.0, 20.0]],
                                                                                  "Ischemia": False,
                                                                             	  "isLV" : False,
                                                                                  "topid" : 'LV base mesh surface ID',
                                                                                  "LVendoid" : 'LV endocardial mesh surface ID',
                                                                                  "RVendoid" : 'RV endocardial mesh surface ID',
                                                                                  "epiid" : 'LV epicardial mesh surface ID',
                                                                                  "abs_tol": 1e-9,
                                                                                  "rel_tol": 1e-9,}}}
    inbuilt_createHeartModel_args_search['heArt.BiV']['convertParameters']['params']={"directory" : [add,'path of folder with case','/'],
                                                                                 "casename" : 'name of mesh', 
                                                                                 "fibre_quad_degree" : 4, 
                                                                                 "outputfolder" : 'subfolder to save files',
                                                                                 "foldername" : 'subfolder to load files',
                                                                                 "state_obj": [heArt.State_Variables,'common communicator', inbuilt_createHeartModel_args_search['heArt.BiV']['convertParameters']['SimDet']],
                                                                                 "common_communicator": 'common communicator',
                                                                                 "MEmesh": [fenics.Mesh],
                                                                                 "isLV": False}
def decode_convertParameters(runlist,refParams):
    if isinstance(runlist,list):
        if callable(runlist[0]):
            args=[]
            kwargs={}
            for n in range(1,len(runlist)):
                if isinstance(runlist[n],dict):
                    temp_kwargs={}
                    for key in runlist[n]:
                        temp_kwargs[key]=decode_convertParameters(runlist[n][key],refParams)
                    if 'kwargs' in runlist[n]:
                        if runlist[n]['kwargs'] is True:
                            kwargs=temp_kwargs
                        else:
                            args.append(temp_kwargs)
                    else:
                        args.append(temp_kwargs)
                else:
                    args.append(decode_convertParameters(runlist[n],refParams))
            return runlist[0](*args,**kwargs)    
        else:
            args=[]
            for n in range(len(runlist)):
                args.append(decode_convertParameters(runlist[n],refParams))
            return args
    elif isinstance(runlist,dict):
        args={}
        for key in runlist:
            args[key]=decode_convertParameters(runlist[key],refParams)
        return args
    elif isinstance(runlist,str):
        if runlist[:4]=='self':
            if len(runlist)==4:
                return refParams
            else:
                var=runlist.split('.')
        elif runlist in refParams:
            return refParams[runlist]
        else:
            return runlist
    else:
        return runlist
def convertParameters(heartModel,params,toConvertParameters,refParams):
    '''
    if inbuilt_convertParameters[key] is list, run inbuilt_convertParameters[key][0](inbuilt_convertParameters[key][1:])
    '''
    new_params={}
    for key in toConvertParameters:
        if key in params:
            if toConvertParameters[key]=='self':
                temp_keys=list(toConvertParameters.keys())
                for find_self_key in temp_keys:
                    if toConvertParameters[find_self_key]=='self':
                        temp_keys.remove(find_self_key)
                new_params[key]=convertParameters(heartModel,temp_keys,toConvertParameters,refParams)
            else:
                new_params[key]=decode_convertParameters(toConvertParameters[key],refParams)
    return new_params
def checkWarnings(funcName,heartModel,params):
    if heartModel=='lcleeHeart':
        if funcName=='createHeartModel':
            if params['fiber relaxation based on phase during FEA']==True:
                logger.warning("Known issue of non-convergence during relaxation")
def functionTemplete(funcName,heartModel,params):
    runFunction=inbulit_function_names[heartModel][funcName]
    runFunction_args_search=globals()['inbuilt_'+funcName+'_args_search'][heartModel]
    runFunction_params_args=convertParameters(heartModel,runFunction_args_search['args'],globals()['inbuilt_'+funcName+'_args_search'][heartModel]['convertParameters'],params)
    runFunction_params_kwargs=convertParameters(heartModel,runFunction_args_search['kwargs'],globals()['inbuilt_'+funcName+'_args_search'][heartModel]['convertParameters'],params)
    logger.debug("Runing "+funcName+" from "+heartModel)
    logger.debug("  Args"+repr(runFunction_params_args))
    logger.debug("  Kwargs"+repr(runFunction_params_kwargs))
    checkWarnings(funcName,heartModel,params)
    pass_args=[]
    pass_kwargs={}
    for key in runFunction_args_search['args']:
        if key in runFunction_params_args:
            pass_args.append(runFunction_params_args[key])
        else:
            raise Exception(repr(key)+' not found in input.')
    for key in runFunction_args_search['kwargs']:
        if key in runFunction_params_kwargs:
            pass_kwargs[key]=runFunction_params_kwargs[key]
    return runFunction(*pass_args,**pass_kwargs)
def generateMesh(heartModel,params):
    functionTemplete('generateMesh',heartModel,params)
    if heartModel=='lcleeHeart':
        return {'LV base mesh surface ID':4,'LV endocardial mesh surface ID':2,'LV epicardial mesh surface ID':1}
    else:
        return {}
def changeFiberAngles(heartModel,params):
    functionTemplete('changeFiberAngles',heartModel,params)
    return {}

class HeartModel:
    def __new__(cls, *args, **kwargs):
        return super().__new__(cls)
    def __init__(self,heartModel,params):
        self.heartModel=heartModel
        self.LV=True
        self.RV=False
        self.LA=False
        self.RA=False
        self.modelSpecificRefference={}
        if heartModel=='heArt.BiV':
            self.RV=True
        elif self.heartModel=='lcleeHeart':
            self.comm=self.heartModelObj.comm
            self.mesh=self.heartModelObj.mesh
        self.heartModelObj=functionTemplete('createHeartModel',heartModel,params)
        if heartModel=='heArt.BiV':
            self.modelSpecificRefference['solver']=self.heartModelObj.Solver()
            self.modelSpecificRefference['state_obj']=self.heartModelObj.state_obj
    def generateMesh(self,params):
        generateMesh(self.heartModel,params)
    def changeFiberAngles(self,params):
        changeFiberAngles(self.heartModel,params)
    def get_CavityNames(self):
        cavity=[]
        if self.LV:
            cavity.append("LV")
        if self.RV:
            cavity.append("RV")
        if self.LA:
            cavity.append("LA")
        if self.RA:
            cavity.append("RA")
        return cavity
    def get_LVVolume(self):
        if self.heartModel=='lcleeHeart':
            return self.heartModelObj.uflforms.cavityvol()
        elif self.heartModel=='heArt.BiV':
            return self.heartModelObj.uflforms.LVcavityvol()
        return
    def get_RVVolume(self):
        if self.heartModel=='heArt.BiV':
            return self.heartModelObj.uflforms.RVcavityvol()
        return
    def get_LAVolume(self):
        return
    def get_RAVolume(self):
        return
    def get_CavityVolume(self):
        cavity=self.get_CavityNames()
        result=[]
        for chamber in cavity:
            result.append(getattr(self,'get_'+chamber+'Volume')())
        return np.array(result)
    def get_LVPressure(self):
        if self.heartModel=='lcleeHeart':
            return self.heartModelObj.uflforms.cavitypressure()*convert_Pa_to_mmHg
        elif self.heartModel=='heArt.BiV':
            return self.heartModelObj.uflforms.LVcavitypressure()*convert_Pa_to_mmHg
        return
    def get_RVPressure(self):
        if self.heartModel=='heArt.BiV':
            return self.heartModelObj.uflforms.RVcavitypressure()*convert_Pa_to_mmHg
        return
    def get_LAPressure(self):
        return
    def get_RAPressure(self):
        return
    def get_CavityPressure(self):
        cavity=self.get_CavityNames()
        result=[]
        for chamber in cavity:
            result.append(getattr(self,'get_'+chamber+'Pressure')())
        return np.array(result)
    def get_Comm(self):
        if self.heartModel=='lcleeHeart':
            return self.heartModelObj.comm
        elif self.heartModel=='heArt.BiV':
            return self.heartModelObj.mesh_me.mpi_comm()
    def get_Mesh(self):
        if self.heartModel=='lcleeHeart':
            return self.heartModelObj.mesh
        elif self.heartModel=='heArt.BiV':
            return self.heartModelObj.Mesh.mesh
        return
    def get_MeshCoordinates(self):
        if self.heartModel=='lcleeHeart':
            return self.heartModelObj.mesh.coordinates()[:]
        elif self.heartModel=='heArt.BiV':
            return self.heartModelObj.Mesh.mesh.coordinates()[:]
        return
    def get_InverseStretchTensorFunction(self):
        if self.heartModel=='lcleeHeart':
            return self.heartModelObj.invF
    def get_InverseStretchTensorSpace(self):
        if self.heartModel=='lcleeHeart':
            return self.heartModelObj.invFspace
    def get_DisplacementFunction(self):
        if self.heartModel=='lcleeHeart':
            return self.heartModelObj.u
    def get_DisplacementSpace(self):
        if self.heartModel=='lcleeHeart':
            return self.heartModelObj.displacementSpace
    def get_DisplacementResult(self):
        if self.heartModel=='lcleeHeart':
            return self.heartModelObj.w.sub(0)
        elif self.heartModel=='heArt.BiV':
            return self.heartModelObj.GetDisplacement()
    def get_RelaxationTimeLength(self):
        if self.heartModel=='lcleeHeart':
            return fenics.project(self.heartModelObj.activeforms.tr(), self.heartModelObj.Q)
    def set_LVVolume(self,vol):
        if self.heartModel=='lcleeHeart':
            if isinstance(vol,(list,np.ndarray)):
                self.heartModelObj.Cavityvol.vol=vol[0]
            else:
                self.heartModelObj.Cavityvol.vol=vol
        elif self.heartModel=='heArt.BiV':
            self.heartModelObj.LVCavityvol.vol=vol
        return
    def set_RVVolume(self,vol):
        if self.heartModel=='heArt.BiV':
            self.heartModelObj.RVCavityvol.vol=vol
        return
    def set_LAVolume(self):
        return
    def set_RAVolume(self):
        return
    def set_CavityVolume(self,vol):
        cavity=self.get_CavityNames()
        if isinstance(vol,(int,float)):
            vol=[vol]
        count=0
        for chamber in cavity:
            getattr(self,'set_'+chamber+'Volume')(vol[count])
            count+=1
    def get_activeTime(self):
        if self.heartModel=='lcleeHeart':
            return self.heartModelObj.t_a.t_a
        elif self.heartModel=='heArt.BiV':
            return self.heartModelObj.t_a.vector()[:][0]
    def set_activeTime(self,activeTime):
        if self.heartModel=='lcleeHeart':
            self.heartModelObj.t_a.t_a=activeTime
        elif self.heartModel=='heArt.BiV':
            self.heartModelObj.t_a.vector()[:] = activeTime
            self.heartModelObj.activeforms.update_activationTime()
    def get_dTime(self):
        if self.heartModel=='lcleeHeart':
            return self.heartModelObj.dt.dt
        elif self.heartModel=='heArt.BiV':
            return self.modelSpecificRefference['state_obj'].dt.dt
    def set_dTime(self,dTime):
        if self.heartModel=='lcleeHeart':
            self.heartModelObj.dt.dt=dTime
        elif self.heartModel=='heArt.BiV':
            self.modelSpecificRefference['state_obj'].dt.dt=dTime
    def solve(self):
        if self.heartModel=='lcleeHeart':
            return self.heartModelObj.solver.solvenonlinear()
        elif self.heartModel=='heArt.BiV':
            self.modelSpecificRefference['solver'].solvenonlinear() 
    def get_CauchyFiberStressTensor(self):
        if self.heartModel=='lcleeHeart':
            cauchy1 =self.heartModelObj.uflforms.Cauchy1() + self.heartModelObj.activeforms.cauchy()
            return fenics.project(cauchy1,fenics.TensorFunctionSpace(self.get_Mesh(), "DG", 1), form_compiler_parameters={"representation":"uflacs"})
    def get_VolumeAverageFiberStress(self):
        cauchy =  self.get_CauchyFiberStressTensor()
        if self.heartModel=='lcleeHeart':
            sigma_fiber_LV = self.heartModelObj.activeforms.CalculateFiberStress(sigma = cauchy, e_fiber = self.heartModelObj.f0, Vol = self.heartModelObj.mesh_volume, Mesh = self.get_Mesh())
            return sigma_fiber_LV 
    def get_StrainEnergyDensity(self):
        if self.heartModel=='lcleeHeart':
            work1= self.heartModelObj.uflforms.strainEnergy()
            return fenics.project(work1,fenics.FunctionSpace(self.get_Mesh(), "CG", 1), form_compiler_parameters={"representation":"uflacs"})
    def get_TotalStrainEnergy(self):
        if self.heartModel=='lcleeHeart':
            return self.heartModelObj.uflforms.strainEnergy(integrate=True)
    def get_FiberDirectionVector(self):
        if self.heartModel=='lcleeHeart':
            return self.heartModelObj.f0
    def get_FiberSheetletDirectionVector(self):
        if self.heartModel=='lcleeHeart':
            return self.heartModelObj.s0
    def get_FiberNormalDirectionVector(self):
        if self.heartModel=='lcleeHeart':
            return self.heartModelObj.n0
    def get_StrainTensorFunction(self):
        if self.heartModel=='lcleeHeart':
            return self.heartModelObj.Emat()
    def get_StretchTensorFunction(self):
        if self.heartModel=='lcleeHeart':
            return self.heartModelObj.Fmat()
    def update_Phase(self):
        if self.heartModel=='lcleeHeart':
            return self.heartModelObj.updatePhase()
