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
  Author: w.x.chan@gmail.com         06Dec2022           - v4.1.0
                                                            - added 'LV geometrical fiber strength variation file'
                                                            - added 'LV geometrical fiber delay variation file'
  Author: w.x.chan@gmail.com         06Dec2022           - v4.1.1
                                                            - debug no key 'convertParameters' in 'heArt.Electro'
  Author: w.x.chan@gmail.com         07FEB2023           - v4.2.0
                                                            - added output for heArt.BiV
  Author: w.x.chan@gmail.com         16FEB2023           - v4.3.0
                                                            - added output for heArt.BiV
  Author: w.x.chan@gmail.com         16FEB2023           - v4.3.1
                                                            - added RVvector for mesh
  Author: w.x.chan@gmail.com         05JUL2024           - v4.4.0
                                                            - added AVHeart
'''
_version='4.4.0'
import logging
logger = logging.getLogger('heartModelLinker v'+_version)
import sys
import numpy as np
import pdb
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
    if len(args)==0:
        return None
    elif args[0] is None:
        return None
    else:
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
                                                            'cstrainFactor',
                                                            'activeFiberStrength',
                                                            'activeFiberDelay',
                                                            'fiberDirectionVec',
                                                            'clipratio',
                                                            'meshsize',
                                                            'meshname',
                                                            'saveSubFolder',
                                                            'RV_vector',
                                                            'fixed_transform'],
                                                  'convertParameters':{'casePath':'path of folder with case',
                                                                       'stlname':'filename of stl',
                                                                       'Laxis':[toNumpyArray,'LV longitudinal axis superior X','LV longitudinal axis superior Y','LV longitudinal axis superior Z'],
                                                                       'endo_angle':'LV endocardial fiber angle in degrees',
                                                                       'epi_angle':'LV epicardial fiber angle in degrees',
                                                                       'fiberSheetletAngle':'LV fiber sheetlet angle in degrees',
                                                                       'fiberSheetletWidth':'LV fiber sheetlet width',
                                                                       'radialFiberAngle':'LV fiber radial angle in degrees',
                                                                       'fiberLength':'LV sarcomere length in um',
                                                                       'cstrainFactor':'LV geometrical overall stiffness constant variation file',
                                                                       'activeFiberStrength':'LV geometrical fiber strength variation file',
                                                                       'activeFiberDelay':'LV geometrical fiber delay variation file',
                                                                       'fiberDirectionVec':'LV geometrical fiber direction variation file',
                                                                       'clipratio':'ratio to clip mesh',
                                                                       'meshsize':'mesh element size factor',
                                                                       'meshname':'name of mesh',
                                                                       'saveSubFolder':'subfolder to save files',
                                                                       'RV_vector':[toNumpyArray,'LV lateral axis X','LV lateral axis Y','LV lateral axis Z'],
                                                                       'fixed_transform':'Planar basal stl file for alignment',}}
    inbuilt_changeFiberAngles_args_search['lcleeHeart']={'args':['casePath',
                                                                 'meshname',
                                                                 'endo_angle',
                                                                 'epi_angle'],
                                                'kwargs':['fiberSheetletAngle',
                                                          'fiberSheetletWidth',
                                                          'radialFiberAngle',
                                                          'cstrainFactor',
                                                          'activeFiberStrength',
                                                          'activeFiberDelay',
                                                          'fiberDirectionVec',
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
                                                                       'fiberLength':'LV sarcomere length in um',
                                                                       'cstrainFactor':'LV geometrical overall stiffness constant variation file',
                                                                       'activeFiberStrength':'LV geometrical fiber strength variation file',
                                                                       'activeFiberDelay':'LV geometrical fiber delay variation file',
                                                                       'fiberDirectionVec':'LV geometrical fiber direction variation file',
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
                                                                       'lr':'LV sarcomere length in um',
                                                                       'l0':'LV sarcomere length threshold where no tension develops in um',#[np.multiply,'LV sarcomere length in um','LV sarcomere stretch ratio threshold where no tension develops'],
                                                                       't0':'time to maximum LV fiber tension in ms',
                                                                       'Ca0':'peak intracellular calcium concentration in uM',
                                                                       'Ca0max':'maximum peak intracellular calcium concentration in uM',
                                                                       'B':'exponential coefficient for relation of peak isometric tension and sarcomere length in um-1',
                                                                       'm':'slope of linear relation of relaxation duration and sarcomere length in ms um-1',
                                                                       'b':'time intercept of linear relation of relaxation duration and sarcomere length in ms',
                                                                       'T0_LV':'maximum LV fiber tension in Pa',
                                                                       "StrainEnergyDensityFunction_Coef":"strain energy density function coefficient in Pa",
                                                                       "StrainEnergyDensityFunction_Cff":"strain energy density function exponential coefficient in fiber direction",
                                                                       "StrainEnergyDensityFunction_Css":"strain energy density function exponential coefficient in fiber sheetlet direction",
                                                                       "StrainEnergyDensityFunction_Cnn":"strain energy density function exponential coefficient in fiber normal direction",
                                                                       "StrainEnergyDensityFunction_Cns":"strain energy density function exponential coefficient in cross fiber normal and fiber sheetlet direction",
                                                                       "StrainEnergyDensityFunction_Cfs":"strain energy density function exponential coefficient in cross fiber and fiber sheetlet direction",
                                                                       "StrainEnergyDensityFunction_Cfn":"strain energy density function exponential coefficient in cross fiber and fiber normal direction",
                                                                       'EDV_LV':'LV end diastolic volume in mL',
                                                                       'EDP_LV':'LV end diastolic pressure in mmHg',
                                                                       'ESV_LV':'LV end systolic volume in mL',
                                                                       'ES_time':'duration of LV diastole in ms',
                                                                       'pressure_scale':'Solver pressure scaling normalization',
                                                                       'length_scale':'Solver length scaling normalization',
                                                                       'pinned_movement':'Solver fixed shifting factor',}}
    inbuilt_createHeartModel_args_search['lcleeHeart']['convertParameters']['runParameters']=inbuilt_createHeartModel_args_search['lcleeHeart']['convertParameters'].copy()
    
    inbulit_function_names['AVHeart']={'generateMesh':lcleeHeart.generateMesh,
                                        'changeFiberAngles':lcleeHeart.changeFiberAngles,
                                        'createHeartModel':lcleeHeart.heart_AV}
    inbuilt_generateMesh_args_search['AVHeart']={'args':['casePath',
                                                          'stlname',
                                                          'Laxis'],
                                                  'kwargs':['endo_angle',
                                                            'epi_angle',
                                                            'gmsh_surf_id',
                                                            'fiberSheetletAngle',
                                                            'fiberSheetletWidth',
                                                            'radialFiberAngle',
                                                            'fiberLength',
                                                            'cstrainFactor',
                                                            'activeFiberStrength',
                                                            'activeFiberDelay',
                                                            'fiberDirectionVec',
                                                            'clipratio',
                                                            'meshsize',
                                                            'meshname',
                                                            'saveSubFolder',
                                                            'RV_vector',
                                                            'fixed_transform',
                                                            'split_avj_extend',
                                                            'split_avj_refinement',
                                                            'model'],
                                                  'convertParameters':{'casePath':'path of folder with case',
                                                                       'stlname':'filename of stl',
                                                                       'Laxis':[toNumpyArray,'LV longitudinal axis superior X','LV longitudinal axis superior Y','LV longitudinal axis superior Z'],
                                                                       'endo_angle':'LV endocardial fiber angle in degrees',
                                                                       'epi_angle':'LV epicardial fiber angle in degrees',
                                                                       'gmsh_surf_id':[toNumpyArray,'gmsh LA endocardial surface ID','gmsh LV endocardial surface ID','gmsh LA epicardial surface ID','gmsh LV epicardial surface ID','gmsh AV juction surface ID'],
                                                                       'fiberSheetletAngle':'LV fiber sheetlet angle in degrees',
                                                                       'fiberSheetletWidth':'LV fiber sheetlet width',
                                                                       'radialFiberAngle':'LV fiber radial angle in degrees',
                                                                       'fiberLength':'LV sarcomere length in um',
                                                                       'cstrainFactor':'LV geometrical overall stiffness constant variation file',
                                                                       'activeFiberStrength':'LV geometrical fiber strength variation file',
                                                                       'activeFiberDelay':'LV geometrical fiber delay variation file',
                                                                       'fiberDirectionVec':'LV geometrical fiber direction variation file',
                                                                       'clipratio':'ratio to clip mesh',
                                                                       'meshsize':'mesh element size factor',
                                                                       'meshname':'name of mesh',
                                                                       'saveSubFolder':'subfolder to save files',
                                                                       'RV_vector':[toNumpyArray,'LV lateral axis X','LV lateral axis Y','LV lateral axis Z'],
                                                                       'fixed_transform':'Planar basal stl file for alignment',
                                                                       'split_avj_extend':'Distance to extend Mitral Plane to create Mitral chamber in mm',
                                                                       'split_avj_refinement':[toNumpyArray,'gmsh AV juction refinement min size factor','gmsh AV juction refinement max size factor','gmsh AV juction refinement min distance factor','gmsh AV juction refinement max distance factor'],
                                                                       'model':'AVsplit'}}
    inbuilt_changeFiberAngles_args_search['AVHeart']={'args':['casePath',
                                                                 'meshname',
                                                                 'endo_angle',
                                                                 'epi_angle'],
                                                'kwargs':['fiberSheetletAngle',
                                                          'fiberSheetletWidth',
                                                          'radialFiberAngle',
                                                          'cstrainFactor',
                                                          'activeFiberStrength',
                                                          'activeFiberDelay',
                                                          'fiberDirectionVec',
                                                          'fiberLength',
                                                          'saveSubFolder',
                                                          'loadSubFolder',
                                                          'model',
                                                          'surfID'],
                                                  'convertParameters':{'casePath':'path of folder with case',
                                                                       'meshname':'name of mesh',
                                                                       'endo_angle':'LV endocardial fiber angle in degrees',
                                                                       'epi_angle':'LV epicardial fiber angle in degrees',
                                                                       'fiberSheetletAngle':'LV fiber sheetlet angle in degrees',
                                                                       'fiberSheetletWidth':'LV fiber sheetlet width',
                                                                       'radialFiberAngle':'LV fiber radial angle in degrees',
                                                                       'fiberLength':'LV sarcomere length in um',
                                                                       'cstrainFactor':'LV geometrical overall stiffness constant variation file',
                                                                       'activeFiberStrength':'LV geometrical fiber strength variation file',
                                                                       'activeFiberDelay':'LV geometrical fiber delay variation file',
                                                                       'fiberDirectionVec':'LV geometrical fiber direction variation file',
                                                                       'clipratio':'ratio to clip mesh',
                                                                       'meshsize':'mesh element size factor',
                                                                       'saveSubFolder':'subfolder to save files',
                                                                       'loadSubFolder':'subfolder to load files',
                                                                       'model':'AVsplit',
                                                                       'surfID':['LA epicardial mesh surface ID',
                                                                                 'LV epicardial mesh surface ID',
                                                                                 'LA endocardial mesh surface ID',
                                                                                 'LV endocardial mesh surface ID',
                                                                                 'LA myocardial boundary mesh surface ID',
                                                                                 'LV myocardial boundary mesh surface ID',
                                                                                 'AV juction mesh surface ID',]}}
    inbuilt_createHeartModel_args_search['AVHeart']={'args':['casename',
                                                                'meshname',
                                                                'runParameters'],
                                                'kwargs':['inverseHeart',
                                                          'trackphase',
                                                          'AVJ_control'],
                                                  'convertParameters':{'casename':'path of folder with case',
                                                                       'meshname':'name of mesh',
                                                                       'trackphase':'fiber relaxation based on phase during FEA',
                                                                       'AVJ_control':'Mitral valve channel control',
                                                                       'inverseHeart':'heart to solve for inverse',
                                                                       'avjid':'AV juction mesh surface ID',
                                                                       'endoAtrid':'LA endocardial mesh surface ID',
                                                                       'epiAtrid':'LA epicardial mesh surface ID',
                                                                       'endoVenid':'LV endocardial mesh surface ID',
                                                                       'epiVenid':'LV epicardial mesh surface ID',
                                                                       'endoAtrAVJid':'LA junction mesh surface ID',
                                                                       'endoVenAVJid':'LV junction mesh surface ID',
                                                                       'BCL':'duration of one cardiac cycle in ms',
                                                                       'Kspring_constant':'spring constant for LV pericardial cavity in Pa',
                                                                       'lr':'LV sarcomere length in um',
                                                                       'l0':'LV sarcomere length threshold where no tension develops in um',#[np.multiply,'LV sarcomere length in um','LV sarcomere stretch ratio threshold where no tension develops'],
                                                                       't0':'time to maximum LV fiber tension in ms',
                                                                       'Ca0':'peak intracellular calcium concentration in uM',
                                                                       'Ca0max':'maximum peak intracellular calcium concentration in uM',
                                                                       'B':'exponential coefficient for relation of peak isometric tension and sarcomere length in um-1',
                                                                       'm':'slope of linear relation of relaxation duration and sarcomere length in ms um-1',
                                                                       'b':'time intercept of linear relation of relaxation duration and sarcomere length in ms',
                                                                       'T0_LV':'maximum LV fiber tension in Pa',
                                                                       "StrainEnergyDensityFunction_Coef":"strain energy density function coefficient in Pa",
                                                                       "StrainEnergyDensityFunction_Cff":"strain energy density function exponential coefficient in fiber direction",
                                                                       "StrainEnergyDensityFunction_Css":"strain energy density function exponential coefficient in fiber sheetlet direction",
                                                                       "StrainEnergyDensityFunction_Cnn":"strain energy density function exponential coefficient in fiber normal direction",
                                                                       "StrainEnergyDensityFunction_Cns":"strain energy density function exponential coefficient in cross fiber normal and fiber sheetlet direction",
                                                                       "StrainEnergyDensityFunction_Cfs":"strain energy density function exponential coefficient in cross fiber and fiber sheetlet direction",
                                                                       "StrainEnergyDensityFunction_Cfn":"strain energy density function exponential coefficient in cross fiber and fiber normal direction",
                                                                       'EDV_LV':'LV end diastolic volume in mL',
                                                                       'EDP_LV':'LV end diastolic pressure in mmHg',
                                                                       'ESV_LV':'LV end systolic volume in mL',
                                                                       'ES_time':'duration of LV diastole in ms',
                                                                       'pressure_scale':'Solver pressure scaling normalization',
                                                                       'length_scale':'Solver length scaling normalization',
                                                                       'pinned_movement':'Solver fixed shifting factor',
                                                                       "solver type" : "FEniCS solver type"}}
    inbuilt_createHeartModel_args_search['AVHeart']['convertParameters']['runParameters']=inbuilt_createHeartModel_args_search['AVHeart']['convertParameters'].copy()
        
if 'heArt' in locals():
    inbulit_function_names['heArt.BiV']={'generateMesh':lcleeHeart.generateMesh,
                                        'changeFiberAngles':lcleeHeart.changeFiberAngles,
                                        'createHeartModel':heArt.MEmodel}
    inbuilt_generateMesh_args_search['heArt.BiV']=inbuilt_generateMesh_args_search['lcleeHeart']
    inbuilt_changeFiberAngles_args_search['heArt.BiV']=inbuilt_changeFiberAngles_args_search['lcleeHeart']
    inbuilt_createHeartModel_args_search['heArt.BiV']={'args':['params',
                                                               'SimDet'],
                                                'kwargs':[],
                                                  'convertParameters':{'SimDet':{"HeartBeatLength": 'duration of one cardiac cycle in ms',
                                                                                  "dt": 1.0,
                                                                                  "writeStep": 10.0,
                                                                                  "GiccioneParams" : {"ParamsSpecified" : True, 
                                                                                             		  "Passive model": {"Name": "Guccione"},
                                                                                             		  "Passive params": {"Cparam": [fenics.Constant,"strain energy density function coefficient in Pa"], 
                                                                                        		                     "bff"  : [fenics.Constant,"strain energy density function exponential coefficient in fiber direction"],
                                                                                        		                     "bfx"  : [fenics.Constant,"strain energy density function exponential coefficient in fiber sheetlet direction"],
                                                                                                          		     "bxx"  : [fenics.Constant,"strain energy density function exponential coefficient in cross fiber normal and fiber sheetlet direction"]},
                                                                                             		  "Active model": {"Name": "Guccione"},
                                                                                             		  "Active params": {'m':'slope of linear relation of relaxation duration and sarcomere length in ms um-1',
                                                                                                                        'b':'time intercept of linear relation of relaxation duration and sarcomere length in ms',
                                                                                                                        "B" : 'maximum peak intracellular calcium concentration in uM', 
                                                                                                                        "t0" : 'time to maximum LV fiber tension in ms', 
                                                                                                                        "l0" : 'LV sarcomere length threshold where no tension develops in um',
                                                                                                                        "Tmax" : [fenics.Constant,'maximum LV fiber tension in Pa'],
                                                                                                                        "Ca0" : 'peak intracellular calcium concentration in uM',
                                                                                                                        "Ca0max" : 'maximum peak intracellular calcium concentration in uM',
                                                                                                                        "lr" : 'LV sarcomere length in um'},
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
                                                                                  "abs_tol": "solving tolerance",
                                                                                  "rel_tol": "solving tolerance",}}}
    inbuilt_createHeartModel_args_search['heArt.BiV']['convertParameters']['params']={"directory" : [add,'path of folder with case','/'],
                                                                                 "casename" : 'name of mesh', 
                                                                                 "fibre_quad_degree" : 4, 
                                                                                 "outputfolder" : 'subfolder to save files',
                                                                                 "foldername" : 'subfolder to load files',
                                                                                 "state_obj": [heArt.State_Variables,'common communicator', inbuilt_createHeartModel_args_search['heArt.BiV']['convertParameters']['SimDet']],
                                                                                 "common_communicator": 'common communicator',
                                                                                 "MEmesh": [fenics.Mesh],
                                                                                 "isLV": False}
    inbulit_function_names['heArt.LV']={'generateMesh':lcleeHeart.generateMesh,
                                        'changeFiberAngles':lcleeHeart.changeFiberAngles,
                                        'createHeartModel':heArt.MEmodel}
    inbuilt_generateMesh_args_search['heArt.LV']=inbuilt_generateMesh_args_search['heArt.BiV']
    inbuilt_changeFiberAngles_args_search['heArt.LV']=inbuilt_changeFiberAngles_args_search['heArt.BiV']
    inbuilt_createHeartModel_args_search['heArt.LV']=inbuilt_createHeartModel_args_search['heArt.BiV']
    
    
    inbulit_function_names['heArt.Electro']={'createHeartModel':heArt.EPmodel}
    inbuilt_createHeartModel_args_search['heArt.Electro']={'args':['params'],
                                                'kwargs':[]}
    inbuilt_createHeartModel_args_search['heArt.Electro']['convertParameters']={'params':{"casename" : 'name of mesh',
                                                                                      "EPmesh":'full path to electromechanical mesh',
                                                                                      "DTI_EP":False,
                                                                                      "dt": 1.0,
                                                                                      "deg": 4,
                                                                                      "ploc":"electro pacing location",
                                                                                      "pacing_timing":"electro pacing timing in ms"}}
    
    
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
        self.Mitral=False
        self.Aortic=False
        self.Tricuspid=False
        self.Pulmonic=False
        self.modelSpecificRefference={}
        self.heartModelObj=functionTemplete('createHeartModel',heartModel,params)
        if self.heartModel=='heArt.BiV':
            self.RV=True
        if self.heartModel=='AVHeart':
            self.LA=True
            if params['Mitral valve channel control']!=0:
                self.Mitral=True
        elif self.heartModel in ['lcleeHeart','AVHeart']:
            self.comm=self.heartModelObj.comm
            self.mesh=self.heartModelObj.mesh
        if self.heartModel in ['lcleeHeart','AVHeart']:
            self.modelsolver=self.heartModelObj.solver
        elif self.heartModel in ['heArt.BiV','heArt.LV']:
            self.modelsolver=self.heartModelObj.Solver()
            self.modelSpecificRefference['state_obj']=self.heartModelObj.state_obj
        elif self.heartModel=='heArt.Electro':
            self.modelsolver=self.heartModelObj.Solver()
            self.modelSpecificRefference['state_obj']=self.heartModelObj.parameters["state_obj"]
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
        if self.Mitral:
            cavity.append("Mitral")
        if self.Aortic:
            cavity.append("Aortic")
        if self.Tricuspid:
            cavity.append("Tricuspid")
        if self.Pulmonic:
            cavity.append("Pulmonic")
        return cavity
    def get_LVVolume(self):
        if self.heartModel=='lcleeHeart':
            return self.heartModelObj.uflforms.cavityvol()
        elif self.heartModel=='AVHeart':
            return self.heartModelObj.uflforms.cavityvol('Ven')
        elif self.heartModel in ['heArt.BiV','heArt.LV']:
            return self.heartModelObj.uflforms.LVcavityvol()
        return
    def get_RVVolume(self):
        if self.heartModel=='heArt.BiV':
            return self.heartModelObj.uflforms.RVcavityvol()
        return
    def get_LAVolume(self):
        if self.heartModel=='AVHeart':
            return self.heartModelObj.uflforms.cavityvol('Atr')
        return
    def get_RAVolume(self):
        return
    def get_MitralVolume(self):
        if self.heartModel=='AVHeart':
            return self.heartModelObj.uflforms.avjvol()
        return
    def get_AorticVolume(self):
        return
    def get_TricuspidVolume(self):
        return
    def get_PulmonicVolume(self):
        return
    def get_CavityVolume(self,verbose=False):
        cavity=self.get_CavityNames()
        if verbose:
            chamberNames=[]
        result=[]
        for chamber in cavity:
            result.append(getattr(self,'get_'+chamber+'Volume')())
            if verbose:
                chamberNames.append(chamber)
        if verbose:
            logger.info("Cavity Volume : "+repr({chamberNames[i]: result[i] for i in range(len(chamberNames))}))
        return np.array(result)
    def get_LVPressure(self):
        if self.heartModel=='lcleeHeart':
            return self.heartModelObj.uflforms.cavitypressure()*convert_Pa_to_mmHg
        elif self.heartModel=='AVHeart':
            return self.heartModelObj.uflforms.cavitypressure('Ven')*convert_Pa_to_mmHg
        elif self.heartModel in ['heArt.BiV','heArt.LV']:
            return self.heartModelObj.uflforms.LVcavitypressure()*convert_Pa_to_mmHg
        return
    def get_RVPressure(self):
        if self.heartModel=='heArt.BiV':
            return self.heartModelObj.uflforms.RVcavitypressure()*convert_Pa_to_mmHg
        return
    def get_LAPressure(self):
        if self.heartModel=='AVHeart':
            return self.heartModelObj.uflforms.cavitypressure('Atr')*convert_Pa_to_mmHg
        return
    def get_RAPressure(self):
        return
    def get_MitralPressure(self):
        if self.heartModel=='AVHeart':
            return self.heartModelObj.uflforms.cavitypressure('AVJ')*convert_Pa_to_mmHg
        return
    def get_AorticPressure(self):
        return
    def get_TricuspidPressure(self):
        return
    def get_PulmonicPressure(self):
        return
    def get_CavityPressure(self,verbose=False):
        cavity=self.get_CavityNames()
        if verbose:
            chamberNames=[]
        result=[]
        for chamber in cavity:
            result.append(getattr(self,'get_'+chamber+'Pressure')())
            if verbose:
                chamberNames.append(chamber)
        if verbose:
            logger.info("Cavity Pressure : "+repr({chamberNames[i]: result[i] for i in range(len(chamberNames))}))
        return np.array(result)
    def get_Comm(self):
        if self.heartModel in ['lcleeHeart','AVHeart']:
            return self.heartModelObj.comm
        elif self.heartModel in ['heArt.BiV','heArt.LV']:
            return self.heartModelObj.mesh_me.mpi_comm()
    def get_Mesh(self):
        if self.heartModel in ['lcleeHeart','AVHeart']:
            return self.heartModelObj.mesh
        elif self.heartModel in ['heArt.BiV','heArt.LV']:
            return self.heartModelObj.Mesh.mesh
        return
    def get_MeshCoordinates(self):
        if self.heartModel in ['lcleeHeart','AVHeart']:
            return self.heartModelObj.mesh.coordinates()[:]
        elif self.heartModel in ['heArt.BiV','heArt.LV']:
            return self.heartModelObj.Mesh.mesh.coordinates()[:]
    def get_InverseStretchTensorFunction(self):
        if self.heartModel in ['lcleeHeart','AVHeart']:
            return self.heartModelObj.invF
    def get_InverseStretchTensorSpace(self):
        if self.heartModel in ['lcleeHeart','AVHeart']:
            return self.heartModelObj.invFspace
    def get_DisplacementFunction(self):
        if self.heartModel in ['lcleeHeart','AVHeart']:
            return self.heartModelObj.u
    def get_DisplacementSpace(self):
        if self.heartModel in ['lcleeHeart','AVHeart']:
            return self.heartModelObj.displacementSpace
    def get_DisplacementResult(self):
        if self.heartModel in ['lcleeHeart','AVHeart']:
            return self.heartModelObj.w.sub(0)
        elif self.heartModel in ['heArt.BiV','heArt.LV']:
            return self.heartModelObj.GetDisplacement()
    def get_ScalarFunctionSpace(self):
        if self.heartModel in ['lcleeHeart','AVHeart']:
            return self.heartModelObj.Q
        elif self.heartModel in ['heArt.BiV','heArt.LV']:
            return fenics.FunctionSpace(self.get_Mesh(), "CG", 1)
    def get_ActiveCalciumStrain(self):
        if self.heartModel in ['lcleeHeart','AVHeart']:
            return self.heartModelObj.activeforms.ECa()
        elif self.heartModel in ['heArt.BiV','heArt.LV']:
            return self.heartModelObj.activeforms.activeforms.ECa()
    def get_RelaxationTimeLength(self):
        if self.heartModel in ['lcleeHeart','AVHeart']:
            return fenics.project(self.heartModelObj.activeforms.tr(), self.get_ScalarFunctionSpace())
        elif self.heartModel in ['heArt.BiV','heArt.LV']:
            return fenics.project(self.heartModelObj.activeforms.activeforms.tr(), self.get_ScalarFunctionSpace(),form_compiler_parameters={"representation":"uflacs","quadrature_degree":4})
    def set_LVVolume(self,vol):
        if self.heartModel=='lcleeHeart':
            if isinstance(vol,(list,np.ndarray)):
                self.heartModelObj.Cavityvol.vol=vol[0]
            else:
                self.heartModelObj.Cavityvol.vol=vol
        elif self.heartModel=='AVHeart':
            logger.debug('set_LVVolume with'+repr(vol))
            if isinstance(vol,(list,np.ndarray)):
                self.heartModelObj.CavityvolVen.vol=vol[0]
            else:
                self.heartModelObj.CavityvolVen.vol=vol
        elif self.heartModel in ['heArt.BiV','heArt.LV']:
            self.heartModelObj.LVCavityvol.vol=vol
        return
    def set_RVVolume(self,vol):
        if self.heartModel=='heArt.BiV':
            self.heartModelObj.RVCavityvol.vol=vol
        return
    def set_LAVolume(self,vol):
        if self.heartModel=='AVHeart':
            logger.debug('set_LAVolume with'+repr(vol))
            if isinstance(vol,(list,np.ndarray)):
                self.heartModelObj.CavityvolAtr.vol=vol[0]
            else:
                self.heartModelObj.CavityvolAtr.vol=vol
        return
    def set_RAVolume(self):
        return
    def set_MitralVolume(self,vol):
        if self.heartModel=='AVHeart':
            logger.debug('set_MitralVolume with'+repr(vol))
            if isinstance(vol,(list,np.ndarray)):
                self.heartModelObj.CavityvolAVJ.vol=vol[0]
            else:
                self.heartModelObj.CavityvolAVJ.vol=vol
        return
    def set_AorticVolume(self):
        return
    def set_TricuspidVolume(self):
        return
    def set_PulmonicVolume(self):
        return
    def set_CavityVolume(self,vol):
        cavity=self.get_CavityNames()
        if isinstance(vol,(int,float)):
            vol=[vol]
        count=0
        for chamber in cavity:
            getattr(self,'set_'+chamber+'Volume')(vol[count])
            count+=1
    def get_ReferenceFiberStrength(self):
        if self.heartModel in ['lcleeHeart','AVHeart']:
            return self.heartModelObj.activeforms.parameters["T0"]
        elif self.heartModel in ['heArt.BiV','heArt.LV']:
            return self.heartModelObj.activeforms.activeforms.parameters["material params"]["Tmax"]
    def set_ReferenceFiberStrength(self,referenceFiberStrength):
        if self.heartModel in ['lcleeHeart','AVHeart']:
            self.heartModelObj.activeforms.parameters["T0"]=referenceFiberStrength
        elif self.heartModel in ['heArt.BiV','heArt.LV']:
            self.heartModelObj.SimDet["GiccioneParams"]["Active model"]["Tmax"] = fenics.Constant(referenceFiberStrength)
            self.heartModelObj.activeforms.activeforms.parameters["material params"]["Tmax"]=self.heartModelObj.SimDet["GiccioneParams"]["Active model"]["Tmax"]
    def get_activeTime(self):
        if self.heartModel in ['lcleeHeart','AVHeart']:
            return self.heartModelObj.t_a.t_a
        elif self.heartModel in ['heArt.BiV','heArt.LV']:
            return self.heartModelObj.t_a.vector()[:][0]
    def set_activeTime(self,activeTime):
        if self.heartModel in ['lcleeHeart','AVHeart']:
            self.heartModelObj.t_a.t_a=activeTime
        elif self.heartModel in ['heArt.BiV','heArt.LV']:
            self.heartModelObj.t_a.vector()[:] = activeTime
            self.heartModelObj.activeforms.update_activationTime()
    def get_dTime(self):
        if self.heartModel in ['lcleeHeart','AVHeart']:
            return self.heartModelObj.dt.dt
        elif self.heartModel in ['heArt.BiV','heArt.LV']:
            return self.modelSpecificRefference['state_obj'].dt.dt
    def set_dTime(self,dTime):
        if self.heartModel in ['lcleeHeart','AVHeart']:
            self.heartModelObj.dt.dt=dTime
        elif self.heartModel in ['heArt.BiV','heArt.LV']:
            self.modelSpecificRefference['state_obj'].dt.dt=dTime
    def solve(self):
        return self.modelsolver.solvenonlinear()
    def get_CauchyFiberStressTensor(self):
        if self.heartModel in ['lcleeHeart','AVHeart']:
            cauchy1 =self.heartModelObj.uflforms.Cauchy1() + self.heartModelObj.activeforms.cauchy()
            return fenics.project(cauchy1,fenics.TensorFunctionSpace(self.get_Mesh(), "DG", 1), form_compiler_parameters={"representation":"uflacs"})
        elif self.heartModel in ['heArt.BiV','heArt.LV']:
            cauchy1 =self.heartModelObj.uflforms.passiveforms.Cauchy1() + self.heartModelObj.activeforms.cauchy()
            return fenics.project(cauchy1,fenics.TensorFunctionSpace(self.get_Mesh(), "DG", 1), form_compiler_parameters={"representation":"uflacs","quadrature_degree":4})
    def get_VolumeAverageFiberStress(self):
        cauchy =  self.get_CauchyFiberStressTensor()
        if self.heartModel in ['lcleeHeart','AVHeart']:
            sigma_fiber_LV = self.heartModelObj.activeforms.CalculateFiberStress(sigma = cauchy, e_fiber = self.heartModelObj.f0, Vol = self.heartModelObj.mesh_volume, Mesh = self.get_Mesh())
            return sigma_fiber_LV 
    def get_StrainEnergyDensity(self):
        if self.heartModel in ['lcleeHeart','AVHeart','heArt.BiV','heArt.LV']:
            work1= self.heartModelObj.uflforms.strainEnergy()
            return fenics.project(work1,fenics.FunctionSpace(self.get_Mesh(), "CG", 1), form_compiler_parameters={"representation":"uflacs"})
    def get_TotalStrainEnergy(self):
        if self.heartModel in ['lcleeHeart','AVHeart']:
            return self.heartModelObj.uflforms.strainEnergy(integrate=True)
        elif self.heartModel in ['heArt.BiV','heArt.LV']:
            return self.heartModelObj.uflforms.strainEnergy(integrate=True)
    def get_FiberDirectionVector(self):
        if self.heartModel in ['lcleeHeart','AVHeart']:
            return self.heartModelObj.f0
        elif self.heartModel in ['heArt.BiV','heArt.LV']:
            return self.heartModelObj.f0_me
    def get_FiberSheetletDirectionVector(self):
        if self.heartModel in ['lcleeHeart','AVHeart']:
            return self.heartModelObj.s0
        elif self.heartModel in ['heArt.BiV','heArt.LV']:
            return self.heartModelObj.s0_me
    def get_FiberNormalDirectionVector(self):
        if self.heartModel in ['lcleeHeart','AVHeart']:
            return self.heartModelObj.n0
        elif self.heartModel in ['heArt.BiV','heArt.LV']:
            return self.heartModelObj.n0_me
    def get_FiberFunctionSpace(self):
        if self.heartModel in ['lcleeHeart','AVHeart']:
            return self.heartModelObj.fiberFS
    def get_LongitudinalDirectionVector(self):
        if self.heartModel in ['lcleeHeart','AVHeart']:
            return self.heartModelObj.long0
        elif self.heartModel in ['heArt.BiV','heArt.LV']:
            return self.heartModelObj.Mesh.long0
    def get_CircumferentialDirectionVector(self):
        if self.heartModel in ['lcleeHeart','AVHeart']:
            return self.heartModelObj.cir0
        elif self.heartModel in ['heArt.BiV','heArt.LV']:
            return self.heartModelObj.Mesh.cir0
    def get_RadialDirectionVector(self):
        if self.heartModel in ['lcleeHeart','AVHeart']:
            return self.heartModelObj.rad0
        elif self.heartModel in ['heArt.BiV','heArt.LV']:
            return self.heartModelObj.Mesh.rad0
    def get_StrainTensorFunction(self):
        if self.heartModel in ['lcleeHeart','AVHeart']:
            return self.heartModelObj.uflforms.Emat()
        elif self.heartModel in ['heArt.BiV','heArt.LV']:
            return self.heartModelObj.uflforms.Emat()
    def get_StretchTensorFunction(self):
        if self.heartModel in ['lcleeHeart','AVHeart']:
            return self.heartModelObj.uflforms.Fmat()
        elif self.heartModel in ['heArt.BiV','heArt.LV']:
            return self.heartModelObj.uflforms.Fmat()
    def get_RightCauchyGreenDeformationTensorFunction(self):
        if self.heartModel in ['lcleeHeart','heArt.BiV','heArt.LV','AVHeart']:
            Fmat=self.get_StretchTensorFunction()
            return (Fmat.T*Fmat)
    def get_ActivationDelay(self):
        if self.heartModel in ['heArt.BiV','heArt.LV']:
            return self.heartModelObj.activeforms.t_init
    def update_Phase(self):
        if self.heartModel in ['lcleeHeart','AVHeart']:
            return self.heartModelObj.updatePhase()
        elif self.heartModel in ['heArt.BiV','heArt.LV','heArt.Electro']:
            return self.heartModelObj.UpdateVar()
    def reset_Phase(self):
        if self.heartModel=='heArt.Electro':
            self.heartModelObj.reset()
    def update_activationTime(self,state):
        if self.heartModel in ['heArt.BiV','heArt.LV','heArt.Electro']:
            potential_me = fenics.Function(fenics.FunctionSpace(self.get_Mesh(),'CG',1))
            potential_me.vector()[:] = state.vector().array()[:]
            self.heartModelObj.activeforms.update_activationTime(potential_n = potential_me, comm = self.get_Comm())
    def get_interpolate_potential_ep2me_phi(self,mesh=None):
        if mesh is None:
            mesh=self.get_Mesh()
        if self.heartModel=='heArt.Electro':
            potential_ref = self.heartModelObj.interpolate_potential_ep2me_phi(V_me = fenics.Function(fenics.FunctionSpace(mesh,'CG',1)))
            potential_ref.rename("v_ref", "v_ref")
        else:
            raise Exception("get_interpolate_potential_ep2me_phi is only usable with heArt.Electro.")
        return potential_ref
    def get_State_object(self):
        if self.heartModel=='heArt.Electro':
            return self.heartModelObj.parameters["state_obj"]
        elif self.heartModel in ['heArt.BiV','heArt.LV']:
            return self.heartModelObj.state_obj
        else:
            raise Exception("No State Oject for this model.")
    def replace_State_object(self,state_obj):
        if self.heartModel=='heArt.Electro':
            self.heartModelObj.parameters["state_obj"]= state_obj
            self.modelSpecificRefference['state_obj']=self.heartModelObj.parameters["state_obj"]
        elif self.heartModel in ['heArt.BiV','heArt.LV']:
            self.heartModelObj.state_obj= state_obj
            self.modelSpecificRefference['state_obj']=self.heartModelObj.state_obj
        else:
            raise Exception("No State Oject for this model.")
    def get_current_fiber_vector(self):
        if self.heartModel in ['lcleeHeart','AVHeart']:
            fibervec=self.heartModelObj.uflforms.currentFiberVector()
            #return fenics.project(fibervec,fenics.VectorFunctionSpace(self.get_Mesh(),  "CG", 1), form_compiler_parameters={"representation":"uflacs"})
            return fibervec
    def set_debug_call(self,debug_call,*args):
        if self.heartModel in ['lcleeHeart','AVHeart']:
            self.heartModelObj.debug_call.append([debug_call]+list(args))
            #return fenics.project(fibervec,fenics.VectorFunctionSpace(self.get_Mesh(),  "CG", 1), form_compiler_parameters={"representation":"uflacs"})
            return "Added debug call. Total number of debug call: "+str(len(self.heartModelObj.debug_call))
    def get_debug_call(self):
        if self.heartModel in ['lcleeHeart','AVHeart']:
            return self.heartModelObj.debug_call
