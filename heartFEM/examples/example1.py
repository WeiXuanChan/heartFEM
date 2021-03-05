import heartFEM
import medImgProc.processFunc as pf
LVsim=heartFEM.LVclosed(defaultParameters={t_0:1})
LVsim.meshname = "t0"
LVsim.casename = '/home/case'
LVsim.generateMesh(Laxis,endo_angle ,epi_angle)
LVsim()

cost=cost_function()
cost.addCost('CavityPressureMeanSquare',1)
cost.addCost('CavityVolumeMeanSquare',10)

opt_link=optimiser_linker(LVsim,cost)
opt_link.addVariable(['Kspring_constant','endo_angle'])
init_variables=np.array([90.,83.])
dVariables=np.array([0.1,0.1])

gain=1.
f_error_factor=10**-4
opt_descent=pf.gradient_descent(opt_link,init_variables,gain=gain,errThreshold=dVariables,f_error=f_error_factor,limitRun=10000,finetune_space=0,normalize_para=True)
result=opt_descent.run(report=1)
