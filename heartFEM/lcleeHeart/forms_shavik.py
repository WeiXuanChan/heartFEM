from dolfin import *
from ufl import cofac
import dolfin as dolfin
import numpy as np
from ufl import k
from ufl import i
from ufl import j
#from mpi4py import MPI					 
import sys

class Forms(object):

	def __init__(self, params):		

		 	self.parameters = params

	def Fmat(self,inverse=False):
		u = self.parameters["displacement_variable"]
		d = u.geometric_dimension()
		I = Identity(d)
		if inverse:
			invF=self.parameters["invF_variable"]
			F=invF*(I + grad(u)) - I
		else:
		    F = I + grad(u)
		return F
	
	def Emat(self,inverse=False):
	 
		u = self.parameters["displacement_variable"]
		d = u.geometric_dimension()
		I = Identity(d)
		F = self.Fmat(inverse=inverse)
		return 0.5*(F.T*F-I)


	def J(self,inverse=False):
		F = self.Fmat(inverse=inverse)
		return det(F)	
	
	
	def cavityvol(self):
	
		u = self.parameters["displacement_variable"]
		N = self.parameters["facet_normal"]
		mesh = self.parameters["mesh"]
		X = SpatialCoordinate(mesh)
		ds = dolfin.ds(subdomain_data=self.parameters["facetboundaries"])
			
		F = self.Fmat()
			
		vol_form = -Constant(1.0/3.0) * inner(det(F)*dot(inv(F).T, N), X + u)*ds(self.parameters["endoid"])
			
		return assemble(vol_form)

	def cavitypressure(self):
	
		W = self.parameters["mixedfunctionspace"]
		w = self.parameters["mixedfunction"]
		mesh = self.parameters["mesh"]

		comm = W.mesh().mpi_comm()
		dofmap =  W.sub(2).dofmap()
		val_dof = dofmap.cell_dofs(0)[0]

			# the owner of the dof broadcasts the value
		own_range = dofmap.ownership_range()
	
		try:
			val_local = w.vector()[val_dof]				  
		except IndexError:
			val_local = 0.0


		pressure = MPI.sum(comm, val_local)

	
  
		# Serial ###############################################
				#dofmap =  W.sub(2).dofmap()
			#val_dof = dofmap.cell_dofs(0)[0]
			#pressure = w.vector()[val_dof][0]
		########################################################
	
		return pressure


	def PassiveMatSEF(self,inverse=False):
		Ea = self.Emat()
		Ea_inv = self.Emat(inverse=inverse)
		f0 = self.parameters["fiber"]
		s0 = self.parameters["sheet"]
		n0 = self.parameters["sheet-normal"]
		p = self.parameters["pressure_variable"]
        
		Cff=self.parameters["StrainEnergyDensityFunction_Cff"]
		Css=self.parameters["StrainEnergyDensityFunction_Css"]
		Cnn=self.parameters["StrainEnergyDensityFunction_Cnn"]
		Cns=self.parameters["StrainEnergyDensityFunction_Cns"]
		Cfs=self.parameters["StrainEnergyDensityFunction_Cfs"]
		Cfn=self.parameters["StrainEnergyDensityFunction_Cfn"]

		Eff = inner(f0, Ea*f0)
		Ess = inner(s0, Ea*s0)
		Enn = inner(n0, Ea*n0)
		Efs = inner(f0, Ea*s0)
		Efn = inner(f0, Ea*n0)
		Ens = inner(n0, Ea*s0)
        
		
		Eff_inv = inner(f0, Ea_inv*f0)
		Ess_inv = inner(s0, Ea_inv*s0)
		Enn_inv = inner(n0, Ea_inv*n0)
		Efs_inv = inner(f0, Ea_inv*s0)
		Efn_inv = inner(f0, Ea_inv*n0)
		Ens_inv = inner(n0, Ea_inv*s0)
		
		
		QQ = Cff*pow(Eff,2.0) + Css*pow(Ess,2.0)+ Cnn*pow(Enn,2.0)+ Cns*pow(Ens,2.0) + Cfs*pow(Efs,2.0) + Cfn*pow(Efn,2.0) #Original
		if inverse:  
		    #pendo = self.parameters["volconst_variable"] 
		    #V0= self.parameters["constrained_vol"] 
		    #QQ_inv = 29.9*pow(V0*Eff_inv,2.0) + 13.3*(pow(V0*Ess_inv,2.0)+ pow(V0*Enn_inv,2.0)+ 2.0*pow(V0*Ens_inv,2.0)) + 26.6*(2.0*pow(V0*Efs_inv,2.0) + 2.0*pow(V0*Efn_inv,2.0))
		    Wp=[Eff_inv,Ess_inv,Enn_inv,Efs_inv,Efn_inv,Ens_inv] #inverse
		else:
		    Wp =(100*(exp(QQ) -  1.0) - p*(self.J() - 1.0))*dx(self.parameters["mesh"]) #original
		return Wp

	def strainEnergy(self, integrate=False):
	    Ea = self.Emat()
	    f0 = self.parameters["fiber"]
	    s0 = self.parameters["sheet"]
	    n0 = self.parameters["sheet-normal"]
	    p = self.parameters["pressure_variable"]
        
	    Cff=self.parameters["StrainEnergyDensityFunction_Cff"]
	    Css=self.parameters["StrainEnergyDensityFunction_Css"]
	    Cnn=self.parameters["StrainEnergyDensityFunction_Cnn"]
	    Cns=self.parameters["StrainEnergyDensityFunction_Cns"]
	    Cfs=self.parameters["StrainEnergyDensityFunction_Cfs"]
	    Cfn=self.parameters["StrainEnergyDensityFunction_Cfn"]

	    Eff = inner(f0, Ea*f0)
	    Ess = inner(s0, Ea*s0)
	    Enn = inner(n0, Ea*n0)
	    Efs = inner(f0, Ea*s0)
	    Efn = inner(f0, Ea*n0)
	    Ens = inner(n0, Ea*s0)
        
		
	    QQ = Cff*pow(Eff,2.0) + Css*pow(Ess,2.0)+ Cnn*pow(Enn,2.0)+ Cns*pow(Ens,2.0) + Cfs*pow(Efs,2.0) + Cfn*pow(Efn,2.0) #Original
	    if integrate:
		    Wp = 0.0 
		    Wp = assemble((100*(exp(QQ) -  1.0))*dx(self.parameters["mesh"]),form_compiler_parameters={"representation":"uflacs"})
	    else:
		    Wp = 100*(exp(QQ) -  1.0)
	    return Wp

	def EmatECC(self):
		Ea = self.Emat()
		CC0 = self.parameters["circumferentialvector"]
		EmatECC = inner(CC0, Ea*CC0)
		return EmatECC

	def EmatELL(self):
		Ea = self.Emat()
		LL0 = self.parameters["longitudinalvector"]
		EmatELL = inner(LL0, Ea*LL0)
		return EmatELL

	def V0constrainedE(self,inverse=False):


		mesh = self.parameters["mesh"]
		u = self.parameters["displacement_variable"]
		ds = dolfin.ds(subdomain_data=self.parameters["facetboundaries"])
		dsendo = ds(self.parameters["endoid"], domain = self.parameters["mesh"], subdomain_data = self.parameters["facetboundaries"])
		pendo = self.parameters["volconst_variable"] 
		V0= self.parameters["constrained_vol"] 

		X = SpatialCoordinate(mesh)
		x = u + X


		F = self.Fmat()
		N = self.parameters["facet_normal"]
		n = cofac(F)*N

		area = assemble(Constant(1.0) * dsendo)
		if inverse:
			V_u = - Constant(0.0/3.0) * inner(x, n)
			Wvol = (Constant(0.0/area) * pendo  * V0 * dsendo) - (pendo * V_u *dsendo)
		else:
		    V_u = - Constant(1.0/3.0) * inner(x, n)
		    Wvol = (Constant(1.0/area) * pendo  * V0 * dsendo) - (pendo * V_u *dsendo)
		return Wvol

    
	def Cauchy1(self): #Joy added this
	  
		u = self.parameters["displacement_variable"]
	   
		d = u.geometric_dimension()
		I = Identity(d)
		F = I + grad(u)
	   
		f0 = self.parameters["fiber"]
		s0 = self.parameters["sheet"]
		n0 = self.parameters["sheet-normal"]
		bff = 29.9 #self.parameters["bff"]
		bfx = 13.3#self.parameters["bfx"]
		bxx = 26.6 #self.parameters["bxx"]
		#Kappa = self.parameters["Kappa"]
		#isincomp = self.parameters["incompressible"]

	  
		p = self.parameters["pressure_variable"]

		C = 200 #self.parameters["C_param"]

		F = dolfin.variable(F)
		J = det(F)

		Ea = 0.5*(as_tensor(F[k,i]*F[k,j] - I[i,j], (i,j)))

		Eff = f0[i]*Ea[i,j]*f0[j]
		Ess = s0[i]*Ea[i,j]*s0[j]
		Enn = n0[i]*Ea[i,j]*n0[j]
		Efs = f0[i]*Ea[i,j]*s0[j]
		Efn = f0[i]*Ea[i,j]*n0[j]
		Ens = n0[i]*Ea[i,j]*s0[j]
	
		
		
		#QQ = bff*pow(Eff,2.0) + bfx*(pow(Ess,2.0)+ pow(Enn,2.0)+ 2.0*pow(Ens,2.0)) + bxx*(2.0*pow(Efs,2.0) + 2.0*pow(Efn,2.0))
		QQ = bff*Eff**2.0 + bfx*(Ess**2.0 + Enn**2.0 + 2.0*Ens**2.0) + bxx*(2.0*Efs**2.0 + 2.0*Efn**2.0)
		Wp = C/2.0*(exp(QQ) -  1.0) - p*(J - 1.0)
		#E = dolfin.variable(Ea)
		sigma =  (1.0/J)*dolfin.diff(Wp,F)*F.T

		
	   
		return sigma


#	   def Cauchy(self):
#	   
#		   d = u.geometric_dimension()
#		   I = Identity(d)
#		   F = I + grad(u)
#		   J = det(F)
#		   Stensor = Stress(u)
#		   sigma = (1.0/J)*F*Stensor*F.T - p*I
#	   
#		   return sigma
#	   
#	   def Stress(u):
#	   
#		   E = Emat(u)
#		  # eQ = exp(E[0,0]*E[0,0] + E[1,1]*E[1,1] + E[2,2]*E[2,2] +  2*(E[0,1]*E[0,1] +  E[0,2]*E[0,2]  +  E[1,2]*E[1,2]))
#		   E11 = E[0,0]
#		   E22 = E[1,1]
#		   E33 = E[2,2]
#		   E12 = E[0,1]
#		   E13 = E[0,2]
#		   E23 = E[1,2]
#		   p1 = 29.9
#		   p2 = 13.3
#		   p3 = 13.5
#		   eQ = exp(29.9*E11*E11 + 13.3*E22*E22 + 13.3*E33*E33 +  2.0*(13.5*E12*E12 + 13.5*E13*E13  +  13.3*E23*E23))
#		   S11 = 2.0*p1*E11*eQ
#		   S22 = 2.0*p2*E22*eQ
#		   S33 = 2.0*p2*E33*eQ
#		   S12 = 2.0*p3*E12*eQ
#		   S13 = 2.0*p3*E13*eQ
#		   S23 = 2.0*p2*E23*eQ
#	   #	S11 = 2*E[0,0]*eQ
#	   #	S22 = 2*E[1,1]*eQ
#	   #	S33 = 2*E[2,2]*eQ
#	   #	S12 = 2*E[0,1]*eQ
#	   #	S13 = 2*E[0,2]*eQ
#	   #	S23 = 2*E[1,2]*eQ
#	   
#		   Smatout = as_matrix([[S11, S12, S13], [S12, S22, S23], [S13, S23, S33]])
#		   return  Smatout





