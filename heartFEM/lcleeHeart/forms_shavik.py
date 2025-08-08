from dolfin import *
from ufl import cofac
import dolfin as dolfin
import numpy as np
from ufl import k
from ufl import i
from ufl import j
#from mpi4py import MPI					 
import sys
import logging
logger = logging.getLogger('forms_shavik')
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
		Cstrain=self.parameters["StrainEnergyDensityFunction_Coef"]
		if self.parameters["cstrainfactor"] is not None:
		    Cstrain=Cstrain*self.parameters["cstrainfactor"].sub(0)
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
		    Wp =(Cstrain*(exp(QQ) -  1.0) - p*(self.J() - 1.0))*dx(self.parameters["mesh"]) #original
		return Wp

	def strainEnergy(self, integrate=False):
	    Ea = self.Emat()
	    f0 = self.parameters["fiber"]
	    s0 = self.parameters["sheet"]
	    n0 = self.parameters["sheet-normal"]
	    p = self.parameters["pressure_variable"]
        
	    Cstrain=self.parameters["StrainEnergyDensityFunction_Coef"]
	    if self.parameters["cstrainfactor"] is not None:
		    Cstrain=Cstrain*self.parameters["cstrainfactor"].sub(0)
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
		    Wp = assemble((Cstrain*(exp(QQ) -  1.0))*dx(self.parameters["mesh"]),form_compiler_parameters={"representation":"uflacs"})
	    else:
		    Wp = Cstrain*(exp(QQ) -  1.0)
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
		Cff=self.parameters["StrainEnergyDensityFunction_Cff"]
		Css=self.parameters["StrainEnergyDensityFunction_Css"]
		Cnn=self.parameters["StrainEnergyDensityFunction_Cnn"]
		Cns=self.parameters["StrainEnergyDensityFunction_Cns"]
		Cfs=self.parameters["StrainEnergyDensityFunction_Cfs"]
		Cfn=self.parameters["StrainEnergyDensityFunction_Cfn"]
		#Kappa = self.parameters["Kappa"]
		#isincomp = self.parameters["incompressible"]

	  
		p = self.parameters["pressure_variable"]

		Cstrain=self.parameters["StrainEnergyDensityFunction_Coef"]
		if self.parameters["cstrainfactor"] is not None:
		    Cstrain=Cstrain*self.parameters["cstrainfactor"].sub(0)

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
		QQ = Cff*Eff**2.0 + Css*Ess**2.0 + Cnn*Enn**2.0 + Cns**Ens**2.0 + Cfs**Efs**2.0 + Cfn*Efn**2.0
		Wp = Cstrain*(exp(QQ) -  1.0) - p*(J - 1.0)
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


class Forms_AV(Forms):
    def __init__(self, params):
        super().__init__(params)
    # def V0constrainedE(self,inverse=False,addstr=''):
    #     mesh = self.parameters["mesh"]
    #     u = self.parameters["displacement_variable"]
    #     ds = dolfin.ds(subdomain_data=self.parameters["facetboundaries"])
    #     dsendo = ds(self.parameters["endo"+addstr+"id"], domain = self.parameters["mesh"], subdomain_data = self.parameters["facetboundaries"])
    #     dsendo2 = ds(self.parameters["endo"+addstr+"AVJid"], domain = self.parameters["mesh"], subdomain_data = self.parameters["facetboundaries"])
        
    #     pendo = self.parameters["vol"+addstr+"const_variable"] 
    #     V0= self.parameters["constrained_vol"+addstr] 

    #     X = SpatialCoordinate(mesh)
    #     x = u + X


    #     F = self.Fmat()
    #     N = self.parameters["facet_normal"]
    #     n = cofac(F)*N
        
    #     area = assemble(Constant(1.0) * dsendo+Constant(1.0) * dsendo2)
    #     if inverse:
    #         V_u = - Constant(0.0/3.0) * inner(x, n)
    #         Wvol = (Constant(0.0/area) * pendo  * V0 * dsendo) - (pendo * V_u *dsendo)
    #     else:
    #         logging.debug("constrained_vol"+addstr+", V0 : "+repr(V0)+', cavityvol : '+repr(self.cavityvol(addstr)))
    #         #V_u = - Constant(1.0/3.0) * inner(x, n)##trying same formulation as self.cavityvol !!debug
    #         V_u = - Constant(1.0/3.0) *inner(det(F)*dot(inv(F).T, N),x)
    #         Wvol = (Constant(1.0/area) * pendo  * V0 * dsendo)+(Constant(1.0/area) * pendo  * V0 * dsendo2) - (pendo * V_u *dsendo) - (pendo * V_u *dsendo2)
    #     return Wvol
    def V0constrainedE(self,inverse=False):
        mesh = self.parameters["mesh"]
        u = self.parameters["displacement_variable"]
        ds = dolfin.ds(subdomain_data=self.parameters["facetboundaries"])
        dsAtrendo = ds(self.parameters["endoAtrid"], domain = self.parameters["mesh"], subdomain_data = self.parameters["facetboundaries"])
        dsAtrendo2 = ds(self.parameters["endoAtrAVJid"], domain = self.parameters["mesh"], subdomain_data = self.parameters["facetboundaries"])
        dsVenendo = ds(self.parameters["endoVenid"], domain = self.parameters["mesh"], subdomain_data = self.parameters["facetboundaries"])
        dsVenendo2 = ds(self.parameters["endoVenAVJid"], domain = self.parameters["mesh"], subdomain_data = self.parameters["facetboundaries"])
        
        pendoAtr = self.parameters["volAtrconst_variable"] 
        V0Atr= self.parameters["constrained_volAtr"] 
        pendoVen = self.parameters["volVenconst_variable"] 
        V0Ven= self.parameters["constrained_volVen"] 
        
        AVJ_control= self.parameters["AVJ_control"]
        X = SpatialCoordinate(mesh)
        x = u + X
        

        F = self.Fmat()
        N = self.parameters["facet_normal"]
        n = cofac(F)*N
        
        
        if inverse:
            V_u = - Constant(0.0/3.0) * inner(x, n)
            Wvol = (Constant(0.0/area) * pendo  * V0 * dsendo) - (pendo * V_u *dsendo)
        else:
            
            logging.debug("constrained_volAtr, V0 : "+repr(V0Atr)+', cavityvol : '+repr(self.cavityvol('Atr')))
            logging.debug("constrained_volVen, V0 : "+repr(V0Ven)+', cavityvol : '+repr(self.cavityvol('Ven')))
            
            ##trying same formulation as self.cavityvol !!debug
            #V_u = - Constant(1.0/3.0) *inner(det(F)*dot(inv(F).T, N),x)
            if AVJ_control==3:
                pendoAVJ=self.parameters["volAVJconst_variable"]
                V0AVJ= self.parameters["constrained_volAVJ"] 
                c11= self.parameters["c11"] 
                V_u = - Constant(1.0/3.0) * inner(x, n)
                areaAtr = assemble(Constant(1.0) * dsAtrendo)
                areaVen = assemble(Constant(1.0) * dsVenendo)
                areaAVJ_Atr = assemble(Constant(1.0) * dsAtrendo2)
                areaAVJ_Ven = assemble(Constant(1.0) * dsVenendo2)
                V_u_atr = - Constant(1.0/3.0) * (inner(x, n)-Constant(2.0/areaAVJ_Atr)*(c11[5]*n[0]+c11[6]*n[1]+c11[7]*n[2]))
                V_u_ven = - Constant(1.0/3.0) * (inner(x, n)-Constant(2.0/areaAVJ_Ven)*(c11[8]*n[0]+c11[9]*n[1]+c11[10]*n[2]))
                V_u_AVJatr = - Constant(1.0/2.0) * (inner(x, n)-Constant(1.0/areaAVJ_Atr)*(c11[5]*n[0]+c11[6]*n[1]+c11[7]*n[2]))
                V_u_AVJven = - Constant(1.0/2.0) * (inner(x, n)-Constant(1.0/areaAVJ_Ven)*(c11[8]*n[0]+c11[9]*n[1]+c11[10]*n[2]))
                WvolAtr = (Constant(1.0/areaAtr) * pendoAtr  * V0Atr * dsAtrendo) - (pendoAtr * V_u_atr *dsAtrendo)
                WvolVen = (Constant(1.0/areaVen) * pendoVen  * V0Ven * dsVenendo) - (pendoVen * V_u_ven *dsVenendo) 
                WvolAVJ = (Constant(1.0/(areaAVJ_Atr+areaAVJ_Ven)) * pendoAVJ  * V0AVJ * dsAtrendo2)+(Constant(1.0/(areaAVJ_Atr+areaAVJ_Ven)) * pendoAVJ  * V0AVJ * dsVenendo2) - (pendoAVJ * V_u_AVJatr *dsAtrendo2) - (pendoAVJ * V_u_AVJven *dsVenendo2)
                Wvol=WvolAtr+WvolVen+WvolAVJ
            elif AVJ_control==2:
                pendoAVJ=self.parameters["volAVJconst_variable"]
                V0AVJ= self.parameters["constrained_volAVJ"] 
                
                V_u = - Constant(1.0/3.0) * inner(x, n)
                areaAtr = assemble(Constant(1.0) * dsAtrendo)
                areaVen = assemble(Constant(1.0) * dsVenendo)
                areaAVJ_Atr = assemble(Constant(1.0) * dsAtrendo2)
                areaAVJ_Ven = assemble(Constant(1.0) * dsVenendo2)
                AVJ_atr_centroid0=Constant(1.0/areaAVJ_Atr)*assemble(x[0]*dsAtrendo2)
                AVJ_atr_centroid1=Constant(1.0/areaAVJ_Atr)*assemble(x[1]*dsAtrendo2)
                AVJ_atr_centroid2=Constant(1.0/areaAVJ_Atr)*assemble(x[2]*dsAtrendo2)
                AVJ_ven_centroid0=Constant(1.0/areaAVJ_Ven)*assemble(x[0]*dsVenendo2)
                AVJ_ven_centroid1=Constant(1.0/areaAVJ_Ven)*assemble(x[1]*dsVenendo2)
                AVJ_ven_centroid2=Constant(1.0/areaAVJ_Ven)*assemble(x[2]*dsVenendo2)
                V_u_atr = - Constant(1.0/3.0) * (inner(x, n)-Constant(2.0)*(AVJ_atr_centroid0*n[0]+AVJ_atr_centroid1*n[1]+AVJ_atr_centroid2*n[2]))
                V_u_ven = - Constant(1.0/3.0) * (inner(x, n)-Constant(2.0)*(AVJ_ven_centroid0*n[0]+AVJ_ven_centroid1*n[1]+AVJ_ven_centroid2*n[2]))
                V_u_AVJatr = - Constant(1.0/2.0) * (inner(x, n)-(AVJ_atr_centroid0*n[0]+AVJ_atr_centroid1*n[1]+AVJ_atr_centroid2*n[2]))
                V_u_AVJven = - Constant(1.0/2.0) * (inner(x, n)-(AVJ_ven_centroid0*n[0]+AVJ_ven_centroid1*n[1]+AVJ_ven_centroid2*n[2]))
                WvolAtr = (Constant(1.0/areaAtr) * pendoAtr  * V0Atr * dsAtrendo) - (pendoAtr * V_u_atr *dsAtrendo)
                WvolVen = (Constant(1.0/areaVen) * pendoVen  * V0Ven * dsVenendo) - (pendoVen * V_u_ven *dsVenendo) 
                WvolAVJ = (Constant(1.0/(areaAVJ_Atr+areaAVJ_Ven)) * pendoAVJ  * V0AVJ * dsAtrendo2)+(Constant(1.0/(areaAVJ_Atr+areaAVJ_Ven)) * pendoAVJ  * V0AVJ * dsVenendo2) - (pendoAVJ * V_u_AVJatr *dsAtrendo2) - (pendoAVJ * V_u_AVJven *dsVenendo2)
                Wvol=WvolAtr+WvolVen+WvolAVJ
            elif AVJ_control:
                pendoAVJ=self.parameters["volAVJconst_variable"]
                V0AVJ= self.parameters["constrained_volAVJ"] 
                V_u = - Constant(1.0/3.0) * inner(x, n)
                areaAtr = assemble(Constant(1.0) * dsAtrendo)
                areaVen = assemble(Constant(1.0) * dsVenendo)
                areaAVJ_Atr = assemble(Constant(1.0) * dsAtrendo2)
                areaAVJ_Ven = assemble(Constant(1.0) * dsVenendo2)
                WvolAtr = (Constant(1.0/areaAtr) * pendoAtr  * V0Atr * dsAtrendo) - (pendoAtr * V_u *dsAtrendo)
                WvolVen = (Constant(1.0/areaVen) * pendoVen  * V0Ven * dsVenendo) - (pendoVen * V_u *dsVenendo) 
                WvolAVJ = (Constant(1.0/(areaAVJ_Atr+areaAVJ_Ven)) * pendoAVJ  * V0AVJ * dsAtrendo2)+(Constant(1.0/(areaAVJ_Atr+areaAVJ_Ven)) * pendoAVJ  * V0AVJ * dsVenendo2) - (pendoAVJ * V_u *dsAtrendo2) - (pendoAVJ * V_u *dsVenendo2)
                Wvol=WvolAtr+WvolVen+WvolAVJ
            else:
                V_u = - Constant(1.0/3.0) * inner(x, n)
                areaAtr = assemble(Constant(1.0) * dsAtrendo+Constant(1.0) * dsAtrendo2)
                areaVen = assemble(Constant(1.0) * dsVenendo+Constant(1.0) * dsVenendo2)
                WvolAtr = (Constant(1.0/areaAtr) * pendoAtr  * V0Atr * dsAtrendo)+(Constant(1.0/areaAtr) * pendoAtr  * V0Atr * dsAtrendo2) - (pendoAtr * V_u *dsAtrendo) - (pendoAtr * V_u *dsAtrendo2)
                WvolVen = (Constant(1.0/areaVen) * pendoVen  * V0Ven * dsVenendo)+(Constant(1.0/areaVen) * pendoVen  * V0Ven * dsVenendo2) - (pendoVen * V_u *dsVenendo) - (pendoVen * V_u *dsVenendo2)
                Wvol=WvolAtr+WvolVen
        return Wvol
    def cavityvol(self,addstr=''):
        
        u = self.parameters["displacement_variable"]
        N = self.parameters["facet_normal"]
        mesh = self.parameters["mesh"]
        X = SpatialCoordinate(mesh)
        ds = dolfin.ds(subdomain_data=self.parameters["facetboundaries"])

        F = self.Fmat()
        n=det(F)*dot(inv(F).T, N)
        x=X + u
        if self.parameters["AVJ_control"]==3:
            c11= self.parameters["c11"] 
            areaAVJ = assemble(Constant(1.0) * ds(self.parameters["endo"+addstr+"AVJid"], domain = self.parameters["mesh"], subdomain_data = self.parameters["facetboundaries"]))
            if addstr=='Atr':
                centroid=as_vector([c11[5],c11[6],c11[7]])
            elif addstr=='Ven':
                centroid=as_vector([c11[8],c11[9],c11[10]])
            vol_form = -Constant(1.0/3.0) * (inner(n, x)-Constant(2.0/areaAVJ)*inner(centroid,n))*ds(self.parameters["endo"+addstr+"id"])
        elif self.parameters["AVJ_control"]==2:
            areaAVJ = assemble(Constant(1.0) * ds(self.parameters["endo"+addstr+"AVJid"], domain = self.parameters["mesh"], subdomain_data = self.parameters["facetboundaries"]))
            AVJ_centroid0=Constant(1.0/areaAVJ)*assemble(x[0]*ds(self.parameters["endo"+addstr+"AVJid"], domain = self.parameters["mesh"], subdomain_data = self.parameters["facetboundaries"]))
            AVJ_centroid1=Constant(1.0/areaAVJ)*assemble(x[1]*ds(self.parameters["endo"+addstr+"AVJid"], domain = self.parameters["mesh"], subdomain_data = self.parameters["facetboundaries"]))
            AVJ_centroid2=Constant(1.0/areaAVJ)*assemble(x[2]*ds(self.parameters["endo"+addstr+"AVJid"], domain = self.parameters["mesh"], subdomain_data = self.parameters["facetboundaries"]))
            vol_form = -Constant(1.0/3.0) * (inner(n, x)-Constant(2.0)*(AVJ_centroid0*n[0]+AVJ_centroid1*n[1]+AVJ_centroid2*n[2]))*ds(self.parameters["endo"+addstr+"id"])
        elif self.parameters["AVJ_control"]:
            vol_form = -Constant(1.0/3.0) * inner(n, x)*ds(self.parameters["endo"+addstr+"id"])
        else:
            vol_form =-Constant(1.0/3.0) * inner(n, x)*ds(self.parameters["endo"+addstr+"id"])-Constant(1.0/3.0) * inner(n,x)*ds(self.parameters["endo"+addstr+"AVJid"])
        return assemble(vol_form)
    
    def avjvol(self):
        
        u = self.parameters["displacement_variable"]
        N = self.parameters["facet_normal"]
        mesh = self.parameters["mesh"]
        X = SpatialCoordinate(mesh)
        ds = dolfin.ds(subdomain_data=self.parameters["facetboundaries"])
			
        F = self.Fmat()
        n=det(F)*dot(inv(F).T, N)
        x=X + u
        if self.parameters["AVJ_control"]==3:
            c11= self.parameters["c11"] 
            areaAVJ_Atr = assemble(Constant(1.0) * ds(self.parameters["endoAtrAVJid"], domain = self.parameters["mesh"], subdomain_data = self.parameters["facetboundaries"]))
            areaAVJ_Ven = assemble(Constant(1.0) * ds(self.parameters["endoVenAVJid"], domain = self.parameters["mesh"], subdomain_data = self.parameters["facetboundaries"]))
            vol_form = -Constant(1.0/2.0) * (inner(n, x)-Constant(1.0/areaAVJ_Atr)*(c11[5]*n[0]+c11[6]*n[1]+c11[7]*n[2]))*ds(self.parameters["endoAtrAVJid"])
            vol_form2 = -Constant(1.0/2.0) * (inner(n, x)-Constant(1.0/areaAVJ_Ven)*(c11[8]*n[0]+c11[9]*n[1]+c11[10]*n[2]))*ds(self.parameters["endoVenAVJid"])
        elif self.parameters["AVJ_control"]==2:
            areaAVJ_Atr = assemble(Constant(1.0) * ds(self.parameters["endoAtrAVJid"], domain = self.parameters["mesh"], subdomain_data = self.parameters["facetboundaries"]))
            areaAVJ_Ven = assemble(Constant(1.0) * ds(self.parameters["endoVenAVJid"], domain = self.parameters["mesh"], subdomain_data = self.parameters["facetboundaries"]))
            AVJ_atr_centroid0=Constant(1.0/areaAVJ_Atr)*assemble(x[0]*ds(self.parameters["endoAtrAVJid"], domain = self.parameters["mesh"], subdomain_data = self.parameters["facetboundaries"]))
            AVJ_atr_centroid1=Constant(1.0/areaAVJ_Atr)*assemble(x[1]*ds(self.parameters["endoAtrAVJid"], domain = self.parameters["mesh"], subdomain_data = self.parameters["facetboundaries"]))
            AVJ_atr_centroid2=Constant(1.0/areaAVJ_Atr)*assemble(x[2]*ds(self.parameters["endoAtrAVJid"], domain = self.parameters["mesh"], subdomain_data = self.parameters["facetboundaries"]))
            AVJ_ven_centroid0=Constant(1.0/areaAVJ_Ven)*assemble(x[0]*ds(self.parameters["endoVenAVJid"], domain = self.parameters["mesh"], subdomain_data = self.parameters["facetboundaries"]))
            AVJ_ven_centroid1=Constant(1.0/areaAVJ_Ven)*assemble(x[1]*ds(self.parameters["endoVenAVJid"], domain = self.parameters["mesh"], subdomain_data = self.parameters["facetboundaries"]))
            AVJ_ven_centroid2=Constant(1.0/areaAVJ_Ven)*assemble(x[2]*ds(self.parameters["endoVenAVJid"], domain = self.parameters["mesh"], subdomain_data = self.parameters["facetboundaries"]))
            vol_form = -Constant(1.0/2.0) * (inner(n, x)-(AVJ_atr_centroid0*n[0]+AVJ_atr_centroid1*n[1]+AVJ_atr_centroid2*n[2]))*ds(self.parameters["endoAtrAVJid"])
            vol_form2 = -Constant(1.0/2.0) * (inner(n, x)-(AVJ_ven_centroid0*n[0]+AVJ_ven_centroid1*n[1]+AVJ_ven_centroid2*n[2]))*ds(self.parameters["endoVenAVJid"])
        else:
            vol_form = -Constant(1.0/3.0) * inner(n, x)*ds(self.parameters["endoAtrAVJid"])
            vol_form2 = -Constant(1.0/3.0) * inner(n, x)*ds(self.parameters["endoVenAVJid"])
        return assemble(vol_form+vol_form2)
    
    def cavitypressure(self,addstr=''):
        if addstr=='Atr':
            idnum=2
        elif addstr=='Ven':
            idnum=3
        elif addstr=='AVJ':
            idnum=4
        else:
            raise Exception("Indicate Atr or Ven or AVJ for cavitypressure.")
            
        W = self.parameters["mixedfunctionspace"]
        w = self.parameters["mixedfunction"]
        mesh = self.parameters["mesh"]

        comm = W.mesh().mpi_comm()
        dofmap =  W.sub(idnum).dofmap()
        val_dof = dofmap.cell_dofs(0)[0]

			# the owner of the dof broadcasts the value
        own_range = dofmap.ownership_range()
	
        try:
            val_local = w.vector()[val_dof]				  
        except IndexError:
            val_local = 0.0


        pressure = MPI.sum(comm, val_local)
        return pressure
    def currentFiberVector(self):
        F = self.Fmat()
        f0 = self.parameters["fiber"]
        return dot(F, f0)