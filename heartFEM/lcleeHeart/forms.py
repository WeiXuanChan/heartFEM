from dolfin import *
from ufl import cofac
import dolfin as dolfin
import numpy as np
#from mpi4py import MPI
import sys
import math



class Forms(object):

	def __init__(self, params):

		self.parameters = params

	def Fmat(self):

		u = self.parameters["displacement_variable"]
		d = u.geometric_dimension()
		I = Identity(d)
		F = I + grad(u)
		return F

	def Emat(self):

		u = self.parameters["displacement_variable"]
		d = u.geometric_dimension()
		I = Identity(d)
		F = self.Fmat()
		return 0.5*(F.T*F-I)

	def J(self):
		F = self.Fmat()
		return det(F)

	def cavityvol(self):

		u = self.parameters["displacement_variable"]
		N = self.parameters["facet_normal"]
		mesh = self.parameters["mesh"]
		X = SpatialCoordinate(mesh)
		ds = dolfin.ds(subdomain_data=self.parameters["facetboundaries"])

		F = self.Fmat()

		vol_form = -Constant(1.0/3.0) * inner(det(F) *
											  dot(inv(F).T, N), X + u)*ds(self.parameters["endoid"])

		return assemble(vol_form)

	def cavitypressure(self):

		W = self.parameters["mixedfunctionspace"]
		w = self.parameters["mixedfunction"]
		mesh = self.parameters["mesh"]

		comm = W.mesh().mpi_comm()
		val_dof = dofmap.cell_dofs(0)[0]
#		comm = W.mesh().MPI.comm_world
		dofmap = W.sub(2).dofmap()

		# the owner of the dof broadcasts the value
		own_range = dofmap.ownership_range()

		try:
			val_local = w.vector()[val_dof]
			#val_local = np.array(w.vector())[val_dof][0]
						#	val_local = w.vector()[val_dof][0] original
		except IndexError:
			val_local = 0.0

		pressure = MPI.sum(comm, val_local)

		# Serial ###############################################
		#dofmap =  W.sub(2).dofmap()
		#val_dof = dofmap.cell_dofs(0)[0]
		#pressure = w.vector()[val_dof][0]
		########################################################

		return pressure

	def PassiveMatSEF(self):

		Ea = self.Emat()
		f0 = self.parameters["fiber"]
		s0 = self.parameters["sheet"]
		n0 = self.parameters["sheet-normal"]
		p = self.parameters["pressure_variable"]

		Eff = inner(f0, Ea*f0)
		Ess = inner(s0, Ea*s0)
		Enn = inner(n0, Ea*n0)
		Efs = inner(f0, Ea*s0)
		Efn = inner(f0, Ea*n0)
		Ens = inner(n0, Ea*s0)

		# QQ = 29.9*pow(Eff,2.0) + 13.3*(pow(Ess,2.0)+ pow(Enn,2.0)+ 2.0*pow(Ens,2.0)) + 26.6*(2.0*pow(Efs,2.0) + 2.0*pow(Efn,2.0)) #Original
		# QQ = 9.2*pow(Eff,2.0) + 2.0*(pow(Ess,2.0)+ pow(Enn,2.0)+ 2.0*pow(Ens,2.0)) + 2*3.7*(2.0*pow(Efs,2.0) + 2.0*pow(Efn,2.0)) #Original
		# Wp =(33*(exp(QQ) -  1.0) - p*(self.J() - 1.0))*dx(self.parameters["mesh"]) #original
		QQ = 29.9*pow(Eff, 2.0) + 13.3*(pow(Ess, 2.0) + pow(Enn, 2.0) + 2.0 *
										pow(Ens, 2.0)) + 26.6*(2.0*pow(Efs, 2.0) + 2.0*pow(Efn, 2.0))  # Original
		Wp = (100*(exp(QQ) - 1.0) - p*(self.J() - 1.0)) * \
			dx(self.parameters["mesh"])  # original
		# QQ = 9.2*pow(Eff,2.0) + 2.0*(pow(Ess,2.0)+ pow(Enn,2.0)+ 2.0*pow(Ens,2.0)) + 2*3.7*(2.0*pow(Efs,2.0) + 2.0*pow(Efn,2.0)) #McCuloch
		# Wp =(0.5*330*(exp(QQ) -  1.0) + 350*(self.J()-1)*ln(self.J())/2)*dx(self.parameters["mesh"]) #McCuloch chiwei intrepretation
		#Wp =(0.5*330*(exp(QQ) -  1.0) + 350*(self.J()-1)*ln(self.J())/2)
		#QQ = 9.2*pow(Eff,2.0) + 2.0*(pow(Ess,2.0)+ pow(Enn,2.0)+ 2.0*pow(Ens,2.0)) + 3*3.7*(2.0*pow(Efs,2.0) + 2.0*pow(Efn,2.0))
		#Wp =(330*(exp(QQ) -  1.0) - p*(self.J() - 1.0))*dx(self.parameters["mesh"])

		return Wp

	def V0constrainedE(self):

		mesh = self.parameters["mesh"]
		u = self.parameters["displacement_variable"]
		ds = dolfin.ds(subdomain_data=self.parameters["facetboundaries"])
		dsendo = ds(self.parameters["endoid"], domain=self.parameters["mesh"],
					subdomain_data=self.parameters["facetboundaries"])
		pendo = self.parameters["volconst_variable"]
		V0 = self.parameters["constrained_vol"]

		X = SpatialCoordinate(mesh)
		x = u + X

		F = self.Fmat()
		N = self.parameters["facet_normal"]
		n = cofac(F)*N

		area = assemble(Constant(1.0) * dsendo)
		V_u = - Constant(1.0/3.0) * inner(x, n)
		Wvol = (Constant(1.0/area) * pendo * V0 *
				dsendo) - (pendo * V_u * dsendo)

		return Wvol


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
