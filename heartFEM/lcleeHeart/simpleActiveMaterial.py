# Guccione et al. 1993: Simple time varying elastance type active contraction model

from dolfin import *
from ufl import i
from ufl import j
import math

class SimpleActiveMaterial(object):

    def __init__(self, params):
        self.parameters = params
        self.default_parameters={"Ca0":4.35,"Ca0max":4.35,"B":4.75,"t0" : 200.5,"l0" : 1.58,"lr" : 1.85,"m" : 1048*0.5,	"b" : -1600*0.5}
    def Fmat(self):
        u = self.parameters["displacement_variable"]
        d = u.geometric_dimension()
        I = Identity(d)
        F = I + grad(u)
        return F

    def set_default_parameters(self,editParameters=None):
        if editParameters is None:
            return
        else:
            for key in self.default_parameters:
                if key in editParameters:
                    self.default_parameters[key]=editParameters[key]

    def PK1StressTensor(self):
        F = self.Fmat()
        f0 = self.parameters["fiber"]
        Mij = f0[i]*f0[j]
        Pact = self.PK1Stress()
        Pact_tensor = Pact*as_tensor(Mij, (i,j)) 
        return Pact_tensor
    def PK1Stress(self):
        Ca0 = self.default_parameters["Ca0"]
        Tmax = self.parameters["T0"]
        Ct = self.Ct()
        ECa = self.ECa()
        Pact = (Tmax*Ct*Ca0**2.0)/(Ca0**2.0 + ECa**2.0)
        return Pact

    def ECa(self):
        F = self.Fmat()
        f0 = self.parameters["fiber"]
        Ca0max = self.default_parameters["Ca0max"]
        B = self.default_parameters["B"]
        l0 = self.default_parameters["l0"]
        t_a = self.parameters["t_a"]
        lr = self.default_parameters["lr"]

        Cmat = F.T*F
        lmbda = sqrt(dot(f0, Cmat*f0))
        ls = lmbda*lr
        ls_l0 = conditional(le(ls, l0+0.003), 0.003, ls - l0)
        #	ls_l0 = conditional(le(ls, l0), 0.1*(l0 - ls), ls - l0)
        denom = sqrt(exp(B*(ls_l0)) - 1)

        ECa = Ca0max/denom
        return ECa

    def tr(self):
        F = self.Fmat()
        f0 = self.parameters["fiber"]
        b = self.default_parameters["b"]
        m = self.default_parameters["m"]
        l0 = self.default_parameters["l0"]
        lr = self.default_parameters["lr"]
        Cmat = F.T*F
        lmbda = sqrt(dot(f0, Cmat*f0))
        ls = lmbda*lr
        
        tr = m*ls + b
        
        return tr

    def Ct(self):
        
        F = self.Fmat()
        f0 = self.parameters["fiber"]
        t0 = self.default_parameters["t0"]
        b = self.default_parameters["b"]
        m = self.default_parameters["m"]
        t_a = self.parameters["t_a"]
        lr = self.default_parameters["lr"]
        
        Cmat = F.T*F
        lmbda = sqrt(dot(f0, Cmat*f0))
        ls = lmbda*lr
        
        tr = self.tr()
        xp1 = conditional(lt(t_a,t0), 1.0, 0.0)
        w1 = xp1*pi*t_a/t0
        xp2 = conditional(le(t0,t_a), 1.0, 0.0)
        xp3 = conditional(lt(t_a,t0+tr), 1.0, 0.0)
        w2 = xp2*xp3*pi*(t_a - t0 + tr)/tr
        
        Ct = 0.5*(1 - cos(w1+w2))
        
        return Ct
    def cauchy(self): #Joy added this
        F = self.Fmat()

        Pact = self.PK1StressTensor() 
        J = det(F)
	
        sigma = (1.0/J)*Pact*F.T
        return sigma


    def CalculateFiberStress(self, sigma, e_fiber, Vol, Mesh): #Joy added this

        Fmat = self.Fmat()

        n = (Fmat*e_fiber)/ sqrt(inner(Fmat*e_fiber, Fmat*e_fiber))
        sigma_fiber = inner(sigma*n, n)
        d_x = dx(Mesh)
        sigma_fiber_LV = 0.0 
        sigma_fiber_LV = assemble(sigma_fiber * d_x,form_compiler_parameters={"representation":"uflacs"})/Vol
        print( sigma_fiber_LV)

        return sigma_fiber_LV
