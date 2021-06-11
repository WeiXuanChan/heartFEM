# A two dimesional finite element implementation of the Nash-Panfilov coupled
# electro-mechanics model in fenics as per their 2004 paper:
#
# "Electromechanical model of excitable tissue to study reentrant cardiac
#   arrythmias", Nash, MP & Panfilov, AV, Progress in Biophysics and Molecular
#   Biology, 85(4):501-522:2004
#
from fenics import *
import math
import os

#
# Simulation and output control
#
# Simulation duration and time step control 
end_time = 160.0  # dimensionless
dt = 0.1          # dimensionless
num_steps = math.ceil(end_time/dt)
# Location & names for solution output files
output_folder_name = "./NashPanfilov_Dynamic_Strong2/"
output_file_name = "solution_phi.pvd"
# How often to write out the solution, in # steps
output_frequency = 10
# This sets up the output files that we will write the solution to. It will
# write each mesh subdomain to a different file if you use mpirun. In paraview
# just load the *.pvd file in either case
vtkfile_phi  = File(output_folder_name + output_file_name)
# Remove the existing solution files, if any
os.system("rm "+ output_folder_name + "*.pvd")
os.system("rm "+ output_folder_name + "*.vtu")
os.system("rm "+ output_folder_name + "*.pvtu")

import matlab.engine
eng = matlab.engine.start_matlab()
tf = eng.isprime(37)
print(tf)

# 
# Create a 2D rectangular mesh
# 
# Note that the crossed option means we will end up with a resolution 
# twice as high as we specify, so 50x50 -> 100x100 nodes
mesh = RectangleMesh(Point(0,0), Point(150,150), 50, 50, 'crossed')
dimension = 2

# 
# Set up the electrics problem
# 
# Linear Lagrange triangles
P1 = FiniteElement('P', triangle, 1)
# Three unknowns (phi, r, ta) so three spaces
element = MixedElement([P1,P1,P1])
V = FunctionSpace(mesh, element)
# Integral weighting functions for finite elements
v1, v2, v3 = TestFunctions(V)
# Create representations of our three unknowns at the new time step...
u = Function(V)
phi, r, Ta  = split(u)
# ...and at the old time step
u_n = Function(V)
phi_n, r_n, Ta_n = split(u_n)
# Directions for Jacobian calculations
du = TrialFunction(V)

# Set all of the initial conditions to be zero over the mesh
zero_at_time_zero = Constant(0.0)
phi_0 = interpolate(zero_at_time_zero, V.sub(0).collapse())
r_0 = interpolate(zero_at_time_zero, V.sub(1).collapse())
Ta_0 = interpolate(zero_at_time_zero, V.sub(2).collapse())
assign(u_n, [phi_0, r_0, Ta_0])

# Here we are assuming no flux all around the boundary of the electrical 
# problem so there is no need to set any explicit boundary conditions as
# zero flux is the default

# The diffusion tensor from Equation 22(a) of Nash & Panfilov (2004)
D = 1.0
D1 = D*Identity(dimension)
# The threshold potential (relative to rest being 0 and plateau 1), Eqn 22(a)
alpha = Constant(0.01)
# Scale factor for the ionic currents, Equation 22(a)  
c = Constant(8.0)
# The next four constants are recovery parameters from Equation 22(b)
g = Constant(0.002)
b = Constant(0.15)
mu1 = Constant(0.2) 
mu2 = Constant(0.3)

# S1-S2 current stimulus protocol for a plane wave and a secondary stimulus
# during the recovery period at the back of the plane wave to cause re-entry
i_s1 = 0  # initial condition
stimulus_1 = Expression(
    'x[0] <= 1 ? i_s1 : 0.0', degree = 0, i_s1 = i_s1)
i_s2 = 0  # initial condition
stimulus_2 = Expression(
    '(x[0] >= 100 && x[0] <= 102) && (x[1] >= 0 and x[1] <= 75) ? i_s2 : 0.0', 
    degree = 0, i_s2 = i_s2)

# Parameters from the that link the voltage to the active stress
kTa = Constant(47.9)
# Equation 23 from Nash & Panfilov (2004)
e_phi = conditional(ge(phi, 0.05), 1.0, 10.0)

# 
# Set up the mechanics problem
# 
# You can set order to 1 for linear Lagrange or 2 for quadratic
order = 1
W = VectorFunctionSpace(mesh,'Lagrange', order) 
# Integral weighting function for finite elements
v4 = TestFunction(W)
# Create a variable to store our displacements...
w = Function(W)
# ...and another to get back to our reference configuration
minus_w = Function(W)
# Directions for Jacobian calculations
dw = TrialFunction(W)

# Boundary conditions. 
# Here we will fix the upper and lower surfaces so they cannot deform
lower =  CompiledSubDomain("near(x[1], side) && on_boundary", side = 0.0)
upper = CompiledSubDomain("near(x[1], side) && on_boundary", side = 150.0) 
zero_two_d = Constant((0.0, 0.0))
# Create the the boundary conditions on W
bc_lower = DirichletBC(W, zero_two_d, lower)
bc_upper = DirichletBC(W, zero_two_d, upper)
bcs_mech = [bc_lower, bc_upper]
# In this case we don't have any surface tractions (T), but they are included
# in the formulation so we will set them to be zero everywhere
T = Constant((0.0, 0.0))
# The formulation also includes a provision for body forces (B), but in this
# case we are not considering them so they are also zero
B = Constant((0.0, 0.0))

# Tensors and invariants for the finite elasticity formulation
I = Identity(dimension)
F = I + grad(w)    # deformation tensor
C = F.T*F          # right Cauchy-Green deformation tensor
Ic = tr(C)         # first invariant
Ic2 = Constant(0.5)*(Ic*Ic - tr(C*C))  # second invariant (for Mooney-Rivlin)
J = det(F)         # Jacobian

# Pull-back of the conductivity tensor to the reference coordinates
F_inv = inv(F)     # F^-1
F_inv_t = F_inv.T  # F^-T
D1_ref = F_inv_t*D1*F_inv
# Note: the Jacobian can be cancelled from all terms of the variational form,
# and so does not appear in D1_ref or the final formulation

# Active stress as a tensor
T_acti = Ta_n*Identity(dimension)

# Strain energy function (Compressible Neo Hookean - 2D)
E = 10.0
nu = 0.3
mu = Constant(E/(2.0*(1.0 + nu)))
lmbda = Constant(E*nu/((1.0 + nu)*(1.0 - nu)))
psi = (mu/2)*(Ic - 2) - mu*ln(J) + (lmbda/2)*(ln(J))**2

# Strain energy function (Compressible Neo Hookean - 3D)
#E = 10.0
#nu = 0.3
#mu = Constant(E/(2.0*(1.0 + nu)))
#lmbda = Constant(E*nu/((1.0 + nu)*(1.0 - 2.0*nu)))
#psi = (mu/2)*(Ic - 3) - mu*ln(J) + (lmbda/2)*(ln(J))**2

# Strain energy function (Mooney-Rivlin - 3D)
#c1 = Constant(2.0)
#c2 = Constant(6.0)
#psi = c1*(Ic - 3) + c2*(Ic2 - 3)

#
# Potential energy minimisation formulation for the mechanics part:
#
# For further reading, this is essentially the same as the code listing in 
# Figure 26.6 in Chapter 26 of the Fenics Book - applications in solid 
# mechanics. 
#
# Potential energy = strain energy - work done by external forces
# The first second is for a body force (B); the third for surface tractions (T)
Pi = psi*dx - dot(B, w)*dx - dot(T, w)*ds
# Minimise potential energy by taking the derivative (here) in each of the 
# directions (v4) and setting them all to zero (later in the time loop)
F_elas = derivative(Pi, w, v4)
# MB: I wonder if this is the correct form for the energy of active stress
# MB: I also wonder what the 0.05 is doing here
F_acti = Constant(0.05)*inner(T_acti, grad(v4))*dx
# The final variational form as the sum of the passive (F_elas) and active 
# (F_acti) components
F_mech = F_elas + F_acti
# If you want to pass in the Jacobian, you can calculate it like this, 
# otherwise you can let Fenics calculate it at each iteration
J_mech = derivative(F_mech, w, dw)

#
# Variational formulation for the electrics part:
#
# For further reading, this is similar to Section 3.5 from the Fenics Tutorial
# volume I - A system of advection-diffusion-reaction equations.
#
# F_v1 is Equation 22(a) from Nash & Panfilov (2004) where their V = phi and
# the time dependent term has been discretised using forward Euler, and Cm = 1
F_v1 = ((phi - phi_n)/Constant(dt))*v1*dx \
    + dot(D1_ref*grad(phi), grad(v1))*dx \
    + c*phi*(phi - alpha)*(phi - Constant(1.0))*v1*dx + r*phi*v1*dx \
    - stimulus_1*v1*dx - stimulus_2*v1*dx
# F_v2 is Equation 22(b) from Nash & Panfilov (2004) - the recovery variable 
F_v2 = ((r - r_n)/Constant(dt))*v2*dx \
    - (g + (mu1*r)/(mu2 + phi))*(-r - c*phi*(phi - b - Constant(1.0)))*v2*dx
# F_v3 is Equation 22(c) from Nash & Panfilov (2004) - the active stress
F_v3 = ((Ta - Ta_n)/Constant(dt))*v3*dx - e_phi*(kTa*phi - Ta)*v3*dx
# The final variational formulation is the sum of these three spaces
F_elec = F_v1 + F_v2 + F_v3
# If you want to pass in the Jacobian, you can calculate it like this, 
# otherwise you can let Fenics calculate it at each iteration
J_elec = derivative(F_elec, u, du)

# 
# Main time loop
#
for n_step in range(num_steps):
    t = n_step*dt

    # The s1 stimulus to create a plane wave at the start
    if t <= 0.5:
        stimulus_1.i_s1 = 5.0
    else:
        stimulus_1.i_s1 = 0.0   
    
    # The s2 stimulus at the back of the action potential to generate re-entry
    if t >= 69.8 and t <= 70.3:
        stimulus_2.i_s2 = 5.0
    else:
        stimulus_2.i_s2 = 0.0   
    
    # Solve the electrics problem
    print('Solving Electrics: t=%0.2f' %t)
    solve(F_elec == 0, u, J = J_elec)
    # If you want to calculate J_elec at every step, you can use this instead
    #solve(F_elec == 0, u)

    # Solve the mechanics problem
    print('Solving Mechanics: t=%0.2f' %t)
    solve(F_mech == 0, w, bcs_mech, J = J_mech)    
 
    # Copy the new u (phi, r, ta) over the old u. Note that we don't need to
    # copy w because it is not time dependent
    u_n.assign(u)  
    
    # Output the solution to a file (or files)
    if abs(n_step % output_frequency) == 0:
        # Deform the mesh using the calculated displacements
        ALE.move(mesh, w)
        # Write the solution to the file (membrane potential, deformed mesh)
        phi_, r_dummy, Ta_dummy = u.split()
        vtkfile_phi << (phi_, t)
        # Reinstate the original mesh (as all calculations are relative to the 
        # reference coordinate system)
        minus_w.vector()[:] = -1.0*w.vector()
        ALE.move(mesh, minus_w) 
