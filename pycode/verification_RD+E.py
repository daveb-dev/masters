from dolfin import *
from dolfin_adjoint import *
from numpy import fliplr, linspace
from os.path import join as osjoin
from scipy.io import loadmat as sc_io_loadmat
from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import colors

set_log_level(ERROR) 

def param_update(nonlinear_solver):
    prm1 = nonlinear_solver.parameters
    prm1['newton_solver']['absolute_tolerance'] = 1E-7
    prm1['newton_solver']['relative_tolerance'] = 1E-6
    prm1['newton_solver']['maximum_iterations'] = 51
    prm1['newton_solver']['relaxation_parameter'] = 1.0
    prm1['newton_solver']['linear_solver'] = 'gmres'
    prm1['newton_solver']['preconditioner'] = 'ilu'
    prm1['newton_solver']['krylov_solver']['absolute_tolerance'] = 1E-8
    prm1['newton_solver']['krylov_solver']['relative_tolerance'] = 1E-6
    prm1['newton_solver']['krylov_solver']['maximum_iterations'] = 1000
    prm1['newton_solver']['krylov_solver']['nonzero_initial_guess'] = True
    return None

def forward(initial_p, name, record=False,  annotate=False):
    """ 
        Here, we define the forward problem. 
    """
    '''
        - E(u) returns the Green-Lagrange strain tensor
        - sigma(...) returns the actual stress tensor
        - sigma_form(...) returns the stress tensor based on the cells (phi), elasticity coefficients, and a coefficient beta
        - vonmises(...) calculates the von Mises stress based on the actual stress tensor
    '''
    global t
    
    ## Define functions
    I = Identity(2)  # Identity tensor
    def E(u):
        return 0.5*(nabla_grad(u) + nabla_grad(u).T)
    def vonmises(u):
        s = sigma(u) - (1./2)*tr(sigma(u))*I  # deviatoric stress
        von_Mises = sqrt(3./2*inner(s, s))
        return project(von_Mises, V)

    # Set up problem
    U    = VectorFunctionSpace(mesh,'Lagrange',1)
    def boundary(x, on_boundary):
        return on_boundary
    bc   = DirichletBC(U, Constant((0.,0.)), boundary)
    p_n  = interpolate(initial_p,V)
    v    = TestFunction(U)
    
    ffc_options = {"quadrature_degree": 2, "cpp_optimize": True}
    parameters['form_compiler']['quadrature_degree'] = 2
    parameters['form_compiler']['cpp_optimize'] = True
    parameters['krylov_solver']['nonzero_initial_guess'] = True
    
    if lin_hyp == 0:
        def sigma(u):
            s = 2*mu*E(u)+lmbda*tr(E(u))*I
            return s
        u    = TrialFunction(U)
        a = inner(2*mu*E(u)+lmbda*tr(E(u))*I,E(v))*dx
        L = inner(2*beta*p_n*I*(mu+lmbda),E(v))*dx
        u    = Function(U)
        def mech():
            solve(a == L, u, bc, 
                      form_compiler_parameters=ffc_options)
            return u
    else:
        def sigma(u):
            F = I + grad(u)             # Deformation gradient
            B = F*F.T
            C = F.T*F
            J = det(F)
            I1 = tr(C)
            s = lmbda*(J-1)*I+mu*(B-1./2*I1*I)/(J**(5./3))
            return s
        def sigma_form(u, phi):
            F = I + grad(u)             # Deformation gradient
            Fs = F/(1+beta*phi)
            Bs = Fs*Fs.T
            Js  = det(Fs)
            return 1/(1+beta*phi)*(mu/(Js**(5./3))*(Bs-1./2*tr(Bs)*I)+lmbda*(Js-1)*I)
        u           = Function(U)
        du          = TrialFunction(U)
        F_HE        = inner(sigma_form(u, p_n), E(v))*dx
        J_HE        = derivative(F_HE,u,du)
        problem_HE  = NonlinearVariationalProblem(F_HE, u, bc, J=J_HE,form_compiler_parameters=ffc_options)
        solver_HE   = NonlinearVariationalSolver(problem_HE)
        param_update(solver_HE)
        def mech():
            solver_HE.solve()
            return u

    # First iteration solving for displacement, and using the von mises stress field for D
    disp = mech()
    vm   = vonmises(disp)
    D    = project(D0*exp(-gammaD*vm),V)
    
    # Set up reaction-diffusion problem
    dp   = TrialFunction(V)
    p    = Function(V)
    q    = TestFunction(V)
    F_RD = (1/dt)*(p - p_n)*q*dx + D*dot(grad(q),grad(p))*dx - k*p*(1 - p_n)*q*dx  
    J_RD = derivative(F_RD,p) 

    Problem_RD  = NonlinearVariationalProblem(F_RD, p, J=J_RD,form_compiler_parameters=ffc_options)
    solver_RD   = NonlinearVariationalSolver(problem_RD)
    param_update(solver_RD)
        
    for n in range(num_steps):
        # Update current time
        solver_RD.solve(annotate=annotate)
        p_n.assign(p)
        
        # Update previous solution
        disp = mech()
        vm   = vonmises(disp)
        D    = project(D0*exp(-gammaD*vm),V)
        t += dt
    return p

#########################################################################
# MAIN   
########################################################################
# Prepare a mesh
mesh = Mesh(input_dir+"gmsh.xml")
V    = FunctionSpace(mesh, 'CG', 2)
lin_hyp  = 0
D0       = 1.
gammaD   = 1.
k0       = 1.
beta     = 1.
dt       = .1

# Model parameters
t          = 0.               # initial time
T          = 2.0              # final time 
num_steps  = int(T/dt)        # number of time steps
theta      = 50970.           # carrying capacity - normalize data by this
mu         = .42              # kPa, bulk shear modulus
nu         = .45              # poisson's ratio
lmbda      = 2*mu*nu/(1-2*nu) # lame parameter
k          = project(k0, V)

u = .1*(cos(pi*x/50)+cos(pi*y/50)+e^(.88*t)+2)
R = .088*e^(.88*t)+.1*(pi/50)^2*(cos(pi*x/50)+cos(pi*y/50)) - .1*u*(1-u)
pde = Dt-Dxx-Dyy-k-R = 0
error = pde-u
# Load initial tumor condition data
initial_p = interp(input_dir+"ic.mat","ic")
initial_p.rename('initial','tumor at day 0')

# run the forward model
forward(initial_p, case) 

