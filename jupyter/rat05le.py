from __future__ import print_function
from dolfin import *
from dolfin_adjoint import *
from numpy import fliplr, linspace
from os.path import join as osjoin
from scipy.io import loadmat as sc_io_loadmat
from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import colors
from xdmf_parser import xparse as xp

set_log_level(PROGRESS) 

class InterpolatedParameter(Expression):
    '''
        Class to get tumor cell distributions by interpolating based off matrices of tumor cell data
    '''
    def __init__(self,X,Y,image,**kwargs):
        self.X = X # A numpy array giving the X-spacing of the image
        self.Y = Y # Same for Y
        self.image = image # The image of measured material property    
    def eval_cell(self,values,x,cell):
        interp_handle = RegularGridInterpolator((self.X,self.Y),self.image)
        values[0] = interp_handle(x)
        
def interp(file_loc,mat_name):
    """
        Function to accept matlab .mat file with tumor data and interpolate values onto mesh
    """
    print("Interpolating "+mat_name+"...")
    mat = sc_io_loadmat(file_loc)[mat_name]
    mat = fliplr(mat.T)/theta  # Needs to be adjusted to fit the mesh correctly
    x,y = mat.shape[0], mat.shape[1]
    mat_interp = InterpolatedParameter(linspace(1,x,x),linspace(1,y,y),mat,degree=1)
    return interpolate(mat_interp,V)

def vis_obs(initial_p,target_p,title1,title2):
    '''
        Compare two quantity plots, for example initial vs. target cellularity
        Accepts titles for each plot
    '''
    print("Plotting "+title1+" and "+title2)
    cm1 = cm.get_cmap('jet')
    plt.figure()
    plt.subplot(1,2,1)
    plt.title(title1)
    plot(quantity1,cmap=cm1)
 
    plt.subplot(1,2,2)
    plt.title(title2)
    plot(quantity2,cmap=cm1)
    
    savefig(title1+'_and_'+title2+'.png')
    if take_diff:
        diff = Function(V,annotation=False)
        diff.rename('diff','diff tumor fraction')
        diff.vector()[:] = quantity1.vector()-quantity2.vector()
        file_results.write(diff)

def forward(initial_p, record=False, annotate=False):
    """ 
        Here, we define the forward problem. 
    """
    '''
        - E(u) returns the Green-Lagrange strain tensor
        - sigma(...) returns the actual stress tensor
        - sigma_form(...) returns the stress tensor based on the cells (phi), elasticity coefficients, and a coefficient beta
        - vonmises(...) calculates the von Mises stress based on the actual stress tensor
    '''
    def E(u):
        return 0.5*(nabla_grad(u) + nabla_grad(u).T)
    def sigma(u):
        return 2*mu*E(u)+lmbda*tr(E(u))*Identity(2)
    def vonmises(u):
        s         = sigma(u) - (1./2)*tr(sigma(u))*Identity(2)  # deviatoric stress
        von_Mises = sqrt(3./2*inner(s, s))
        return project(von_Mises, V,annotate=annotate)
    def sigma_form(u,phi):
        return 2*mu*(E(u)-beta*phi*Identity(2))+lmbda*(tr(E(u))-2*beta*phi)*Identity(2)
    
    U           = VectorFunctionSpace(mesh,'Lagrange',1)
    def boundary(x, on_boundary):
        return on_boundary
    bc          = DirichletBC(U, Constant((0.,0.)), boundary)
    p_n         = interpolate(initial_p,V)
    u           = TrialFunction(U)
    v           = TestFunction(U)
    F_LE        = inner(sigma_form(u, p_n), E(v))*dx 
    a, L        = lhs(F_LE), rhs(F_LE)
    u           = Function(U,annotate=annotate)
    parameters["form_compiler"]["quadrature_degree"] = 2
    parameters["form_compiler"]["cpp_optimize"] = True
    parameters['krylov_solver']['nonzero_initial_guess'] = True
    ffc_options = {"quadrature_degree": 2}

    # First iteration solving for displacement, and using the von mises stress field for D
    solve(a == L, u, bc, form_compiler_parameters=ffc_options,annotate=annotate)
    vm          = vonmises(u)
    D           = project(D0*exp(-gammaD*vm),V,annotate=annotate)

    # Set up reaction-diffusion problem
    dp          = TrialFunction(V)
    p           = Function(V,annotate=annotate)
    q           = TestFunction(V)
    F_RD        = (1/dt)*(p - p_n)*q*dx + D*dot(grad(q),grad(p))*dx - k*p*(1 - p)*q*dx  
    J_RD        = derivative(F_RD,p,dp)
    problem_RD  = NonlinearVariationalProblem(F_RD, p, J=J_RD,form_compiler_parameters=ffc_options)
    solver_RD   = NonlinearVariationalSolver(problem_RD)
    solver_RD.parameters['newton_solver']['krylov_solver']['nonzero_initial_guess'] = True

    # Prepare the solution
    t = dt
        
    for n in range(num_steps):
        if (n%rtime == 0):
            print("Solved reaction diffusion for time = "+str(t))
            if record:        # save the current solution, k field, displacement, and diffusion
                u.rename('u','displacement')
                p_n.rename('phi_T','tumor fraction')
                vm.rename("vm","Von Mises")
                D.rename("D","diffusion coefficient")
                k.rename('k','k field')          
                file_results.write(p_n,t)
                file_results.write(u,t)
                file_results.write(k,t)
                file_results.write(vm,t)
                file_results.write(D,t)
        # Update current time
        solver_RD.solve(annotate=annotate)
        p_n.assign(p)
        
        # Update previous solution
        solve(a == L, u, bc, form_compiler_parameters=ffc_options,annotate=annotate)
        vm   = vonmises(u)
        D    = project(D0*exp(-gammaD*vm),V,annotate=annotate)

        t += dt
    return p

# Callback function for the optimizer
# Writes intermediate results to a logfile
def eval_cb(j, m):
    """ The callback function keeping a log """
    print("objective = %15.10e " % j)

def objective(p, target_p, r_coeff1, r_coeff2):
    return assemble(inner(p-target_p, p-target_p)*dx) + r_coeff1*assemble(k*k*dx) + r_coeff2*assemble(dot(grad(k),grad(k))*dx)

def optimize(dbg=False):
    """ The optimization routine """
    print("Optimizing...")
    
    # Define the control
    m = [Control(k), Control(D0)]
    
    # Execute first time to annotate and record the tape
    p = forward(initial_p, True, True)

    J = objective(p, target_p, r_coeff1, r_coeff2)

    # Prepare the reduced functional
    rf = ReducedFunctional(J,m,eval_cb_post=eval_cb)
    
    # upper and lower bound for the parameter field
    k_lb, k_ub = Function(V,annotate=False), Function(V,annotate=False)
    k_lb.vector()[:] = 0.
    k_ub.vector()[:] = 4.
    D_lb = 0.
    D_ub = 4.
    bnds = [[k_lb,D_lb],[k_ub,D_ub]]

    # Run the optimization
    m_opt = minimize(rf,method='L-BFGS-B', bounds=bnds, tol=1.0e-4,options={"disp":True,"gtol":1.0e-4})
    
    return m_opt

#########################################################################
# MAIN 
########################################################################
t1         = time()
case       = 0
r_coeff1   = 0.01
r_coeff2   = 0.01
input_dir  = "../rat-data/rat05/"
output_dir = './output/rat05/le'

# Prepare output file
file_results = XDMFFile(osjoin(output_dir,'le.xdmf'))
file_results.parameters["flush_output"] = True
file_results.parameters["functions_share_mesh"] = True
rtime = 10 # How often to record results

# Prepare a mesh
mesh = Mesh(input_dir+"gmsh.xml")
V    = FunctionSpace(mesh, 'CG', 1)

# Model parameters
T             = 2.               # final time 
num_steps     = 100              # number of time steps
dt            = T/num_steps      # time step size
theta         = 50970.           # carrying capacity - normalize cell data by this 
mu            = .42              # kPa, bulk shear modulus
nu            = .45
lmbda         = 2*mu*nu/(1-2*nu)
beta          = 1.
gammaD        = 2.     # initial guess of gamma_D

# Parameters to be optimized
D0     = Constant(1.)            # mobility or diffusion coefficient
k0     = Constant(2.)            # growth rate initial guess
k      = project(k0,V)

# Load initial tumor condition data
initial_p = interp(input_dir+"ic.mat","ic")
initial_p.rename('initial','tumor at day 0')
# target_p  = interp(input_dir+"tumor_t2.mat","tumor")  
# target_p.rename('target','tumor at day 2')

annotate=False
target_p = forward(initial_p,False,False)

# Visualize initial cellularity and target cellularity
vis_obs(initial_p,target_p,'initial','target') 

# Initial guesses
D0     = Constant(2.)     # mobility or diffusion coefficient
k0     = Constant(1.5)    # growth rate initial guess
k      = project(k0,V)

# Optimization module
[k, D0] = optimize() # optimize the k field, gammaD, and D0 using the adjoint method provided by adjoint_dolfin
print('Elapsed time is ' + str((time()-t1)/60) + ' minutes')

model_p = forward(initial_p,False,False) # run the forward model using the optimized k field
vis_obs(model_p,target_p,'model','actual')

print('J_opt = '+str(objective(model_p, target_p, r_coeff1, r_coeff2)))
print('J_opt (without regularization) = '+str(objective(model_p, target_p, 0., 0.)))
print('D0 = '+str(D0.values()[0]))

xp('/output/rat05/le/le.xdmf')
