from __future__ import print_function
from dolfin import *
from dolfin_adjoint import *
from numpy import fliplr, linspace, inf
from os.path import join as osjoin
from scipy.io import loadmat as sc_io_loadmat
from scipy.interpolate import RegularGridInterpolator

set_log_level(DEBUG) 

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
    f_log.write("Interpolating "+mat_name+"...\n")
    mat = sc_io_loadmat(file_loc)[mat_name]
    mat = fliplr(mat.T)/theta  # Needs to be adjusted to fit the mesh correctly; also scaled
    x,y = mat.shape[0], mat.shape[1]
    mat_interp = InterpolatedParameter(linspace(1,x,x),linspace(1,y,y),mat,degree=1)
    return interpolate(mat_interp,V)

def forward(initial_p, name, record=False, annotate=False):    
    """ 
        Here, we define the forward problem. 
    """
    '''
        - E(u) returns the Green-Lagrange strain tensor
        - sigma(...) returns the actual stress tensor
        - sigma_form(...) returns the stress tensor based on the cells (phi), elasticity coefficients, and a coefficient beta
        - vonmises(...) calculates the von Mises stress based on the actual stress tensor
    '''
    global D
    
    ## Define functions
    def E(u):
        return 0.5*(nabla_grad(u) + nabla_grad(u).T)
    def sigma(u):
        I = Identity(2)             # Identity tensor
        F = I + grad(u)             # Deformation gradient
        B = F*F.T
        C = F.T*F
        J = det(F)
        I1 = tr(C)
        sigma = lmbda*(J-1)*I+mu*(B-1./2*I1*I)/(J**(5./3))
        return sigma
    def vonmises(u):
        s = sigma(u) - (1./2)*tr(sigma(u))*Identity(2)  # deviatoric stress
        von_Mises = sqrt(3./2*inner(s, s))
        return project(von_Mises, V, annotate=annotate)
    def sigma_form(u, phi):
        I = Identity(2)             # Identity tensor
        F = I + grad(u)             # Deformation gradient
        Fs = F/(1+beta*phi)
        Bs = Fs*Fs.T
        Js  = det(Fs)
        return 1/(1+beta*phi)*(mu/(Js**(5./3))*(Bs-1./2*tr(Bs)*I)+lmbda*(Js-1)*I)

    # Set up hyperelasticity problem
    U           = VectorFunctionSpace(mesh,'Lagrange',1)
    def boundary(x, on_boundary):
        return on_boundary
    bc          = DirichletBC(U, Constant((0.,0.)), boundary)
    du          = TrialFunction(U)
    u           = Function(U, annotate=annotate) 
    v           = TestFunction(U)
    p_n         = interpolate(initial_p,V)
    F_HE        = inner(sigma_form(u, p_n), E(v))*dx
    J_HE = derivative(F_HE,u)
    ffc_options = {"quadrature_degree": 2, "cpp_optimize": True}
    problem_HE  = NonlinearVariationalProblem(F_HE, u, bc,J=J_HE, form_compiler_parameters=ffc_options)
    solver_HE   = NonlinearVariationalSolver(problem_HE)
    solver_HE.parameters['newton_solver']['krylov_solver']['nonzero_initial_guess'] = True
    # First iteration solving for displacement, and using the von mises stress field for D
    solver_HE.solve(annotate=annotate)
    vm   = vonmises(u)
    D    = project(D0*exp(-gammaD*vm),V,annotate=annotate)
    
    # Set up reaction-diffusion problem
    dp          = TrialFunction(V)
    p           = Function(V,annotate=annotate)
    q           = TestFunction(V)
    F_RD        = (1/dt)*(p - p_n)*q*dx + D*dot(grad(q),grad(p))*dx - k*p*(1 - p)*q*dx  
    J_RD = derivative(F_RD,p)
    problem_RD  = NonlinearVariationalProblem(F_RD, p,J=J_RD,form_compiler_parameters=ffc_options)
    solver_RD   = NonlinearVariationalSolver(problem_RD)
    solver_RD.parameters['newton_solver']['krylov_solver']['nonzero_initial_guess'] = True
    
    # Prepare the solution
    t = 0.
    if record:        # save the current solution, k field, displacement, and diffusion
        u.rename('u_'+name,'displacement')
        p_n.rename('phi_T_'+name,'tumor fraction')
        vm.rename('vm_'+name,'Von Mises')
        D.rename('D_'+name,'diffusion coefficient')
        k.rename('k_'+name,'k field')          
        f_timeseries.write(p_n,t)
        f_timeseries.write(u,t)
        f_timeseries.write(k,t)
        f_timeseries.write(vm,t)
        f_timeseries.write(D,t)
        
    for n in range(num_steps+1):
        t += dt
        
        # Update current time and Compute solution
        solver_RD.solve(annotate=annotate)
        p_n.assign(p)
        
        # Update previous solution
        solver_HE.solve(annotate=annotate)
        vm   = vonmises(u)
        D    = project(D0*exp(-gammaD*vm),V,annotate=annotate)
        
        if (n%rtime == 0):
            print("Solved reaction diffusion for time = "+str(t))
            if record:        # save the current solution, k field, displacement, and diffusion
                u.rename('u_'+name,'displacement')
                p_n.rename('phi_T_'+name,'tumor fraction')
                vm.rename('vm_'+name,'Von Mises')
                D.rename('D_'+name,'diffusion coefficient')
                k.rename('k_'+name,'k field')  
                f_timeseries.write(p_n,t)
                f_timeseries.write(u,t)
                f_timeseries.write(k,t)
                f_timeseries.write(vm,t)
                f_timeseries.write(D,t)
                
    return p

# Callback function for the optimizer; Writes intermediate results to a logfile
def eval_cb(j, m):
    """ The callback function keeping a log """
    f_log.write("objective = %15.10e \n" % j)

def objective(p, target_p, r_coeff1, r_coeff2):
    return assemble(inner(p-target_p, p-target_p)*dx) + \
        r_coeff1*assemble(k*k*dx) + \
        r_coeff2*assemble(dot(grad(k),grad(k))*dx) + \
        r_coeff1*assemble(D*D*dx) + \
        r_coeff2*assemble(dot(grad(D),grad(D))*dx)

def optimize(dbg=False):
    """ The optimization routine """
    f_log.write("Optimizing...\n")
    
    # Define the control
    m = [Control(k), Control(D0), Control(gammaD), Control(beta)]
    
    # Execute first time to annotate and record the tape
    p = forward(initial_p, 'annt', True, True)

    Obj = objective(p, target_p, r_coeff1, r_coeff2)

    # Prepare the reduced functional
    rf = ReducedFunctional(Obj,m,eval_cb_post=eval_cb)

    # upper and lower bound for the parameter field
    k_lb, k_ub = Function(V,annotate=False), Function(V,annotate=False)
    k_lb.vector()[:] = 0.
    k_ub.vector()[:] = 5.
    D_lb = 0.
    D_ub = 1.
    gD_lb = 0
    gD_ub = 1.
    beta_lb = 0
    beta_ub = 1.
    bnds = [[k_lb,D_lb, gD_lb, beta_lb],[k_ub,D_ub, gD_ub, beta_ub]]
    
    # Run the optimization
    m_opt = minimize(rf,method='L-BFGS-B', bounds=bnds, tol=1.0e-8,options={"disp":True,"gtol":1.0e-8})
        
    return m_opt

#########################################################################
# MAIN 
########################################################################
t1         = time()
case       = 0
r_coeff1   = 0.01
r_coeff2   = 0.01
input_dir  = "../rat-data/rat05/"
output_dir = './output/rat05hecoarse'

# Prepare output file
f_timeseries = XDMFFile(osjoin(output_dir,'timeseries.xdmf'))
f_timeseries.parameters["flush_output"] = True
f_timeseries.parameters["functions_share_mesh"] = True
f_notime     = XDMFFile(osjoin(output_dir,'notime.xdmf'))
f_notime.parameters["flush_output"] = True
f_notime.parameters["functions_share_mesh"] = True
f_log = open(osjoin(output_dir,'info.log'),'w+')
rtime = 5 # How often to record results 

# Prepare a mesh
mesh = Mesh(input_dir+"gmsh.xml")
V    = FunctionSpace(mesh, 'CG', 1)

# Model parameters
T             = 2.0              # final time 
num_steps     = 100              # number of time steps
dt            = T/num_steps      # time step size
theta         = 50970.           # carrying capacity - normalize data by this
mu            = .42              # kPa, bulk shear modulus
nu            = .45              # poisson's ratio
lmbda         = 2*mu*nu/(1-2*nu) # lame parameter

# Load initial tumor condition data
initial_p = interp(input_dir+"ic.mat","ic")
initial_p.rename('initial','tumor at day 0')
target_p  = interp(input_dir+"tumor_t2.mat","tumor")  
target_p.rename('target','tumor at day 2')
f_notime.write(target_p)

annotate=False

# Parameters to be optimized
D0     = Constant(1.)     # mobility or diffusion coefficient
gammaD = Constant(1.)     # initial guess of gamma_D
beta   = Constant(1.)     # force coefficient for HE
k0     = Constant(2.)     # growth rate initial guess
k      = project(k0,V)

# Optimization module
[k, D0, gammaD, beta] = optimize() # optimize the k field, gammaD, and D0 using the adjoint method provided by adjoint_dolfin
f_log.write('Elapsed time is ' + str((time()-t1)/60) + ' minutes\n')

model_p = forward(initial_p,'opt',True, False) # run the forward model using the optimized k field

f_log.write('J_opt = '+str(objective(model_p, target_p, r_coeff1, r_coeff2))+'\n')
f_log.write('J_opt (without regularization) = '+str(objective(model_p, target_p, 0., 0.))+'\n')
f_log.write('D0 = '+str(D0.values()[0])+'\n')
f_log.write('gammaD = '+str(gammaD.values()[0])+'\n')
f_log.write('beta = '+str(beta.values()[0])+'\n')
f_log.close()
