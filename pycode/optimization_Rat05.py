from dolfin import *
from dolfin_adjoint import *
from numpy import fliplr, linspace, inf
from os.path import join as osjoin
from scipy.io import loadmat as sc_io_loadmat
from scipy.interpolate import RegularGridInterpolator

set_log_level(ERROR) 

class InterpolatedParameter(Expression):
    '''
        Class to get tumor cell distributions by interpolating 
        based off matrices of tumor cell data
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
        Function to accept matlab .mat file with tumor data 
        and interpolate values onto mesh
    """
    mat = sc_io_loadmat(file_loc)[mat_name]
    mat = fliplr(mat.T)/theta  # Needs to be adjusted to fit the mesh correctly
    x,y = mat.shape[0], mat.shape[1]
    mat_interp = InterpolatedParameter(linspace(1,x,x),linspace(1,y,y),mat,degree=1)
    return interpolate(mat_interp,V)

def set_nonlinear_params(param):
    param['newton_solver']['absolute_tolerance'] = 1E-7
    param['newton_solver']['relative_tolerance'] = 1E-6
    param['newton_solver']['maximum_iterations'] = 51
    param['newton_solver']['relaxation_parameter'] = 1.0
    param['newton_solver']['linear_solver'] = 'gmres'
    param['newton_solver']['preconditioner'] = 'ilu'
    param['newton_solver']['krylov_solver']['absolute_tolerance'] = 1E-8
    param['newton_solver']['krylov_solver']['relative_tolerance'] = 1E-6
    param['newton_solver']['krylov_solver']['maximum_iterations'] = 1000
    param['newton_solver']['krylov_solver']['nonzero_initial_guess'] = True
    

def forward(initial_p, name, record=False,  annotate=False):
    """ 
        Here, we define the forward problem with mechanical functions
        
        -E(u) returns the Green-Lagrange strain tensor
        -sigma(...) returns the actual stress tensor
        -sigma_form(...) returns the stress tensor based on 
          the cells (phi), elasticity coefficients, and a 
          coefficient beta 
        -vonmises(...) calculates the von Mises stress based 
          on the actual stress tensor
    """
    I = Identity(2)  # Identity tensor
    def E(u):
        return 0.5*(nabla_grad(u) + nabla_grad(u).T)
    def vonmises(u):
        s = sigma(u) - (1./2)*tr(sigma(u))*I  # deviatoric stress
        von_Mises = sqrt(3./2*inner(s, s))
        return project(von_Mises, V, annotate=annotate)

    #Set up linear elasticity problem
    U   = VectorFunctionSpace(mesh,'Lagrange',1)
    def boundary(x, on_boundary):
        return on_boundary
    bc  = DirichletBC(U, Constant((0.,0.)), boundary)
    p_n = interpolate(initial_p,V)
    v   = TestFunction(U)
    
    parameters["form_compiler"]["quadrature_degree"] = 2
    parameters["form_compiler"]["cpp_optimize"] = True
    parameters['krylov_solver']['nonzero_initial_guess'] = True
    ffc_options = {"quadrature_degree": 2}
    
    if lin_hyp == 0:
        def sigma(u):
            s = 2*mu*E(u)+lmbda*tr(E(u))*I
            return s
        u    = TrialFunction(U)
        a = inner(2*mu*E(u)+lmbda*tr(E(u))*I,E(v))*dx
        L = inner(2*beta*p_n*I*(mu+lmbda),E(v))*dx
        u    = Function(U, annotate=annotate)
        def mech():
            solve(a == L, u, bc, 
                      form_compiler_parameters=ffc_options,
                      annotate=annotate)
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
        
        u           = Function(U, annotate=annotate)
        du          = TrialFunction(U)
        F_HE        = inner(sigma_form_he(u, p_n), E(v))*dx
        J_HE        = derivative(F_HE,u,du)
        problem_HE  = NonlinearVariationalProblem(F_HE, u, bc,
                                      J=J_HE,
                                      form_compiler_parameters=ffc_options)
        solver_HE   = NonlinearVariationalSolver(problem_HE)
        param_HE = solver_HE.parameters
        set_nonlinear_params(param_HE)
        def mech():
            solver_HE.solve(annotate=annotate)
            return u

    # First iteration solving for displacement, 
    #  and using the von mises stress field for D
    disp = mech()
    vm   = vonmises(disp)
    D    = project(D0*exp(-gammaD*vm),V,annotate=annotate)

    # Set up reaction-diffusion problem
    dp   = TrialFunction(V)
    p    = Function(V,annotate=annotate)
    q    = TestFunction(V)
    F_RD = (1/dt)*(p - p_n)*q*dx + D*dot(grad(q),grad(p))*dx - k*p*(1 - p)*q*dx  
    J_RD = derivative(F_RD,p,dp) 
    
    # Prepare the solution
    t = 0.                
    for n in range(num_steps+1):
        if record and (n%rtime == 0): 
            # Rename parameters for saving
            u.rename('u_'+name,'displacement')
            p_n.rename('phi_T_'+name,'tumor fraction')
            vm.rename('vm_'+name,'Von Mises')
            D.rename('D_'+name,'diffusion coefficient')
            k.rename('k_'+name,'k field')
            f_timeseries.write(u,t)          
            f_timeseries.write(p_n,t)
            f_timeseries.write(vm,t)
            f_timeseries.write(D,t)
            f_timeseries.write(k,t)
        
        # Solve reaction diffusion
        t += dt
        problem_RD  = NonlinearVariationalProblem(F_RD, p,
                                                  J=J_RD,
                                                  form_compiler_parameters=ffc_options)
        solver_RD   = NonlinearVariationalSolver(problem_RD)
        param_RD = solver_RD.parameters
        set_nonlinear_params(param_RD)
        solver_RD.solve(annotate=annotate)
        p_n.assign(p)
        
        # Solve for displacement and vonmises stress
        disp = mech()
        vm   = vonmises(disp)
        D    = project(D0*exp(-gammaD*vm),V,annotate=annotate)
        
    if record and (n%rtime == 0): 
        # Rename parameters for saving
        u.rename('u_'+name,'displacement')
        p_n.rename('phi_T_'+name,'tumor fraction')
        vm.rename('vm_'+name,'Von Mises')
        D.rename('D_'+name,'diffusion coefficient')
        k.rename('k_'+name,'k field')
        f_timeseries.write(u,t)          
        f_timeseries.write(p_n,t)
        f_timeseries.write(vm,t)
        f_timeseries.write(D,t)
        f_timeseries.write(k,t)
    return p

# Callback function for the optimizer
# Writes intermediate results to a logfile
def eval_cb(j, m):
    """ The callback function keeping a log """
    f_log.write("objective = %15.10e \n" % j)

def objective(p, target_p, r_coeff1, r_coeff2, r_coeff3):
    return assemble(inner(p-target_p, p-target_p)*dx,annotate=False) \
        + r_coeff1*assemble(k*k*dx,annotate=False) \
        + r_coeff2*assemble(dot(grad(k),grad(k))*dx,annotate=False) \
        + r_coeff3*assemble(sqrt(dot(grad(p),grad(p))+eps)*dx)
        # assemble(inner(p-target_p, p-target_p)*dx) + r_coeff1*assemble(k*k*dx) + r_coeff2*assemble(dot(grad(k),grad(k))*dx)

def optimize(dbg=False):
    # Define the control
    m = [Control(k), Control(D0)]
    
    # Execute first time to annotate and record the tape
    p = forward(initial_p, 'annt', True, True)

    Obj = objective(p, target_p, r_coeff1, r_coeff2, r_coeff3)

    # Prepare the reduced functional
    rf = ReducedFunctional(Obj,m,eval_cb_post=eval_cb)
    
    # upper and lower bound for the parameter field
    k_lb, k_ub = Function(V,annotate=False), Function(V,annotate=False)
    k_lb.vector()[:] = 0.
    k_ub.vector()[:] = 10.
    D_lb = 0.
    D_ub = 5.
    gD_lb = 0
    gD_ub = 5.
    beta_lb = 0.
    beta_ub = 5.
    #bnds = [[k_lb,D_lb, gD_lb, beta_lb],[k_ub,D_ub, gD_ub, beta_ub]]
    bnds = [[k_lb,D_lb],[k_ub,D_ub]]

    # Run the optimization
    m_opt = minimize(rf,method='L-BFGS-B', bounds=bnds, 
                     options={"disp":True,"gtol":2.0e-4,"ftol":2.0e-6})
    
    return m_opt

#########################################################################
# MAIN 
########################################################################
if __name__ == "__main__":
    # call the function with:
    # python <this file> case r_coeff1 r_coeff2
    if(len(sys.argv) != 6):
        print("wrong number of inputs, should be:\n ")
        print("Syntax: python <this file's name> [0=LE/1=HE] D0 gammaD beta k0 ")
        quit()
    lin_hyp  = int(sys.argv[1])
    D0       = float(sys.argv[2])
    gammaD   = float(sys.argv[3])
    beta     = float(sys.argv[4])
    k0       = float(sys.argv[5])
    
    t1         = time()
    r_coeff1   = 0.01
    r_coeff2   = 0.01
    r_coeff3 = 1.
    eps = .01
    
    input_dir  = "../rat-data/rat05/"
    if lin_hyp == 0:
        output_dir = './output/rat05le'
    else:
        output_dir = './output/rat05he'

    # Prepare output file
    f_timeseries = XDMFFile(osjoin(output_dir,'timeseries.xdmf'))
    f_timeseries.parameters["flush_output"] = True
    f_timeseries.parameters["functions_share_mesh"] = True
    f_nosteps    = XDMFFile(osjoin(output_dir,'nosteps.xdmf'))
    f_nosteps.parameters["flush_output"] = True
    f_nosteps.parameters["functions_share_mesh"] = True
    f_log = open(osjoin(output_dir,'log.txt'),'w+')
    rtime = 2 # How often to record results

    # Prepare a mesh
    mesh = Mesh(input_dir+"gmsh.xml")
    V    = FunctionSpace(mesh, 'CG', 1)

    # Model parameters
    T             = 2.               # final time 
    num_steps     = 20              # number of time steps
    dt            = T/num_steps      # time step size
    theta         = 50970.           # carrying capacity - normalize cell data by this 
    mu            = .42              # kPa, bulk shear modulus
    nu            = .45
    lmbda         = 2*mu*nu/(1-2*nu)

    # Load initial tumor condition data
    initial_p = interp(input_dir+"ic.mat","ic")
    initial_p.rename('initial','tumor at day 0')
    target_p  = interp(input_dir+"tumor_t2.mat","tumor")  
    target_p.rename('p_day2','tumor at day 2')

    annotate=False

    # Initial guesses
    D0     = Constant(D0)     # mobility or diffusion coefficient
    gammaD = Constant(gammaD)     # initial guess of gamma_D
    beta   = Constant(beta)     # force coefficient for HE
    k0     = Constant(k0)     # growth rate initial guess
    k      = project(k0,V, annotate=False)

    # Optimization module
    [k, D0] = optimize() # optimize the k field, gammaD, and D0 using the adjoint method provided by adjoint_dolfin
    
    # Record time and optimized values
    f_log.write('Elapsed time is ' + str((time()-t1)/60) + ' minutes\n')
    f_log.write('D0 = '+str(D0.values()[0])+'\n')
    f_log.write('gammaD = '+str(gammaD.values()[0])+'\n')
    f_log.write('beta = '+str(beta.values()[0])+'\n')
    
    # Compare optimized tumor growth to actual at several time points
    model_p = initial_p  # Initialize
    day = 0
    for T in [2,2,1,1,3]:
        day      += T
        num_steps = T*10              # number of time steps
        dt        = T/float(num_steps)      # time step size
        
        # Run forward model using optimized values
        model_p = forward(model_p,'opt',True,False) 
        model_p.rename('opt_p','optimized tumor')
        f_nosteps.write(model_p,day)
        
        # Save actual tumor for comparison
        target_p = interp(input_dir+"tumor_t"+str(day)+".mat","tumor")
        target_p.rename('true_p','actual tumor')
        f_nosteps.write(target_p,day)
        
        # Save J_opt
        f_log.write('J_opt day '+str(day)+' = '+str(objective(model_p, target_p,
                                                              r_coeff1, r_coeff2,
                                                              r_coeff3))+'\n')
        f_log.write('J_opt day '+str(day)+' (no regularization) = ' 
                    +str(objective(model_p, target_p, 0., 0.,0.))+'\n')
        
    f_log.close()
    
    