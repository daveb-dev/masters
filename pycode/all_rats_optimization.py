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

def forward(initial_p, name, record=False,  annt=False):
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
    global t
    I = Identity(2)  # Identity tensor
    def E(u):
        return 0.5*(nabla_grad(u) + nabla_grad(u).T)
    def vonmises(u):
        s = sigma(u) - (1./2)*tr(sigma(u))*I  # deviatoric stress
        von_Mises = sqrt(3./2*inner(s, s))
        return project(von_Mises, V, annotate=annt)

    #Set up linear elasticity problem
    U   = VectorFunctionSpace(mesh,'Lagrange',1)
    def boundary(x, on_boundary):
        return on_boundary
    bc  = DirichletBC(U, Constant((0.,0.)), boundary)
    p_n = interpolate(initial_p,V)
    v   = TestFunction(U)
    
    parameters['form_compiler']['quadrature_degree'] = 2
    parameters['form_compiler']['cpp_optimize'] = True
    parameters['krylov_solver']['nonzero_initial_guess'] = True
    ffc_options = {"quadrature_degree": 2, 'cpp_optimize': True}
    
    if lin_hyp == 0:
        def sigma(u):
            s = 2*mu*E(u)+lmbda*tr(E(u))*I
            return s
        u    = TrialFunction(U)
        a = inner(2*mu*E(u)+lmbda*tr(E(u))*I,E(v))*dx
        L = inner(2*beta*p_n*I*(mu+lmbda),E(v))*dx
        u    = Function(U, annotate=annt)
        def mech():
            solve(a == L, u, bc, 
                      form_compiler_parameters=ffc_options,
                      annotate=annt)
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
        
        u           = Function(U, annotate=annt)
        du          = TrialFunction(U)
        F_HE        = inner(sigma_form(u, p_n), E(v))*dx
        J_HE        = derivative(F_HE,u,du)
        problem_HE  = NonlinearVariationalProblem(F_HE, u, bc,
                                      J=J_HE,
                                      form_compiler_parameters=ffc_options)
        solver_HE   = NonlinearVariationalSolver(problem_HE)
        param_HE = solver_HE.parameters
        set_nonlinear_params(param_HE)
        def mech():
            solver_HE.solve(annotate=annt)
            return u

    # First iteration solving for displacement, 
    #  and using the von mises stress field for D
    disp = mech()
    vm   = vonmises(disp)
    D    = project(D0*exp(-gammaD*vm),V,annotate=annt)
    #k    = project(k0*exp(-gammaK*vm),V,annotate=annotate)
    
    if record: 
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
            
    # Set up reaction-diffusion problem
    dp   = TrialFunction(V)
    p    = Function(V,annotate=annt)
    q    = TestFunction(V)
    F_RD = (1/dt)*(p - p_n)*q*dx + D*dot(grad(q),grad(p))*dx - k*p*(1 - p)*q*dx  
    J_RD = derivative(F_RD,p,dp) 
    
    for n in range(num_steps):
        # Solve reaction diffusion
        t += dt
        problem_RD  = NonlinearVariationalProblem(F_RD, p,
                                                  J=J_RD,
                                                  form_compiler_parameters=ffc_options)
        solver_RD   = NonlinearVariationalSolver(problem_RD)
        param_RD = solver_RD.parameters
        set_nonlinear_params(param_RD)
        solver_RD.solve(annotate=annt)
        p_n.assign(p)
   
        # Solve for displacement and vonmises stress
        disp = mech()
        vm   = vonmises(disp)
        D    = project(D0*exp(-gammaD*vm),V,annotate=annt)
        #k    = project(k0*exp(-gammaK*vm),V,annotate=annotate)
        
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

def objective(p, target_p):
    return assemble(inner(p-target_p, p-target_p)*dx) \
       
def optimize(dbg=False):
    # Define the control
    m = [Control(D0), Control(gammaD), Control(k), Control(beta)]
    
    # Execute first time to annotate and record the tape
    p = forward(initial_p, 'annt', True, True)

    Obj = objective(p, target_p)

    # Prepare the reduced functional
    rf = ReducedFunctional(Obj,m,eval_cb_post=eval_cb)
    
    # upper and lower bound for the parameter field
    D_lb = 0.
    D_ub = 5.
    k_lb, k_ub = Function(V,annotate=False), Function(V,annotate=False)
    k_lb.vector()[:] = 0.
    k_ub.vector()[:] = 5.
    gD_lb = .01
    gD_ub = 1.
    beta_lb = 0.01
    beta_ub = 1.
    bnds = [[D_lb, gD_lb, k_lb, beta_lb],[D_ub, gD_ub, k_ub, beta_ub]]

    # Run the optimization
    m_opt = minimize(rf,method='L-BFGS-B', bounds=bnds, 
                     options={"disp":True,
                              "gtol":2.0e-5,
                              "ftol":2.0e-7,
                              "maxiter":100,
                              "maxls": 15})
    
    return m_opt

#########################################################################
# MAIN 
# call the function with: python <this file> 
########################################################################
    
rat_nums = ["01","02","05","06","09","12"]
for index, rat_num in enumerate(rat_nums):
    
    input_dir  = "../rat-data/rat"+rat_num+"/"
    output_dir = "./output/"
    f_log = open(osjoin(output_dir,'log.txt'),'w+')
    
    # Days data and time steps
    alldata  = sc_io_loadmat("../rat-data/finaldata.mat",)
    days     = alldata['rat'][0][index][3][0]
    days[:]  = [x-days[0] for x in days];
    steps    = []
    for i in range(1,len(days)):
        steps.append(days[i]-days[i-1])
        
    # Constant inputs for optimization
    D0       = 1.
    gammaD   = .5
    k0       = 1.
    beta     = .5
    theta    = 50970.       # carrying capacity - normalize cell data by this 
    mu       = .42          # kPa, bulk shear modulus
    nu       = .45
    lmbda    = 2*mu*nu/(1-2*nu)
    t1       = time()
    rtime    = 1 # How often to record results
    
    #fl = h5py.File("./output/rat05le_day2/notime.h5", "r")
    
    # Prepare a mesh
    mesh = Mesh(input_dir+"gmsh.xml")
    V    = FunctionSpace(mesh, 'CG', 1)
    
    # Load initial tumor condition data
    initial_p = interp(input_dir+"tumor_t0.mat","tumor")
    initial_p.rename('initial','tumor at day 0')
    
    for day in days[1:4]:
        
        for lin_hyp in [0,1]:
            
            # Model parameters
            t         = 0.           # initial time 
            T         = day           # final time 
            num_steps = T*10           # number of time steps
            dt        = T/float(num_steps)  # time step size

            # Load tumor condition data for day to optimize
            target_p  = interp(input_dir+"tumor_t"+str(day)+".mat","tumor")  
            target_p.rename('p_day'+str(day),'tumor at day '+str(day))

            # Prepare output file
            rat_id = rat_num+"-"+str(day)
            f_timeseries = XDMFFile(osjoin(output_dir,rat_id+"timeseries.xdmf"))
            f_timeseries.parameters["flush_output"] = True
            f_timeseries.parameters["functions_share_mesh"] = True
            f_nosteps    = XDMFFile(osjoin(output_dir,rat_id+"nosteps.xdmf"))
            f_nosteps.parameters["flush_output"] = True
            f_nosteps.parameters["functions_share_mesh"] = True
            f_notime     = XDMFFile(osjoin(output_dir,rat_id+"notime.xdmf"))
            f_notime.parameters["flush_output"] = True
            f_notime.parameters["functions_share_mesh"] = True

            ####################### D0 as a constant and k as a field ##########################
            # Initial guesses
            D0     = Constant(D0)
            gammaD = Constant(gammaD)     # initial guess of gamma_D
            k      = project(Constant(k0),V,annotate=False)     # growth rate initial guess
            beta   = Constant(beta)

            # Optimization 
            [D0, gammaD, k, beta] = optimize() # optimize these params using the adjoint method provided by adjoint_dolfin

            # Record time and optimized values
            f_log.write('Rat'+rat_num)
            f_log.write('Linear(0) or Hyper(1): '+str(lin_hyp))
            f_log.write('Day used for optimization: '+str(day))
            f_log.write('Elapsed time is ' + str((time()-t1)/60) + ' minutes\n')   
            f_log.write('gammaD = '+str(gammaD.values()[0])+'\n')
            f_log.write('D0 = '+str(D0.values()[0])+'\n')
            f_log.write('beta = '+str(beta.values()[0])+'\n')
            f_log.write('\n')

            k.rename('k0','diffusion field')
            f_notime.write(k,0.)

            # Compare optimized tumor growth to actual at several time points
            t   = 0.

            model_p = initial_p  # Initialize
            model_p.rename('opt_p','optimized tumor')
            f_nosteps.write(model_p, 0.)

            target_p = initial_p  # Initialize
            target_p.rename('true_p', 'optimized tumor')
            f_nosteps.write(target_p, 0.)

            for index2, step in steps:
                num_steps = step*10              # number of time steps
                dt        = step/float(num_steps)      # time step size

                # Run forward model using optimized values
                model_p = forward(model_p,'opt',True,False) 
                model_p.rename('opt_p','optimized tumor')
                f_nosteps.write(model_p,float(days[index2+1]))

                # Save actual tumor for comparison
                target_p = interp(input_dir+"tumor_t"+str(days[index2+1])+".mat","tumor")
                target_p.rename('true_p','actual tumor')
                f_nosteps.write(target_p,float(days[index2+1]))

                # Save J_opt
                f_log.write('J_opt day '+str(days[index2+1])+' = '
                            +str(objective(model_p, target_p))+'\n')

    f_log.close()

