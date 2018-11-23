from dolfin import *
from numpy import fliplr, linspace, inf
from os.path import join as osjoin
from scipy.io import loadmat as sc_io_loadmat
from scipy.interpolate import RegularGridInterpolator
import sys

set_log_level(ERROR) 

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
    mat = sc_io_loadmat(file_loc)[mat_name]
    mat = fliplr(mat.T)/theta  # Needs to be adjusted to fit the mesh correctly; also scaled
    x,y = mat.shape[0], mat.shape[1]
    mat_interp = InterpolatedParameter(linspace(1,x,x),linspace(1,y,y),mat,degree=1)
    return interpolate(mat_interp,V)

def forward(initial_p, name):    
    """ 
        Here, we define the forward problem. 
    """
    '''
        - E(u) returns the Green-Lagrange strain tensor
        - sigma(...) returns the actual stress tensor
        - sigma_form(...) returns the stress tensor based on the cells (phi), elasticity coefficients, and a coefficient beta
        - vonmises(...) calculates the von Mises stress based on the actual stress tensor
    '''
    
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
        return project(von_Mises, V)
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
    u           = Function(U) 
    v           = TestFunction(U)
    p_n         = interpolate(initial_p,V)
    F_HE        = inner(sigma_form(u, p_n), E(v))*dx
    J_HE        = derivative(F_HE,u,du)
    
    ffc_options = {"quadrature_degree": 2, "cpp_optimize": True}
    parameters['form_compiler']['quadrature_degree'] = 2
    parameters['krylov_solver']['nonzero_initial_guess'] = True
    
    problem_HE  = NonlinearVariationalProblem(F_HE, u, bc, J=J_HE,form_compiler_parameters=ffc_options)
    solver_HE   = NonlinearVariationalSolver(problem_HE)
    prm1 = solver_HE.parameters
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
    def he():
        solver_HE.solve()
        return u
    
    # First iteration solving for displacement, and using the von mises stress field for D
    disp = he()
    vm   = vonmises(disp)
    D    = project(D0*exp(-gammaD*vm),V)
    
    # Set up reaction-diffusion problem
    dp          = TrialFunction(V)
    p           = Function(V)
    q           = TestFunction(V)
    F_RD        = (1/dt)*(p - p_n)*q*dx + D*dot(grad(q),grad(p))*dx - k*p*(1 - p)*q*dx  
    J_RD = derivative(F_RD,p) 
    
    # Prepare the solution
    t = 0.
    for n in range(num_steps):        
        # Update current time and Compute solution
        t += dt
        problem_RD  = NonlinearVariationalProblem(F_RD, p, J=J_RD,form_compiler_parameters=ffc_options)
        solver_RD   = NonlinearVariationalSolver(problem_RD)
        prm = solver_RD.parameters
        prm['newton_solver']['absolute_tolerance'] = 1E-7
        prm['newton_solver']['relative_tolerance'] = 1E-6
        prm['newton_solver']['maximum_iterations'] = 51
        prm['newton_solver']['relaxation_parameter'] = 1.0
        prm['newton_solver']['linear_solver'] = 'gmres'
        prm['newton_solver']['preconditioner'] = 'ilu'
        prm['newton_solver']['krylov_solver']['absolute_tolerance'] = 1E-8
        prm['newton_solver']['krylov_solver']['relative_tolerance'] = 1E-6
        prm['newton_solver']['krylov_solver']['maximum_iterations'] = 1000
        prm['newton_solver']['krylov_solver']['nonzero_initial_guess'] = True
        solver_RD.solve()
        
        # Update previous solution
        p_n.assign(p)
        disp = he()
        vm   = vonmises(disp)
        D    = project(D0*exp(-gammaD*vm),V)
        
        if (n%2 == 0):
            u.rename('u_'+name,'displacement')
            p_n.rename('phi_T_'+name,'tumor fraction')
            vm.rename('vm_'+name,"Von Mises")
            D.rename('D_'+name,"diffusion coefficient")
            k.rename('k_'+name,'k field') 
            f_notime.write(p_n,t)
            f_notime.write(u,t)
            f_notime.write(k,t)
            f_notime.write(vm,t)
            f_notime.write(D,t)

#########################################################################
# MAIN 
########################################################################
if __name__ == "__main__":
    # call the function with:
    # python <this file> case r_coeff1 r_coeff2
    if(len(sys.argv) != 6):
        print("wrong number of inputs, should be ")
        print("Syntax: python <this file's name> D0 gammaD beta k0 case ")
        quit()
    D0       = float(sys.argv[1])
    gammaD   = float(sys.argv[2])
    beta     = float(sys.argv[3])
    k0       = float(sys.argv[4])
    case     = sys.argv[5]
    
    t1         = time()
    input_dir  = "../rat-data/rat05/"
    output_dir = './output/rat05hesense'

    # Prepare a mesh
    mesh = Mesh(input_dir+"gmsh.xml")
    V    = FunctionSpace(mesh, 'CG', 1)

    # Model parameters
    T             = 2.0              # final time 
    num_steps     = 20              # number of time steps
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

    # Parameters to be optimized
    D0     = Constant(D0)     # mobility or diffusion coefficient
    gammaD = Constant(gammaD)     # initial guess of gamma_D
    beta   = Constant(beta)     # force coefficient for HE
    k0     = Constant(k0)     # growth rate initial guess
    k      = project(k0,V)

    # Prepare output file
    f_notime     = XDMFFile(osjoin(output_dir,'notime.xdmf'))
    f_notime.parameters["flush_output"] = True
    f_notime.parameters["functions_share_mesh"] = True
    
    # run the forward model
    forward(initial_p, str(case)) 

