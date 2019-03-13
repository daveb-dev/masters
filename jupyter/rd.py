## Reaction Diffusion Formulation in 2D
from fenics import *

# Define and solve variational problem
def rd(dt, V, phi_n, D, alpha):
    phi = Function(V)
    v = TestFunction(V)
    k = Constant(dt)

    F = (phi - phi_n)*v*dx \
        + k*D*dot(grad(v),grad(phi))*dx \
        - k*alpha*phi*(1 - phi)*v*dx 
    
    solve(F == 0, phi)
    return phi
   