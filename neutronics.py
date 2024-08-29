from scipy.sparse.linalg import spsolve
# from scipy.integrate import solve_ivp
import numpy as np
from numpy.linalg import inv

from parameters import *
from ode_solver import ode_solver

# discretize the spatial domain
# Define the finite difference matrix for the second derivative using Crank-Nicolson method
main_diag = -2 * np.ones(N)
off_diag = np.ones(N - 1)
D2 = (np.diag(main_diag) + np.diag(off_diag, 1) + np.diag(off_diag, -1)) / dz**2
# Boundary conditions: Dirichlet boundary conditions for zero flux at boundaries
D2[0, :] = 0
D2[-1,:] = 0
D2[0, 0] = -1
D2[-1, -1] = -1
# Convert to sparse matrix format
D2_sparse = csc_matrix(D2)
I = np.eye(N)
A = csc_matrix(I - 0.5 * dt * V * D * D2_sparse)
B = csc_matrix(I + 0.5 * dt * V * D * D2_sparse)

def neutronics(y_n, rho, step):
    Keff=1/(np.ones(N)-rho)
    def pde_to_ode_neutronics(t, y):
        
        phi = y[:N]
 
        # Apply Crank-Nicolson method
        lambda_ci = np.zeros(N)
        for i in range(6):
            lambda_ci += lambda_i[i]*y[(i+1)*N:(i+2)*N]
            
        rhs_phi = B @ phi + dt * V * ((-sigma_a + ((1 - Beta) * nu_sigma_f / Keff)) * phi + lambda_ci)
        # Check for invalid values
        # if np.any(np.isnan(rhs_φ)) or np.any(np.isinf(rhs_φ)):
        #     print(f"Invalid values encountered at t={t}")
        #     print(f"φ: {φ}")
        #     print(f"c: {c}")
        #     print(f"rhs_φ: {rhs_φ}")
        #     raise ValueError("Invalid values encountered in rhs_φ")

        # print(rhs_φ.shape)

        phi_new = spsolve(A, rhs_phi)
        dphi_dt = (phi_new - phi) / dt
    
        dci_dt = np.zeros((6, N))
        for i in range(6):
            dci_dt[i]=beta[i]*(nu_sigma_f/Keff)*phi-lambda_i[i]*y[(i+1)*N:(i+2)*N]
        return np.concatenate([dphi_dt, dci_dt[0], dci_dt[1], dci_dt[2], dci_dt[3], dci_dt[4], dci_dt[5]])

    # Initial condition vector
    if step == 0:
        y0 = np.concatenate([phi_0, c0, c0, c0, c0, c0, c0])
    else:
        y0 = y_n
    # Solve the system of ODEs
    # TODO: update the integration method
    # print(y0.shape)
    
    bc=[]
    solution_y_n = ode_solver(y0, bc, pde_to_ode_neutronics)
    # print("testNeutronics")
    # Extract the solution
    phi = solution_y_n.y[:N,-1].T
    q_prime = phi
    # print(q_prime.shape)
    y_n = solution_y_n.y
    # print(y_n.shape)
    
    return y_n, q_prime