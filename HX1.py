import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import diags, csc_matrix
from scipy.sparse.linalg import spsolve

from parameters import *
from ode_solver import ode_solver

# New Form Parameters
C1 = V_he_s
C2 = -U_hx / (M_he_s * c_p_s)
C3 = V_he_ss
C4 = U_hx / (M_he_ss * c_p_ss)

# Discretize the spatial domain
A_HX = (-2*np.diag(np.ones(Nx)) + np.diag(np.ones(Nx - 1), 1) + np.diag(np.ones(Nx - 1), -1)) / dx
A_HX[0, 0] = 1 / dx
A_HX[-1, -1] = 1 / dx



def HX1(y_hx1, Ts_HX1_L, Tss_HX1_0, step):
    
    y_hx1[:Nx,-1]=Ts_HX1_L
    y_hx1[Nx:,0]=Tss_HX1_0
    
    # y is a vector with u and v concatenated
    def pde_to_ode_hx1(t,y):
        u=y[:Nx]
        v=y[Nx:]
        
        du_dt = C1 * (A_HX @ u) + C2 * (u - v)
        dv_dt = C3 * (A_HX @ v) + C4 * (u - v)
        
        # Apply time-varying boundary conditions
        du_dt[0] = u_L - u[0]
        du_dt[-1] = u_H - u[-1]
        dv_dt[0] = v_L - v[0]
        dv_dt[-1] = v_H - v[-1]
        
        dydt = np.concatenate([du_dt, dv_dt])
        return dydt
    
    # Initial condition vector
    if step==0:
        y0 = np.concatenate([u_init, v_init])
    else:
        y0 = y_hx1[:,-1]
        
    # solution_y_hx1 = solve_ivp(pde_to_ode_hx1, (step, step+1), y0, t_eval=t, method='RK45')     
    bc=[]
    solution_y_hx1 = ode_solver(y0, bc, pde_to_ode_hx1)
        
    # y_hx1 is a vector at time step+1
    y_hx1 = solution_y_hx1.y   
    # print("testHX1")
    # print(y_hx1.shape)
    return y_hx1