import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import diags, csc_matrix
from scipy.sparse.linalg import spsolve

from parameters import *
from ode_solver import ode_solver

def HX1(y_hx1, Ts_HX1_L, Tss_HX1_0, params, step):
    
    V_he_s = params['V_he_s']
    U_hx = params['U_hx']
    M_he_s = params['M_he_s']
    c_p_s = params['c_p_s']
    V_he_ss = params['V_he_ss']
    U_hx = params['U_hx']
    M_he_ss = params['M_he_ss']
    c_p_ss = params['c_p_ss']
    Nx = params['Nx']
    dx = params['dx']
    u_L = params['u_L']
    u_H = params['u_H']
    v_L = params['v_L']
    v_H = params['v_H']
    err = params['err']
    u_init = params['u_init']
    v_init = params['v_init']
    h_s = params['h_s']
    k_s = params['k_s']
    h_ss = params['h_ss']
    k_ss = params['k_ss']
    T_ambient = params['T_ambient']
    
    # New Form Parameters
    C1 = V_he_s
    C2 = -U_hx / (M_he_s * c_p_s)
    C3 = V_he_ss
    C4 = U_hx / (M_he_ss * c_p_ss)

    A_HX_sparse = params.get('A_HX_sparse')
    if A_HX_sparse is None:
        A_HX = np.diag(-2 * np.ones(Nx)) + np.diag(np.ones(Nx - 1), 1) + np.diag(np.ones(Nx - 1), -1)
        A_HX[0, 0] = A_HX[-1, -1] = -1
        A_HX[0, 1] = A_HX[-1, -2] = 0
        A_HX_sparse = csc_matrix(A_HX) / dx**2
    
    y_hx1[:Nx,-1]=Ts_HX1_L
    y_hx1[Nx:,0]=Tss_HX1_0
    
    # y is a vector with u and v concatenated
    def pde_to_ode_hx1(t,y):
        u=y[:Nx]
        v=y[Nx:]
        
        du_dt = C1 * (A_HX_sparse @ u) + C2 * (u - v) + err
        dv_dt = C3 * (A_HX_sparse @ v) + C4 * (u - v) + err
        
        # Apply time-varying boundary conditions
        du_dt[0] = u_L - u[0] / dx
        du_dt[-1] = u_H - u[-1] / dx
        dv_dt[0] = v_L - v[0] / dx
        dv_dt[-1] = v_H - v[-1] / dx
        # du_dt[0] = (h_s / k_s) * (T_ambient - u[0]) - u[0] / dx  # Robin condition at the left boundary for `u`
        # du_dt[-1] = (h_s / k_s) * (T_ambient - u[-1]) - u[-1] / dx  # Robin condition at the right boundary for `u`

        # dv_dt[0] = (h_ss / k_ss) * (T_ambient - v[0]) - v[0] / dx  # Robin condition at the left boundary for `v`
        # dv_dt[-1] = (h_ss / k_ss) * (T_ambient - v[-1]) - v[-1] / dx  # Robin condition at the right boundary for `v`

        dydt = np.concatenate([du_dt, dv_dt])
        return dydt
    
    # Initial condition vector
    if step==0:
        y0 = np.concatenate([u_init, v_init])
    else:
        y0 = y_hx1[:,-1]
        
    # solution_y_hx1 = solve_ivp(pde_to_ode_hx1, (step, step+1), y0, t_eval=t, method='RK45')     
    bc=[]
    y_hx1 = ode_solver(y0, bc, pde_to_ode_hx1, params)
    
    # y_hx1 is a vector at time step+1
    # print("testHX1")
    # print(y_hx1.shape)
    return y_hx1
