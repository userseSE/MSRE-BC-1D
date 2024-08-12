from parameters import *
from ode_solver import ode_solver
# New Form Parameters
C1 = V_he2_s
C2 = -U2_hx / (M_he2_s * c_p_ss)
C3 = V_he2_ss
C4 = U2_hx / (M_he2_ss * c_p_sss)

# Discretize the spatial domain
A_HX2 = (-2*np.diag(np.ones(Nx)) + np.diag(np.ones(Nx - 1), 1) + np.diag(np.ones(Nx - 1), -1)) / dx
A_HX2[0, 0] = 1 / dx
A_HX2[-1, -1] = 1 / dx

def HX2(y_hx2, Ts_HX2_L, Tss_HX2_0, step):
    
    y_hx2[:Nx,-1]=Ts_HX2_L
    y_hx2[Nx:,0]=Tss_HX2_0
    
    # y is a vector with u and v concatenated
    def pde_to_ode_hx2(t,y):
        u=y[:Nx]
        v=y[Nx:]
        
        du_dt = C1 * (A_HX2 @ u) + C2 * (u - v)
        dv_dt = C3 * (A_HX2 @ v) + C4 * (u - v)
        
        # Apply time-varying boundary conditions
        du_dt[0] = u2_L - u[0]
        du_dt[-1] = u2_H - u[-1]
        dv_dt[0] = v2_L - v[0]
        dv_dt[-1] = v2_H - v[-1]
        
        dydt = np.concatenate([du_dt, dv_dt])
        return dydt
    
    # Initial condition vector
    if step==0:
        y0 = np.concatenate([u2_init, v2_init])
    else:
        y0 = y_hx2[:,-1]
        
    # solution_y_hx1 = solve_ivp(pde_to_ode_hx1, (step, step+1), y0, t_eval=t, method='RK45')     
    bc=[]
    solution_y_hx2 = ode_solver(y0, bc, pde_to_ode_hx2)
        
    # y_hx1 is a vector at time step+1
    y_hx2 = solution_y_hx2.y   
    # print("testHX1")
    # print(y_hx1.shape)
    return y_hx2