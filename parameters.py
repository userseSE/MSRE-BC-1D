import numpy as np
from scipy.sparse import csc_matrix

# Neutronics
dt=0.2 # fixed time step
L = 172  # Length of the spatial domain, m
N = 200  # Number of spatial points
dz = L / (N - 1)    # Spatial step size
V = 4e5             
D = 0.390016          
sigma_a =0.0835           
nu_sigma_f = 3.33029e-2
# nu_sigma_f = 3.33029
# sigma_f= nu_sigma_f/2.41   
beta = [0.000228, 0.000788, 0.000664, 0.000736, 0.000136, 0.000088] # Delayed neutron fractions
Beta = sum(beta)
delta=Beta*nu_sigma_f
# beta2 = [0.000228, 0.000788, 0.000664, 0.000736, 0.000136, 0.000088]
lambda_i = [0.0126, 0.0337, 0.139, 0.325, 1.13, 2.5]    # Decay constants
# initial condition
phi_0=5226.54 * 1e6 * np.ones(N); #5226.54
c0 = (delta / sum(lambda_i)) * phi_0

# Thermal-Hydraulics
c_p_s = 1983  # Specific heat of primary salt, J/kgK
Vc = 0.2    # Salt velocity in the core, m/s
Ms = 1448   # Fuel salt mass in the core, kg
Mg = 3687   # Mass of graphite in the core, kg
gamma = 0.93    # Fraction of power released in the salt
# U = 36000   # Overall heat transfer coefficient between salt and graphite, W/K
U = 18000
c_p_g = 1757    # Specific heat of graphite, J/kg K
# Temp input for test:
# Amplitude = 4.5e6;
# q_prime=Amplitude * sin(pi * linspace(0, L, N) / L)';
# q_prime=Amplitude;
# Boundary conditions
# bc_s0 = 910
# bc_sL = 958.15
# bc_g0 = 920
# bc_gL = 968.71
bc_s0 = 500
bc_sL = 500
bc_g0 = 800
bc_gL = 800
# bc_g0 = 931.15
# bc_gL = 931.15
# Initial conditions
initialS = (bc_s0 + (bc_sL - bc_s0) * (0.5 + 0.5 * np.sin(np.pi * (np.linspace(0, L, N) ) / L)) * 0.8).T
initialG = (bc_g0 + (bc_gL - bc_g0) * (0.5 + 0.5 * np.sin(np.pi * (np.linspace(0, L, N) ) / L)) * 1.05).T
# referenceS=930
# referenceG=931.15
# Heat Exchanger 1
err = 1e-3
Nx = N  # Number of spatial points
L_HX = 2    # length of the spatial domain
dx = L / (N - 1)    # Spatial step size
# V_s0 = Vc   # Salt velocity in the core m/s
V_he_s = 75.7 * 1e-3 / 23.6 #0.003207627118644068 m/s
V_he_ss = 53.6 * 1e-3 / 23.6    # 0.002271186440677966 m/s
U_hx = 82800    # Heat transfer coefficient between primary and secondary salt, W/K
M_he_s = 342    # Mass of salt in the heat exchanger (fuel side), kg
M_he_ss = 117   # Mass of salt in the heat exchanger (coolant side), kg
c_p_ss = 2416   # Specific heat of secondary salt, J/kg/K
# Initial conditions
u_L = 900.15
u_H = 938
v_L = 824.85
v_H = 866.45
# initial conditions using a sine function
u_init = u_L + (u_H - u_L) * (0.5 + 0.5 * np.sin(np.pi * (np.linspace(0, L_HX, Nx) / L_HX))).T
v_init = v_L + (v_H - v_L) * (0.5 + 0.5 * np.sin(np.pi * (np.linspace(0, L_HX, Nx) / L_HX))* 1.05).T

# Heat Exchanger 2, TODO: change the parameters
L_HX2 = 2    # length of the spatial domain
V_he2_s = 53.6 * 1e-3 / 23.6 #0.003207627118644068 m/s
V_he2_ss = 33.6 * 1e-3 / 23.6    # 0.002271186440677966 m/s
U2_hx = 82800    # Heat transfer coefficient between primary and secondary salt, W/K
M_he2_s = 117    # Mass of salt in the heat exchanger (fuel side), kg
M_he2_ss = 100   # Mass of salt in the heat exchanger (coolant side), kg
c_p_sss = 2416   # Specific heat of secondary salt, J/kg/K
# Initial conditions
u2_L = 824
u2_H = 866
v2_L = 744
v2_H = 786
# initial conditions using a sine function
u2_init = u2_L + (u2_H - u2_L) * (0.5 + 0.7 * np.sin(np.pi * (np.linspace(0, L_HX2, Nx) / L_HX2))).T
v2_init = v2_L + (v2_H - v2_L) * (0.5 + 0.7 * np.sin(np.pi * (np.linspace(0, L_HX2, Nx) / L_HX2))* 1.05).T

# Reactivity
rho_init=50 * np.ones(N)
alpha_f    = -5.904E-5  # U233 (drho/K) fuel salt temperature-reactivity feedback coefficient ORNL-TM-1647 p.3 % -5.904E-05; % ORNL-TM-0728 p. 101 %
alpha_g    = -6.624E-5  # U233  (drho/K) graphite temperature-reactivity feedback coefficient ORNL-TM-1647 p.3 % -6.624E-05; % ORNL-TM-0728 p.101
tau_l  = 16.73  # ORNL-TM-0728 %16.44; % (s)
tau_c  = 8.46   # ORNL-TM-0728 %8.460; % (s)

# Transport Delays
# Pure time delays between components
tau_hx_c = 9 # (sec) delay from hx to core TDAMSRE p.6
tau_c_hx = 4 # (sec) subtracted 1 sec for external loop power generation node resident time; delay from core to fuel hx TDAMSRE p.6
tau_hx_r = 5 # (sec) fertile hx to core TDAMSRE p.6
tau_r_hx = 8 # (sec) core to fertile hx TDAMSRE p.6
tau_r_pp = 6 # TODO:delay from HX2 to the power plant, to be fixed
tau_pp_r = 6 # TODO:delay from the power plant to HX2, to be fixed
# Initial conditions - DONE
Ts_in=bc_s0
Ts_out=bc_sL
Tss_in=v_L
Tss_out=v_H
Tsss_in=v2_L
Tsss_out=v2_H