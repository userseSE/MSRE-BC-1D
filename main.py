import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import diags, csc_matrix
from scipy.sparse.linalg import spsolve
import queue
import matplotlib.pyplot as plt

from parameters import *
from reactivity import reactivity
from neutronics import neutronics
from thermal_hydraulics import thermal_hydraulics
from HX1 import HX1
from HX2 import HX2
from transport_delay import transport_delay

time_span = 1000

# reactivity insertion rod
rho_insertion = 50  # pcm

# initialization
rho=0
# initial variables
y_n=np.zeros((7*N, 1))
q_prime=np.zeros((N, 1))
y_th=np.zeros((2*N, 1))
y_hx1=np.zeros((2*Nx, 1))
Tss_HX2_0=0
Ts_HX1_0=0
Tss_HX1_0=0
buffer_hx_c=[]
buffer_c_hx=[]
buffer_r_hx=[]
buffer_hx_r=[]

# initial plotting matrices
rho_matrix = np.zeros((time_span, 1))
phi_middle_matrix = np.zeros((time_span, 1))
temperature_fuel_middle_matrix = np.zeros((time_span, 1))

for step in range (time_span):
    print (f"Time step: {step}")
    # int(time_span/dt)
     
    y_n, q_prime = neutronics(y_n[:,-1], rho, step)
    # print("test1") 
    phi = y_n[:N ,-1].T
    ci = y_n[N:,-1].T
    phi_middle_matrix[step] = phi[int(N/2)]
    # print(phi.shape)
    # print(ci.shape)
    # print(y_n.shape)
    
    # tau_hx_c
    Ts_core_0=transport_delay(Ts_HX1_0, tau_hx_c, Ts_in, buffer_hx_c, step)
    # print("size of buffer_hx_c: " + str(len(buffer_hx_c)))
    y_th = thermal_hydraulics(y_th, q_prime, Ts_core_0, step)
    temperature_fuel = y_th[:N, -1].T
    temperature_graphite = y_th[N:, -1].T
    Ts_core_L = y_th[-1, -1]
    temperature_fuel_middle_matrix[step] = temperature_fuel[int(N/2)]
    
    # tau_c_hx
    Ts_HX1_L=transport_delay(Ts_core_L,tau_c_hx, Ts_out, buffer_c_hx, step)
    # print("size of buffer_c_hx: " + str(len(buffer_c_hx)))
    # tau_r_hx
    Tss_HX1_0=transport_delay(Tss_HX2_0, tau_r_hx, Tss_in, buffer_r_hx, step)
    # print("size of buffer_r_hx: " + str(len(buffer_r_hx)))
    y_hx1 = HX1(y_hx1, Ts_HX1_L, Tss_HX1_0, step)
    Ts_HX1= y_hx1[:Nx,-1]
    Tss_HX1= y_hx1[Nx:,-1]
    Ts_HX1_0 = Ts_HX1[0]
    Tss_HX1_L = Tss_HX1[-1]
    
    # tau_hx_r
    Tss_HX2_L=transport_delay(Tss_HX1_L, tau_hx_r, Tss_out, buffer_hx_r, step)
    # print("size of buffer_hx_r: " + str(len(buffer_hx_r)))
    Tss_HX2_0=HX2(Tss_HX2_L, step)
    
    rho=reactivity(temperature_fuel, temperature_graphite, step, time_span, rho_insertion)
    rho_matrix[step]=rho

# plot phi to axis N
z = np.linspace(0, L, N)
fig, ax = plt.subplots(2, 2, figsize=(14, 6))
im1 = ax[0, 0].plot(z, phi)
ax[0, 0].set_title('Neutron Flux')

for i in range(6):
    ax[0, 1].plot(z, ci[i*N:(i+1)*N])
ax[0, 1].set_title('Delayed Neutron Precursors')

ax[1, 0].plot(rho_matrix)
ax[1, 0].set_title('Reactivity')

ax[1, 1].plot(z, temperature_fuel)
ax[1, 1].plot(z, temperature_graphite)
ax[1, 1].set_title('Temperature in the core')
plt.show()

fig, ax = plt.subplots(1, 3, figsize=(14, 6))
ax[0].plot(z, Ts_HX1)
ax[0].plot(z, Tss_HX1)
ax[0].set_title('Temperature in the heat exchanger 1')

ax[1].plot(phi_middle_matrix)
ax[1].set_title('Neutron Flux with time in the middle')

ax[2].plot(temperature_fuel_middle_matrix)
ax[2].set_title('Temperature in the core with time in the middle')
plt.show()
