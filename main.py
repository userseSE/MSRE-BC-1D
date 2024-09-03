import numpy as np
import os
# from scipy.integrate import solve_ivp
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
from power_plant import power_plant_temp

time_span = 5000

# reactivity insertion rod
rho_insertion = (rho_init) * np.ones(N)  # pcm

# initialization
rho=rho_init*1e-5
# initial variables
y_n=np.zeros((7*N, 1))
q_prime=np.zeros((N, 1))
y_th=np.zeros((2*N, 1))
y_hx1=np.zeros((2*Nx, 1))
y_hx2=np.zeros((2*Nx, 1))
Tss_HX2_0=0
Ts_HX1_0=0
Tss_HX1_0=0
Tsss_pp_0=0
buffer_hx_c=[]
buffer_c_hx=[]
buffer_r_hx=[]
buffer_hx_r=[]
buffer_r_pp=[]
buffer_pp_r=[]

# initial plotting matrices
rho_matrix = np.zeros((time_span, 1))
phi_middle_matrix = np.zeros((time_span, 1))
ci_middle_matrix = np.zeros((time_span, 6))
temperature_fuel_middle_matrix = np.zeros((time_span, 1))
rho_dt_matrix = np.zeros((time_span, 1))

for step in range (time_span):
    print (f"Time step: {step}")
    # int(time_span/dt)
     
    y_n, q_prime = neutronics(y_n[:,-1], rho, step)
    # q_prime=q_prime*sigma_f/(3.12*1e10)
    # print("test1") 
    phi = y_n[:N ,-1].T
    ci = y_n[N:,-1].T
    phi_middle_matrix[step] = phi[int(N/2)]
    for i in range(6):
        ci_middle_matrix[step, i] = ci[int((i*N+(i+1)*N)/2)]
    # print(phi.shape)
    # print(ci.shape)
    # print(y_n.shape)
    
    # tau_hx_c
    Ts_core_0=transport_delay(Ts_HX1_0, tau_hx_c, Ts_in, buffer_hx_c, step)
    # print("size of buffer_hx_c: " + str(len(buffer_hx_c)))
    y_th = thermal_hydraulics(y_th, q_prime, Ts_core_0, step)
    temperature_fuel = y_th[:N, -1].T
    print("avg of fuel: "+str(np.sum(temperature_fuel/N)))
    temperature_graphite = y_th[N:, -1].T
    print("avg of graphite: "+str(np.sum(temperature_graphite/N)))
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
    # tau_pp_r
    Tsss_HX2_0=transport_delay(Tsss_pp_0, tau_pp_r, Tsss_in, buffer_pp_r, step)
    # print("size of buffer_hx_r: " + str(len(buffer_hx_r)))
    y_hx2=HX2(y_hx2, Tss_HX2_L, Tsss_HX2_0, step)
    Tss_HX2= y_hx2[:Nx,-1]
    Tsss_HX2= y_hx2[Nx:,-1]
    Tss_HX2_0=Tss_HX2[0]
    Tsss_HX2_L=Tsss_HX2[-1]
    
    # tau_r_pp
    Tsss_pp_L=transport_delay(Tsss_HX2_L, tau_r_pp, Tsss_out, buffer_r_pp, step)
    y_pp=power_plant_temp(Tsss_pp_L, step)
    Tsss_pp_0=y_pp
    
    rho=reactivity(temperature_fuel, temperature_graphite, step, time_span, rho_insertion)
    rho_matrix[step]=rho[int(N/2)]
    if step>0:
        rho_dt_matrix[step]=rho[int(N/2)]-rho_matrix[step-1]

# plot phi to axis N
z = np.linspace(0, L, N)
os.chdir('plotting')

fig, ax = plt.subplots(2, 2, figsize=(14, 6))
ax[0, 0].plot(z, phi)
ax[0, 0].set_title('Neutron Flux')

for i in range(6):
    ax[0, 1].plot(z, ci[i*N:(i+1)*N], label=f'Ci{i+1}')
    ax[0, 1].legend()
ax[0, 1].set_title('Delayed Neutron Precursors')

ax[1, 0].plot(rho_matrix * 1e5, label='Reactivity at middle with time(pcm)')
ax[1, 0].set_title('Reactivity')

ax[1, 1].plot(z, temperature_fuel, label='Fuel')
ax[1, 1].plot(z, temperature_graphite, label='Graphite')
ax[1, 1].legend()
ax[1, 1].set_title('Temperature in the core')

plt.tight_layout()
plt.savefig('neutron_flux_and_precursors.png')  # Save the figure
plt.close(fig)  # Close the figure to free memory

# Second plot: Temperature in the heat exchangers
fig, ax = plt.subplots(2, 2, figsize=(14, 6))
ax[0, 0].plot(z, Ts_HX1, label='Salt')
ax[0, 0].plot(z, Tss_HX1, label='Coolant')
ax[0, 0].legend()
ax[0, 0].set_title('Temperature in the heat exchanger 1')

ax[0, 1].plot(z, Tss_HX2, label='Tss')
ax[0, 1].plot(z, Tsss_HX2, label='Tsss')
ax[0, 1].legend()
ax[0, 1].set_title('Temperature in the heat exchanger 2')

ax[1, 0].plot(phi_middle_matrix)
ax[1, 0].set_title('Neutron Flux with time in the middle')

ax[1, 1].plot(temperature_fuel_middle_matrix)
ax[1, 1].set_title('Temperature in the core with time in the middle')

plt.tight_layout()
plt.savefig('temperature_in_heat_exchangers.png')  # Save the figure
plt.close(fig)  # Close the figure to free memory

# Third plot: Delayed Neutron Precursors with time in the middle
fig, ax = plt.subplots(figsize=(14, 6))
for i in range(6):
    ax.plot(ci_middle_matrix[:,i], label=f'Ci{i+1}')
    ax.legend()
ax.set_title('Delayed Neutron Precursors with time in the middle')

plt.tight_layout()
plt.savefig('precursors_with_time_middle.png')  # Save the figure
plt.close(fig)  # Close the figure

fig, ax = plt.subplots(2, 2, figsize=(14, 6))
ax[0, 0].plot(rho_dt_matrix * 1e5, label='Reactivity change (pcm)')
ax[0, 0].set_title('Reactivity_dt')

ax[0, 1].plot(z, rho * 1e5, label='Reactivity')
ax[0, 1].set_title('Reactivity')

plt.tight_layout()
plt.savefig('rho_dt.png')  # Save the figure
plt.close(fig)  # Close the figure