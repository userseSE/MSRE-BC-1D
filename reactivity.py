from parameters import *
import numpy as np

def reactivity(temperature_fuel, temperature_graphite, step, time_span, rho_insertion):
    
    rho_0=sum(beta)
    for i in range(6):
        rho_0=rho_0-beta[i]/(1+(1/(lambda_i[i]*tau_c))*(1-np.exp(-lambda_i[i]*tau_l)))
    # print(rho_0)
    
    # SOURCE INSERTION
    # REACTIVITY INSERTION
    # No reactivity insertion
    # simtime = 100
    # reacttime = [0, 2500]
    # Periodic 60 PCM for 50 seconds
    # simtime = 500;
    # periodic = [0, 0; 50, 6e-4; 100, 0; 150, -6e-4; 200, 0; 250, 6e-4; 300, 0; 350, -6e-4; 400, 0]; 
    # reactdata = periodic(:,2);
    # reacttime = periodic(:,1);
    # Step up 60 pcm 
    # simtime = 1000;
    # reactdata = [0 6e-3];
    # reacttime = [0 300];
    # % Step down -60 pcm for 10 sec
    # simtime = 100;
    # reactdata = [0 -6e-4];
    # reacttime = [0 50];
    # % Pulse 600 pcm for 0.1 sec
    # simtime = 30;
    # reactdata = [0 6e-3 0];
    # reacttime = [0 10 10.1];
    # react = timeseries(reactdata,reacttime);
    
    reactdata = [rho_init*1e-5, rho_insertion*1e-5]
    if step < time_span/2:
        react = reactdata[0]
    else:
        react = reactdata[1]
        
    rho_feedback=(initialS-temperature_fuel)*alpha_f+(initialG-temperature_graphite)*alpha_g
    
    rho=rho_0+rho_feedback+react
    
    return rho