## 00 - converge to zero (14 at 2000 steps)
```python
# Neutronics
dt=0.2 # fixed time step
L = 172  # Length of the spatial domain, m
# L = 22.9
N = 200  # Number of spatial points
dz = L / (N - 1)    # Spatial step size
# V = 4e5
V = 1.103497 * 1e8  # cm/s           
# D = 0.390016 * 8   
D = 0.96343 * 7  # cm      
# sigma_a =0.0835   
sigma_a=1.58430*1e-2/8 # cm^-1        
nu_sigma_f = 3.33029e-2/8 # cm^-1
# nu_sigma_f = 3.33029
# sigma_f= nu_sigma_f/2.41   
beta = [0.000228, 0.000788, 0.000664, 0.000736, 0.000136, 0.000088] # Delayed neutron fractions
Beta = sum(beta)      # 0.00264
```