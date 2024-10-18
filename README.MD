# Work stream
## Done
### update: use `odeint` lib to solve ode and compilable version
### update: get results align what in python program
### update: write ode solver
### update: modify parameters and reasonable results
### update: use struct instead of `extern` in parameters
### update: write a solver to replace `SparseLU`
## Pending
### update: rearrange the structs in params, all initial conditions should be in another struct (input), seperated with computing used parameters 
### update: use array, instead of vectors (static, use `dense` instead of `eigen`, only basic operants can be used)
### update: go back to cpp version 14
### update: test more sets of export temperature from simulations of brayton cycle
### update: align all the results details with article diagrams, and reserve some dataset logs
# Theory
### Crank-Nicolson
### LU solver
### RK45
### variable time step in ode solver
### FDTD (mesh: space with crank-nicolson, time with variable step), semi-implicit

