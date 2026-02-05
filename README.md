Code to produce Lagrangian trajectories of particles in 2D Rayleigh-Benard convection. The particles can be tracers or have varying Stokes numbers. The values of temperature, velocity, etc. are stored along the trajectories. 

These stored values can then be used to nudge temperature fields along the same trajectories. 

Paper reference - https://doi.org/10.1063/5.0079625
Method - Lattice-Boltzmann Method, parallelised with MPI