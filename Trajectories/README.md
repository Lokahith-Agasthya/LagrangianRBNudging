Guide to generating the trajectories

In temppart_dump.c - 
    nx  -  length of domain in x
    ny  - length of domain in y (vertical)
    scratch - 1 for initialising from scratch, 0 to restart from previous run
    nsteps - length of run (timesteps)
    nout - Frequency of field and profile outputs
    noutconfig - Frequency of particle outputs
    relax, relaxg, Tup, Tdown, gravity - Lattice Boltzmann, Rayleigh-Benard parameters
    Particle - 1 to run with particles, 0 to only solve for the flow-field
    Particle_scratch - 1 to initialise particles, 0 to restart from previous run

    The full 2D fields are written by the routine snapshottemp in the directory Tempfolder (needs to be created)
    The particle data is written in Prtcl_Veldata1

In mympi2.h - 
    Nprtcl - Number of particles
    tracer - 1 for tracer, 0 for non-tracer
    heavyprtcl - 0 if tracer, 1 if not
    taulag - Lagrangian time-scale for particle (Maxey-Riley equation)
    beta - parameter (Maxey-Riley equation)
    Prtcl_Writefull - 1 to write stretching, time-evolution terms of the M-R equation

