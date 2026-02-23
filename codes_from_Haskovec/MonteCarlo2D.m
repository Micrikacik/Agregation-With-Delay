%Monte Carlo simulation of the stochastic particle
%aggregation model with memory in 1D

clearvars

N=400;                 % number of agents
T = 1e5;               % number of time steps
dt = 1e-3;             % time step length
nK = 2;                % memory layers K = 1, ..., nK
nMC = 3;               % number of Monte Carlo runs
int_r = 0.025;         % interaction radius


% run over K (number of memory layers)
for K=1:nK
    MonteCarlo2D_K(nMC,N,K,T,dt,int_r);
end

