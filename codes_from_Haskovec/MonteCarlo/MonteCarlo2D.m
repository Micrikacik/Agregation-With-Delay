%Monte Carlo simulation of the stochastic particle
%aggregation model with memory in 1D

clearvars

N=400;                 % number of agents
T = 1e7;               % number of time steps
dt = 1e-3;             % time step length
Ntau = 10;             % number of values of tau
taustep = 30;
nMC = 100;             % number of Monte Carlo runs
int_r = 0.05;          % interaction radius


% run over tau
parfor i=1:Ntau
    tau = i*taustep;
    MonteCarlo2D_tau(nMC,N,tau,T,dt,int_r);
end

