%runs one simulation of the stochastic particle
%aggregation model with memory in 1D

clearvars

%parameters
N=400;                  % number of agents (particles)
T = 1e6;                % number of time steps
dt = 1e-3;              % time step length
K = 2;                  % number of memory layers
int_r = 0.025;          % interaction radius
WhichPlot = 1e3;        % how often plot

%run the simulation and collect the x and y variables of the final timestep
[x,y] = f_1D_AggrMemory(N,K,T,dt,int_r,WhichPlot);
