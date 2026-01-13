%Monte Carlo simulation of the stochastic particle
%aggregation model with memory in 1D

clearvars

N=400;                 % number of agents
T = 1e6;               % number of time steps
dt = 1e-3;             % time step length
nK = 6;                % memory layers K = 1, ..., nK
nMC = 100;               % number of Monte Carlo runs
int_r = 0.025;         % interaction radius

% pre-allocate xdata
xdata = zeros(nK,nMC,N);

% run over K (number of memory layers)
for K=1:nK

    %Monte-Carlo runs
    for i=1:nMC
        %collect the x-coordinates into data
        xdata(K,i,:) = f_1D_AggrMemory(N,K,T,dt,int_r,0);
    end
end

% save data
save('data1D.mat');

