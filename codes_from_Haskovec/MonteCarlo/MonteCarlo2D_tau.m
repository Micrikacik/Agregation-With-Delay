%Monte Carlo simulation of the stochastic particle
%aggregation model with memory in 1D, for a specified
%number of memory layers K
%
% Input parameters:
% nMC = number of Monte Carlo runs
% N = number of particles (N>2)
% tau = delay in timesteps (tau>0)
% T = number of time steps (T>0)
% dt = time step length
% int_r = interaction radius
%
% Output is saved into 'data2D_K.mat'

function MonteCarlo2D_tau(nMC, N,tau,T,dt,int_r)

% pre-allocate data
xdata = zeros(nMC,2,N);
vdata = zeros(nMC,N,5000);

%Monte-Carlo runs
for i=1:nMC
   [xdata(i,:,:), vd] = f_2D_AggRDelay(N,tau,T,dt,int_r,0);
   vdata(i,1:N,1:size(vd,2)) = vd;
end


% save data
filename = sprintf('data2D_%d.mat', tau);
save(filename);

end
