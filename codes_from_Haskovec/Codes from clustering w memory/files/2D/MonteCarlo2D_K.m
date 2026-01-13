%Monte Carlo simulation of the stochastic particle
%aggregation model with memory in 1D, for a specified
%number of memory layers K
%
% Input parameters:
% nMC = number of Monte Carlo runs
% N = number of particles (N>2)
% K = number of memory layers (K>0)
% T = number of time steps (T>0)
% dt = time step length
% int_r = interaction radius
%
% Output is saved into 'data2D_K.mat'

function MonteCarlo2D_K(nMC, N,K,T,dt,int_r)

% pre-allocate data
xdata = zeros(nMC,2,N);
y1data = zeros(nMC,K,N);
y2data = zeros(nMC,K,N);

%Monte-Carlo runs
for i=1:nMC

   [xdata(i,:,:), y1data(i,:,:), y2data(i,:,:), vdata(i,:,:)] = f_2D_AggrMemory(N,K,T,dt,int_r,sc,0);
   hdata(i,:) = hdata(i,:) + hd;
end


% save data
filename = sprintf('data2D_%d.mat', K);
save(filename);

end
