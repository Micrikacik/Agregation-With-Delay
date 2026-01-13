% Particle Aggregation with Memory in 1D
% with periodic boundary conditions
% and a given interaction radius
%with random initial particle positions x
%
% Input parameters:
% N = number of particles (N>2)
% K = number of memory layers (K>0)
% T = number of time steps (T>0)
% dt = time step length
% int_r = interaction radius
% WhichPlot = how often plot (in timesteps); set to zero to switch off plotting
%
% Output:
% x = array(1,N) - x-coordinates of all particles at the final time step
% y = array(K,N) - internal variables y at time the final time step

function [x,y] = f_1D_AggrMemory(N,K,T,dt,int_r,WhichPlot)

alpha = ones(1,K);        % vector of alpha parameters

% Initial condition
x = sort(rand(1,N));
y = zeros(K,N);

%preallocation
y_updt = y;


% solve for t=1:T
for t=1:T
    
    %calculate the SQUARED distances between the agents with periodic BC
    di = pdist(x');
    di = min(di,1-di);
    dist = squareform(di);

    %local density
    rho = sum(dist<=int_r) / (int_r*N);
        
    %diffusivity
    F = exp(-rho);
                
    %calculate the y-update for the first K-1 layers of internal variables
    for k=1:K-1
        y_updt(k,:) = dt*( -alpha(k)*y(k,:) + y(k+1,:) );
    end

    
    %calculate the y-update for the K-th layer of internal variables
    y_updt(K,:) = -dt*alpha(K)*y(K,:) + sqrt(dt)*F.*randn(1,N);

    
    %make one time step
    y = y + y_updt;
    x = x + dt*y(1,:);
    
    %periodic BCs
    x = mod(x,1);


    %plot
    if (~mod(t,WhichPlot))
        plot(x,rho,'o');
        getframe;                
    end
    
end

end
