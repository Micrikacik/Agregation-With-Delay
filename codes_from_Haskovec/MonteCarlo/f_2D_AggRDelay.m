% Particle Aggregation with Reaction-type Delay in 2D
% with periodic boundary conditions
% and a given interaction radius
%with random initial particle positions x
%
% Input parameters:
% N = number of particles (N>2)
% tau = delay in timesteps (tau>0)
% T = number of time steps (T>0)
% dt = time step length
% int_r = interaction radius
% WhichPlot = how often plot (in timesteps); set to zero to switch off plotting
%
% Output:
% x = array(2,N) - x-coordinates of all particles at the final time step
% vdata = array(N,vMax) - collection of datapoints (Nrho(t), Nrho(t-tau))

function [x, vdata] = f_2D_AggRDelay(N,tau,T,dt,int_r,WhichPlot)

% Prepare function W
dim=2;
kappa = pi^(dim/2.) / gamma(dim / 2. + 1);  % volume of a unit d-ball
W_norm = kappa * int_r^dim;    % W_norm to normalize W
%W = @(D) (D <= int_r) ./ W_norm;    % default handle of the weight function


% Prepare initial condition
xBuf = zeros(tau,2,N);
xBuf(1,:,:) = rand(2,N);
for t=2:tau
    xBuf(t,:,:) = xBuf(t-1,:,:) + sqrt(dt)*randn(1,2,N);
end

% Prepare x
x = reshape(xBuf(tau,:,:),dim,N) + sqrt(dt)*randn(2,N);

%preallocation
vdata = zeros(N,N);

% solve for t=1:T
for t=1:T
    
    %plot
    if (~mod(t,WhichPlot))

        scatter(x(1,:),x(2,:),'o'); 
        axis([0 1 0 1]);
        getframe;        
    end

    %take xdelay from the buffer
    ttau = mod(t-1-tau,tau)+1;
    xdelay = reshape(xBuf(ttau,:,:),dim,N);

    %store present x to the buffer
    xBuf(ttau,:,:) = x;

    %calculate the distances (over the torus) from xdelay
    dist = torusDistances(xdelay);

    %local number density
    rhoN = sum(dist<int_r,2);
        
    %collect the velocity data
    if true %(~mod(t-1,WhichV))
        %calculate current distances
        pres_dist = torusDistances(x);
        
        %normalized current local density
        pres_rhoN = sum(pres_dist<int_r,2);
                
        %record vdata for all particles
        %vdata = vdata + full(sparse(pres_rhoN, rhoN, 1, size(vdata,1),
        %size(vdata,2)));
        vdata = vdata + accumarray([pres_rhoN(:), rhoN(:)], 1, size(vdata));
    end

    % Diffusivity
    theta = rhoN / (W_norm*(N-1));
    G_vals = exp(-theta);
            
    % Calculate the update
    updt = repmat(G_vals,[1,dim]).*randn(N,dim);
    
    % Make one time step
    x = x + sqrt(dt)*updt';

    % Periodic BCs
    x = mod(x,1);

end


end



