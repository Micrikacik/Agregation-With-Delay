clear all;
clc;

% coefficient setting:
N = 500;                % number of particles
d = 2;                  % dimension of the space

int_r = 3/sqrt(2*N);      % interaction radius

function kappa = unitBallVolume(d)
    kappa = pi^(d/2.) / gamma(d / 2. + 1);
end

W_norm = unitBallVolume(d) * int_r^d;
weight_params = struct("W_norm", W_norm, "int_r", int_r);

T = 300;               % number of time steps
dt = 1e-2;               % time step length

WhichPlot = 1;          % How often plot; set to zero to switch off plotting


% Initial condition
x = rand(N,2);


% solve for t=1:T
for t=1:T
    
    %calculate the distances (over the toerus) between the agents
    D = torusDistances(x);
    
    %local density
    theta = sum(weightFunction(D,weight_params),2) / (N - 1);
    % theta = sum(dist<=int_r) / (pi*int_r^2*(N-1));
        
    %diffusivity
    G = responseFunction(theta);
            
    %calculate the update
    updt = repmat(G.*randn(N,1),d);
    

    %make one time step
    x = x + sqrt(dt)*updt;

    %periodic BCs
    x = mod(x,1);

    
    %plot
    if (~mod(t-1,WhichPlot))
        scatter(x(:,1),x(:,2),'o'); 
        axis([0 1 0 1]);
        getframe;
        %if (t==1) pause; end;
    end
end

%identify clusters:

%parameter setting for DBSCAN method
epsilon = 1.3/sqrt(2*N);
minpts = 12;

%calculate distance over the torus
dist = torusDistances(x);
        
%identify clusters
idx = dbscan(dist,epsilon,minpts,'Distance','precomputed');

gscatter(x(:,1),x(:,2),idx);
getframe;
