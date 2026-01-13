% 1D mean field model for particle aggregation
% using semi-implicit finite difference method
% with periodic boundary conditions
%
% Input:
% L = domain (interval) length
% Nx = number of grid points
% T = number of time steps
% dt = time step length
% int_r = interaction radius
% IC = initial condition (function defined on [0,L] with periodic BC)
%
% Output:
% rho = array(Nx,T) - the solution over T timesteps
% Gx = array(Nx,1) - equidistant grid

function [rho,Gx] = f_mean_field_1D_aggregation(L,Nx,T,dt,int_r,IC)


% convenience variables
dx = L / Nx;
dtxx = dt / (dx^2);

% Equidistant grid
Gx = (1:Nx)*dx;


% Pre-allocations
rho=zeros(Nx,T);
A=zeros(Nx,Nx);
W=zeros(Nx,1);

%pre-calculate the distances over the periodic domain
dist=squareform(pdist(Gx'));
dist=min(dist,1-dist);

%pre-calculate the kernel W
W = (dist <= int_r);
W = W / (sum(W(1,:)));


% Initial condition for rho
%rhoI = 1+0.01*randn(Nx,1);
%rhoI= 1+0.01*sin(2*pi*Gx');
rhoI = IC(Gx');

%Normalize to unit mass
rhoI = rhoI/(dx*sum(rhoI));

rho(:,1) = rhoI;

% solve for t=1:T
for t=1:T

    %calculate f
    Wrho = W*rho(:,t);
    G = exp(-2*Wrho);    
        
    %compose the matrix for semi-implicit finite differences
    A=diag(1+dtxx*G) - 0.5*dtxx*diag(G(2:Nx),1) - 0.5*dtxx*diag(G(1:Nx-1),-1);
    A(1,Nx) = -0.5*dtxx*G(Nx);
    A(Nx,1) = -0.5*dtxx*G(1);

    %make one time step
    rho(:,t+1) = A\rho(:,t);
        
end

end

