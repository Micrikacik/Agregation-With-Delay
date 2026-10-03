% 1D mean field model for diffusive aggregation, semi-implicit fin. diff.
% periodic boundary conditions
%clearvars;

%random seed
rng(1);

% coefficient setting:
dt = 1e-3;             % time step length
phys_tau = 1;          % "physical" delay
tau = floor(phys_tau/dt);     % delay measured in dt
int_r = getIntRad(1);           % interaction radius
a=1;                    % G(s) = exp(-a*s)
L = 1;                  % interval length

Nx = 400;               % number of gridpoints
dx = L / Nx;

T = 100;               % number of time steps

dtxx = dt / (dx^2);

WhichPlot = 100;         % How often plot; set to zero to switch off plotting


% Equidistant grid
Gx = (1:Nx)*dx;


% Pre-allocations
%rho=zeros(Nx,1);
G=zeros(Nx,Nx);
A=zeros(Nx,Nx);


% Initial condition for rho
rho = rand(Nx,1);

%Normalize such that max(rho) = 1;
rho = rho/(dx*sum(rho));

rhoBuf = zeros(tau,Nx);
rhoBuf(1,:) = rho;
for t=2:tau
    rrho = rand(Nx,1);
    rrho = rrho/(dx*sum(rrho));
    rhoBuf(t,:) = rrho;
end


%pre-calculate the distances over the periodic domain
k = 0:Nx-1;
dist = min(k, Nx-k) * dx;

%pre-calculate the kernel W, normalized such that its integral is 1
% Remark: Don't use simple 'W = (dist <  int_r)' since this gives issues
% with roundoff errors (staircasing) when int_r is a multiple of dx
tol = 100*eps(L);
w = double(dist < int_r - tol);   % strict cutoff, robust near multiples of dx
w = w / sum(w);

W = zeros(Nx,Nx);
for i = 1:Nx
    W(i,:) = circshift(w, [0, i-1]);
end

% solve for t=1:T
for t=1:T
    
    ttau = mod(t-1-tau,tau)+1;
    rhodelay = rhoBuf(ttau,:)';

    %store present rho to the buffer
    rhoBuf(ttau,:,:) = rho;

    %calculate G
    Wrho = W*rhodelay;
    F = exp(-a*Wrho);
    G = (F.^2)/2;
    
    %calculate the L^2 norm
    %E(t)=sum((rho.^2).*G);
    
    %plot
    if (~mod(t-1,WhichPlot) || t == T)
        
        % subplot(2,1,1);
        plot(Gx,rho); axis([0 L 0 max(4,max(rho+0.5))]);
        ttl = sprintf('rho; physical time = %e', (t-1)*dt);
        title(ttl);
        % subplot(2,1,2);
        % plot(Gx,Wrho,Gx,G,'--');
        % title('G, Wrho')
        getframe;
        
        %pause to see the IC
        %if (t==1) pause; end
        
        %fname = sprintf('num_mf1/num_mf1_%f.eps',phys_time);
        %saveas(gcf, fname);
        %fname = sprintf('num_mf1/num_mf1_%f.fig',phys_time);
        %saveas(gcf, fname);        
        
    end
        
    %compose the matrix for semi-implicit finite differences
    A=diag(1+dtxx*G) - 0.5*dtxx*diag(G(2:Nx),1) - 0.5*dtxx*diag(G(1:Nx-1),-1);
    A(1,Nx) = -0.5*dtxx*G(Nx);
    A(Nx,1) = -0.5*dtxx*G(1);

    %make the step
    rho = A\rho;
    
    %print the total mass to check mass conservation
    %sum(rho)
        
end
