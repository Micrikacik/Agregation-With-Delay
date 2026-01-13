% Particle Aggregation with Memory in 2D
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
% x = array(2,N) - x-coordinates of all particles at the final time step
% y1 = array(K,N) - internal variables y_1 at time the final time step
% y2 = array(K,N) - internal variables y_2 at time the final time step
% vdata = array(N,vMax) - collection of datapoints (Nrho, |y^1|)

function [x, y1, y2, vdata] = f_2D_AggrMemory(N,K,T,dt,int_r,WhichPlot)

alpha = ones(1,K);        % vector of alpha parameters

% Initial condition
x = rand(2,N);
y1 = zeros(K,N);
y2 = zeros(K,N);

%preallocation
y_updt1 = y1;
y_updt2 = y2;
vdataCount = zeros(N);
vMax = min(T,100);
vdata = zeros(N,vMax);
WhichV = ceil(T/vMax);

% solve for t=1:T
for t=1:T
    
    %plot
    if (~mod(t,WhichPlot))

        scatter(x(1,:),x(2,:),'o'); 
        axis([0 1 0 1]);
        getframe;        
    end


    %calculate the distances (over the torus) between the agents
    dist = torusDistances(x);

    %normalized local density
    Nrho = sum(dist<=int_r);
    rho = Nrho / (pi*int_r*int_r*N);
        
    %diffusivity
    F = exp(-2*rho);
                
    %calculate the y-update for the first K-1 layers of internal variables
    for k=1:K-1
        y_updt1(k,:) = dt*( -alpha(k)*y1(k,:) + y1(k+1,:) );
        y_updt2(k,:) = dt*( -alpha(k)*y2(k,:) + y2(k+1,:) );
    end

    %calculate the y-update for the K-th layer of internal variables
    y_updt1(K,:) = -dt*alpha(K)*y1(K,:) + sqrt(dt)*F.*randn(1,N);
    y_updt2(K,:) = -dt*alpha(K)*y2(K,:) + sqrt(dt)*F.*randn(1,N);

    
    %make one time step
    y1 = y1 + y_updt1;
    y2 = y2 + y_updt2;
    x(1,:) = x(1,:) + dt*y1(1,:);
    x(2,:) = x(2,:) + dt*y2(1,:);
    
    %periodic BCs
    x = mod(x,1);


    %collect the velocity data
    if (~mod(t-1,WhichV))
        %calculate |y(1,:)|
        v = [y1(1,:)' y2(1,:)'];
        nv = sqrt(sum(v.^2,2));
        
        %record vdata for all particles
        for i=1:N
           vdataCount(Nrho(i)) = vdataCount(Nrho(i))+1;
           vdata(Nrho(i),vdataCount(Nrho(i))) = nv(i);
        end
    end

end


end



