%script for running the 1D mean field model for particle aggregation

%parameter values
L = 1;         % domain length
Nx = 300;      % number of gridpoints
T = 2e4;       % number of time steps
dt = 1e-3;     % time step length
int_r = 0.0501;   % interaction radius

%initial datum
IC = @(x) 1+0.1*sin(2*pi*L*x.*(2-x)).^2;

%run the simulation
[rho, Gx] = mean_field_1D_aggregation(L,Nx,T,dt,int_r,IC);

%plot:
figure;

% Create axes
axes1 = axes;
hold(axes1,'on');

% Create plot
plot(Gx,rho(:,end),'LineWidth',2);

% Create ylabel
ylabel('$\rho$','Interpreter','latex');

% Create xlabel
xlabel({'$x$'},'Interpreter','latex');

% Create title
title({'\rho=\rho(x)'});

box(axes1,'on');
hold(axes1,'off');
% Set the remaining axes properties
set(axes1,'FontSize',24);

