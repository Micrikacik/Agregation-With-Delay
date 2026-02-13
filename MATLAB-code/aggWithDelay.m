function [xRec] = aggWithDelay(expParams)

% Runs the discrete simulation of agregation with a constant delay.
%
%---------------------------------------------------------------------------
%
% INPUT :
%   expParams - struct, which can contain following fields, if an important
%       field is missing a default value is used. If no input is required,
%       set to '{}'.
%
%       POSSIBLE FIELDS:
%       rngSeed (integer) - rng seed to replicate experiments.
%       stepDelay (nonnegative integer) - number of steps used to delay the simulation.
%       dims (float ROW vector) - dimensions of the simulation, i.e., dimensions
%           of the box, in which the agents move.
%           Thic row vector does NOT influence the parameter d, in fact, dims
%           will be either truncated or filled up with 1, to match its
%           length with d.
%       x0 (float matrix) - matrix of initial positions.
%           x(i,:) - position (float vector) in [0,1]^d of the i-th agent.
%           Alternatively, instead of x0, set the dimension d and number of
%           agent N.
%       xInitHist (float matrix) - initial history of the matrix of positions used in
%           calculation of the first few iterations.
%           x(:,:,i) - position matrix i steps into the past, 
%           where 1 <= i <= stepDelay.
%       T (positive integer) - number of time steps.
%       dt (positive float) - time step length.
%       W - handle of the weight function .
%           Should take the distance matrix D from torusDistance.m and return 
%           a weight matrix corresponding to D.
%       G - handle of the response function.
%           Should take the average density vector theta and return a response 
%           vector corresponding to theta.
%       D - handle of the distance function.
%           Should take the matrix of positions x and return a distance matrix
%           which will be used as an input to W.
%       delayType (string) - type of the delay.
%           Should be on of the following strings:
%               "Reaction"
%               "Transmission"
%               "Memory"
%               "None"
%       stepRecMod (positive integer) - simulation records the t-th step if the reminder of t 
%           divided by stepRecMod is zero. Last state is always recorded.
%           If stepRecMod = 0, then only initial and last states are recorded.
%           If stepRecMod = -1, then only the last state is recorded. 
%       waitForConf (logical) - if true, then wait for user to start the simulation.
%       stepPlotMod (positive integer) - simulation plots the t-th step if the reminder of t 
%           divided by stepRecMod is zero. The plots are shown in movie-like sense. 
%           If stepPlotMod = -1, then only the last state is plotted.
%           If stepPlotMod = -2, then no states are plotted.

fprintf("----------------------------------\n\n")

fprintf("Initializing the experiment: Agregation with delay.\n\n")

function result =  IsInteger(x)
    result = isnumeric(x) && x == floor(x);
end

% Set or initialize experiment parameters

% Setting rng seed for replicable experiment
if ~isfield(expParams,"rngSeed") || ~IsInteger(expParams.rngSeed) || expParams.rngSeed < 0
    fprintf("Either no or wrong value for the rng seed 'rngSeed'.\n")
    fprintf("Randomness is uncontrolled.\n\n")
else
    rng(expParams.rngSeed)
    fprintf("Rng seed: %i.\n\n", expParams.rngSeed)
end

% Setting step delay
if ~isfield(expParams,"stepDelay") || ~IsInteger(expParams.stepDelay) || expParams.stepDelay < 0
    fprintf("Either no or wrong value for the step delay 'stepDelay'.\n")
    stepDelay = int16(5);
    fprintf("Setting step delay to %i.\n\n", stepDelay)
else
    stepDelay = expParams.stepDelay;
    fprintf("Step delay: %i.\n\n", stepDelay)
end

% Setting initial positions
if ~isfield(expParams,"x0") || ~isfloat(expParams.x0) || isempty(expParams.x0)
    fprintf("Either no or wrong value for the matrix of initial positions 'x0'.\n")
    fprintf("Searching for the input values for 'N' and 'd'.\n")
    fprintf('   |\n')
    if ~isfield(expParams,"N") || ~IsInteger(expParams.N)
        fprintf("   Either no or wrong value for the number of agents 'N'.\n")
        N = 400;                % default number of particles
        fprintf("   Setting N = %i.\n", N);
    else
        N = single(expParams.N);
        fprintf("   N = %i.\n", N)
    end
    if ~isfield(expParams,"d") || ~IsInteger(expParams.d)
        fprintf("   Either no or wrong value for the dimension 'd'.\n")
        d = 2;                  % default dimension of the space
        fprintf("   Setting d = %i.\n", d);
    else
        d = single(expParams.d);
        fprintf("   d = %i.\n", d)
    end
    fprintf('   |\n')
    fprintf("Initializing random experiment with N = %i, d = %i.\n\n", N, d);
    x = [];   % we will set the initial position later as a random matrix
else
    x = expParams.x0;
    N = size(expParams.x0,1);
    d = size(expParams.x0,2);
    fprintf("Matrix of initial positions accepted, N = %i, d = %i.\n\n", N, d)
end

% Setting dimensions
if ~isfield(expParams,"dims") || ~isfloat(expParams.dims) || isempty(expParams.dims) || ...
        (size(expParams.dims, 1) ~= 1 && size(expParams.dims, 2) ~= 1)
    fprintf("Either no or wrong value for the vector of dimensions 'dims'.\n")
    dims = ones(1,d);  % default number of time steps
    fprintf("Setting all dimensions to 1.\n\n")
else 
    dims = expParams.dims;
    % Making dims to be a row vector
    if size(dims,2) == 1
        dims = dims.';
    end
    % Resizing dims to match d
    l = length(dims);
    if l > d
        dims = dims(1:d);
    else 
        if l < d
            dims(l+1:d) = ones(d-l,1);
        end
    end
    fprintf("dims = %.3d.\n\n", dims)
end

% Setting random initial positions
if isempty(x)
    x = rand(N, d) * diag(dims);
end

% Setting initial history of positions
if ~isfield(expParams,"xInitHist") || ~isfloat(expParams.xInitHist) || ...
        ~isequal(size(expParams.xInitHist),[N,d,stepDelay])
    fprintf("Either no or wrong value for the matrix of initial history of positions 'xInitHist'.\n")
    xHist = repmat(x,[1,1,stepDelay]);  % default number of time steps
    fprintf("Initializing experiment with constant initial history.\n\n")
else
    xHist = expParams.xInitHist;
    fprintf("Initial history of the matrix of position accepted.\n\n")
end

% Setting step count
if ~isfield(expParams,"T") || ~IsInteger(expParams.T) || expParams.T <= 0
    fprintf("Either no or wrong value for the number of time steps T.\n")
    T = 300;                    % default number of time steps
    fprintf("Setting T = %i.\n\n", T)
else
    T = expParams.T;
    fprintf("T = %i.\n\n", T)
end

% Setting time step length
if ~isfield(expParams,"dt") || ~isfloat(expParams.dt) || expParams.dt <= 0
    fprintf("Either no or wrong value for the time step length dt.\n")
    dt = 1e-2;                  % default time step length
    fprintf("Setting dt = %.3d\n\n", dt)
else
    dt = expParams.dt;
    fprintf("dt = %.3d\n\n", dt)
end

% Setting the handle of the weight function
if ~isfield(expParams,"W") || ~isa(expParams.W,"function_handle")
    fprintf("Either no or wrong value for the handle of the weight function W.\n")
    %int_r = 3/sqrt(2*N);        % interaction radius
    int_r = 1/10.;        % interaction radius
    kappa = pi^(d/2.) / gamma(d / 2. + 1);  % volume of a unit d-ball
    W_norm = kappa * int_r^d;   % W_norm to normalize W
    W = @(D) (D <= int_r) ./ W_norm;    % default handle of the weight function
    fprintf("Setting W = %s.\n\n", func2str(W))
else
    W = expParams.W;
    fprintf("W = %s.\n\n", func2str(W))
end

% Setting the handle of the response function
if ~isfield(expParams,"G") || ~isa(expParams.G,"function_handle")
    fprintf("Either no or wrong value for the handle of the response function G.\n")
    G = @(theta) exp(-theta);   % default handle of the response function
    fprintf("Setting G = %s.\n\n", func2str(G))
else
    G = expParams.G;
    fprintf("G = %s.\n\n", func2str(G))
end

% Setting the delay type
if ~isfield(expParams,"delayType") || ~isstring(expParams.delayType)
    fprintf("Either no or wrong value for the delay type.\n")
    delayType = "Reaction";   % default handle of the distance function
    fprintf("Setting delay type to %s.\n\n", delayType)
else
    delayType = expParams.delayType;
    fprintf("Delay type: %s.\n\n", delayType)
end

% Forcing no delay
if stepDelay == 0
    delayType = "None";
end

% Setting step record mod
if ~isfield(expParams,"stepRecMod") || ~IsInteger(expParams.stepRecMod) || ... 
        (expParams.stepRecMod <= 0 && expParams.stepRecMod ~= 0 && expParams.stepRecMod ~= -1)
    fprintf("Either no or wrong value for the record mod.\n")
    stepRecMod = -1;  % default step record mod
    fprintf("Setting step record mod to %i.\n\n", stepRecMod)
else
    stepRecMod = expParams.stepRecMod;
    fprintf("Step record mod: %i.\n\n", stepRecMod)
end

% Setting step plot mod
if ~isfield(expParams,"stepPlotMod") || ~IsInteger(expParams.stepPlotMod) || ...
        ((expParams.stepPlotMod <= 0) && expParams.stepPlotMod ~= -1 && expParams.stepPlotMod ~= -2)
    fprintf("Either no or wrong value for the step plot mod.\n")
    stepPlotMod = 1;  % default step plot mod
    fprintf("Setting step plot mod to %i.\n\n", stepPlotMod)
else
    stepPlotMod = expParams.stepPlotMod;
    fprintf("Step plot mod: %i.\n\n", stepPlotMod)
end


% Set auxiliary variables
histCoeff = stepDelay;

% Set output variables
% Decide the size of xRec
recCount = 0;
if stepRecMod > 0
    recCount = idivide(T, stepRecMod);
    if mod(T, stepRecMod) == 0
        recCount = recCount - 1;
    end
end

% Do we record initial state
if stepRecMod >= 0
    xRec = zeros([N,d,recCount + 1]);
    xRec(:,:,1) = x;
    recIndex = 2;
else
    xRec = zeros([N,d,recCount + 1]);
    recIndex = 1;
end


if ~isfield(expParams, "waitForConf") || expParams.waitForConf == true
    fprintf("----------------------------------\n\n")
    fprintf("Press space to start the simulation.\n\n")
    pause
end


fprintf("----------------------------------\n\n")
fprintf("Starting the simulation.\n\n")

% Saving maximal values of W to normalize 1D graphs
W_max = W(0);

volume = prod(dims);

% Simulate for t=1:T
for t=1:T
    
    % Calculate the distances (over the torus) between the agents
    % Delay is included by taking the last x from x_history
    switch delayType
        case "Reaction"
            D = torusDistances(xHist(:,:,histCoeff),xHist(:,:,histCoeff),dims);
        case "Transmission"
            D = torusDistances(x,xHist(:,:,histCoeff),dims);
        case "Memory"
            D = torusDistances(xHist(:,:,histCoeff),x,dims);
        case "None"
            D = torusDistances(x,x,dims);
        otherwise
            error('Invalid delay type.');
    end
    
    % Local density
    theta = sum(W(D),2) / (N - 1) * volume; 
        
    % Diffusivity
    G_vals = G(theta);
            
    % Calculate the update
    updt = repmat(G_vals,[1,d]).*randn(N,d);
    
    % Make one time step
    x = x + sqrt(dt)*updt;

    % Periodic BCs
    x = mod(x,dims);

    % Update history of x
    if delayType ~= "None"
        xHist(:,:,histCoeff) = x;
        histCoeff = histCoeff - 1;
        if histCoeff <= 0
            histCoeff = stepDelay;
        end
    end

    % Record simulation state
    if stepRecMod ~= -1 && mod(t,stepRecMod) == 0
        xRec(:,:,recIndex) = x;
        recIndex = recIndex + 1;
    end

    % Plot
    if stepPlotMod > 0 && (~mod(t,stepPlotMod))
        plotSimState()
    end
end

function plotSimState()
    switch d
        case 1
            scatter(x(:,1),theta ./ (W_max * volume),'o'); 
            axis([0 dims(1) 0 1]);
            getframe;
        case 2
            scatter(x(:,1),x(:,2),'o'); 
            axis([0 dims(1) 0 dims(2)]);
            getframe;
        case 3
            scatter3(x(:,1),x(:,2),x(:,3),'o'); 
            axis([0 dims(1) 0 dims(2) 0 dims(3)]);
            getframe;
    end
end

% Record final simulation state
xRec(:,:,end) = x;

fprintf("----------------------------------\n\n")
fprintf("Simulation finished.\n\n")

if stepPlotMod ~= -2
    % We need to update theta for the last plot
    switch delayType
        case "Reaction"
            D = torusDistances(xHist(:,:,histCoeff),xHist(:,:,histCoeff),dims);
        case "Transmission"
            D = torusDistances(x,xHist(:,:,histCoeff),dims);
        case "Memory"
            D = torusDistances(xHist(:,:,histCoeff),x,dims);
        case "None"
            D = torusDistances(x,x,dims);
        otherwise
            error('Invalid delay type.');
    end
    
    theta = sum(W(D),2) / (N - 1);

    fprintf("----------------------------------\n\n")
    fprintf("Plotting agregation groups.\n\n")
    plotAgg(x,theta ./ (W_max * volume), dims)
end

end