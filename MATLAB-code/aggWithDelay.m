function [xRec, thetaRec, xHist, rngSetts] = aggWithDelay(expParams)

% Runs the discrete simulation of agregation with a constant delay.
%
%---------------------------------------------------------------------------
%
% INPUT:
%   expParams - struct, which can contain following fields, if an important
%       field is missing, default value is used. If no input is required,
%       set to '{}'.
%
%       POSSIBLE FIELDS:
%
%       RNG:
%       rngSeed (integer) - rng seed to replicate experiments.
%       rngSetts (struct) - struct returned by the rng function, containing
%           random generator settings. This struct is also an output of
%           this function, and is created after the simulation is done.
%           This field's main goal is to enable user to continue in an
%           experiment which already finished, by using the values it
%           returned as initial conditions and random generator settings.
%
%       SPACE & INITIAL CONDITIONS:
%       x0 (float matrix) - matrix of initial positions.
%           x(i,:) - position (float vector) in [0,1]^d of the i-th agent.
%           Alternatively, instead of x0, set the dimension d and number of
%           agents N.
%       dims (float ROW vector) - dimensions of the simulation, i.e., dimensions
%           of the box, in which the agents move.
%           Thic row vector does NOT influence the parameter d, in fact, dims
%           will be either truncated or filled up with 1, to match its
%           length with d.
%
%       TIME & DELAY:
%       T (positive integer) - number of time steps.
%       dt (positive float) - time step length.
%       delayType (string) - type of the delay.
%           Should be one of the following strings:
%               "Reaction"
%               "Transmission"
%               "Memory"
%               "None"
%       stepDelay (nonnegative integer) - number of steps used to delay the simulation.
%       xInitHist (float matrix) - initial history of the matrix of positions used in
%           calculation of the first few iterations.
%           x(:,:,i) - position matrix i steps into the past, 
%           where 1 <= i <= stepDelay.
%
%       MODEL:
%       intRad (positive float) - radius of interactions between agents.
%       W (function handle) - handle of the weight function. If a change of
%           the interaction radius is desired, use 'intRad' field.
%           Should take the distance matrix D from torusDistance.m and return 
%           a weight matrix corresponding to D.
%       G (function handle) - handle of the response function.
%           Should take the average density vector theta and return a response 
%           vector corresponding to theta.
%       boundConds (string) - type of the boundary conditions.
%           Should be one of the following strings:
%               "Periodic"
%               "Ricochet"
%               "Reflective"
%
%       USER:
%       waitForConf (logical) - if true, then wait for user to start the simulation.
%       stepPlotMod (positive integer) - simulation plots the t-th step if the reminder of t 
%           divided by stepRecMod is zero. The plots are shown in movie-like sense. 
%           If stepPlotMod = -1, then only the last state is plotted.
%           If stepPlotMod = -2, then no states are plotted.
%       stepRecMod (positive integer) - simulation records the t-th matrix 
%           of positions if the reminder of t divided by stepRecMod is zero. 
%           Last matrix is always recorded.
%           If stepRecMod = 0, then all matrices (including x0 - initial
%           matrix) are recorded.
%           If stepRecMod = -1, then only the last matrix is recorded.
%       thetaRecMod (positive integer) - simulation records the t-th local
%           density vector (theta) and its delayed value if the reminder of t 
%           divided by stepRecMod is zero.e 
%           Last vector is always recorded.
%           If thetaRecMod = 0, then all vectors (including the one for x0) are recorded.
%           If thetaRecMod = -1, then only the last vector is recorded.
%           If no value is specified, then the value of 'stepRecMod' is
%           used instead.
%       expTitle (string) - title to be printed before the experiment begins
%
% OUTPUT:
%   xRec (float matrix) - 3 dimensional matrix of all recorded matrices of
%       positions. Its dimensions are [N,d,count], where count is the final
%       count of recorded matrices. 
%       xRec(:,:,end) is always the matrix of positions at the end of the simulation. 
%       If input stepRecMod = 0, then xRec(:,:,1) is the matrix of initial
%       positions.
%   thetaRec (float matrix) - 3 dimensional matrix of all recorded local density vectors (theta). 
%       Its dimensions are [N,2,count], where count is the final
%       count of recorded matrices.
%       thetaRec(:,1,i) is the density vector after i-th step of the simulation
%       calculated with NO delay.
%       thetaRec(:,2,i) is the density vector after i-th step of the simulation
%       calculated WITH delay, based on the delay type.
%       thetaRec(:,:,end) is always the density vector at the end of the simulation. 
%       If input thetaRecMod = 0, then thetaRec(:,:,1) is the initial
%       density vector.
%   xHist (float matrix) - 3 dimensional matrix of the last history of
%       positions, which still affect the following simulation steps.
%       Its dimensions are [N,d,stepDelay]. It can be directly plugged in
%       as an input to this function to continue in the experiment.
%   rngSetts (struct) - struct returned by the rng function, containing
%       random generator settings right after the simulation have finished.
%       It can be directly plugged in as an input to this function to 
%       continue in the experiment.


%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------


fprintf("----------------------------------\n\n")

fprintf("Initializing the experiment: Agregation with delay.\n\n")

addpath("auxiliaryFunctions\")

fprintf("----------------------------------\n\n")

function result =  IsInteger(x)
    result = isnumeric(x) && x == floor(x);
end

% Set or initialize experiment parameters

%---------------------------------RNG---------------------------------------

% Random generator settings for continuation of experiment
if ~isfield(expParams,"rngSetts") || ~isstruct(expParams.rngSetts)
    % We do not have all settings, so we check just the seed
    % Rng seed for replicable experiment
    if ~isfield(expParams,"rngSeed") || ~IsInteger(expParams.rngSeed) || expParams.rngSeed < 0
        fprintf("Either no or wrong value for the rng seed 'rngSeed'.\n")
        fprintf("Randomness is uncontrolled.\n\n")
    else
        rng(expParams.rngSeed)
        fprintf("Rng seed: %i.\n\n", expParams.rngSeed)
    end
else
    rng(expParams.rngSetts)
    fprintf("Random generator settings accepted.\n\n")
end

%-----------------------SPACE-&-INITIAL-CONDITIONS------------------------

% Initial positions
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
    N = size(x,1);
    d = size(x,2);
    fprintf("Matrix of initial positions accepted, N = %i, d = %i.\n\n", N, d)
end

% Dimensions
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

%------------------------------TIME-&-DELAY---------------------------------

% Step count
if ~isfield(expParams,"T") || ~IsInteger(expParams.T) || expParams.T < 0
    fprintf("Either no or wrong value for the number of time steps 'T'.\n")
    T = 1000;                    % default number of time steps
    fprintf("Setting T = %i.\n\n", T)
else
    T = expParams.T;
    fprintf("T = %i.\n\n", T)
end

% Time step length
if ~isfield(expParams,"dt") || ~isfloat(expParams.dt) || expParams.dt <= 0
    fprintf("Either no or wrong value for the time step length 'dt'.\n")
    dt = 1e-3;                  % default time step length
    fprintf("Setting dt = %.3d\n\n", dt)
else
    dt = expParams.dt;
    fprintf("dt = %.3d\n\n", dt)
end

% Delay type
if ~isfield(expParams,"delayType") || ~isstring(expParams.delayType)
    fprintf("Either no or wrong value for the delay type 'delayType'.\n")
    delayType = "Reaction";   % default handle of the distance function
    fprintf("Setting delay type to %s.\n\n", delayType)
else
    delayType = expParams.delayType;
    fprintf("Delay type: %s.\n\n", delayType)
end

% Step delay
if ~isfield(expParams,"stepDelay") || ~IsInteger(expParams.stepDelay) || expParams.stepDelay < 0
    fprintf("Either no or wrong value for the step delay 'stepDelay'.\n")
    stepDelay = 5;
    fprintf("Setting step delay to %i.\n\n", stepDelay)
else
    stepDelay = expParams.stepDelay;
    fprintf("Step delay: %i.\n\n", stepDelay)
end

% Forcing no delay
if stepDelay == 0
    fprintf("Zero step delay detected. Forcing delay type to be 'None'.\n\n")
    delayType = "None";
end

% Initial history of positions
if ~isfield(expParams,"xInitHist") || ~isfloat(expParams.xInitHist) || ...
        ~isequal(size(expParams.xInitHist),[N,d,stepDelay])
    fprintf("Either no or wrong value for the matrix of initial history of positions 'xInitHist'.\n")
    xHist = generInitHist(x,dt,stepDelay)  % default initial history
    fprintf("Initializing experiment with 'blind motion' initial history.\n\n")
else
    xHist = expParams.xInitHist;
    fprintf("Initial history of the matrix of positions accepted.\n\n")
end

%----------------------------------MODEL------------------------------------

% Weight function (handle)
if ~isfield(expParams,"W") || ~isa(expParams.W,"function_handle")
    fprintf("Either no or wrong value for the handle of the weight function 'W'.\n")
    fprintf("Using normed characterictic function of d-dimensional ball.\n")
    kappa_1 = 2;  % "volume" of a unit 1-ball
    kappa = pi^(d/2) / gamma(d / 2 + 1);  % volume of a unit d-ball
    fprintf("   |\n")
    if ~isfield(expParams,"intRad") || ~isnumeric(expParams.intRad) || expParams.intRad < 0
        fprintf("   Either no or wrong value for the interaction radius 'intRad'.\n")
        intRad_1 = 0.05^2 / 2 * pi;
        intRad = (kappa_1 / kappa * intRad_1)^(1/d);         % interaction radius
        fprintf("   Setting intRad = %.3d\n", intRad)
    else
        intRad = expParams.intRad;
        fprintf("   intRad = %.3d\n", intRad)
    end
    fprintf("   |\n")
    W_norm = kappa * intRad^d;   % W_norm to normalize W
    W = @(D) (D <= intRad) ./ W_norm;    % default handle of the weight function
    fprintf("Setting W = %s.\n\n", func2str(W))
else
    W = expParams.W;
    fprintf("W = %s.\n\n", func2str(W))
end

% Response function (handle)
if ~isfield(expParams,"G") || ~isa(expParams.G,"function_handle")
    fprintf("Either no or wrong value for the handle of the response function 'G'.\n")
    G = @(theta) exp(-theta);   % default handle of the response function
    fprintf("Setting G = %s.\n\n", func2str(G))
else
    G = expParams.G;
    fprintf("G = %s.\n\n", func2str(G))
end

% Boundary conditions
if ~isfield(expParams,"boundConds") || ~isstring(expParams.boundConds)
    fprintf("Either no or wrong value for the boundary conditions 'boundConds'.\n")
    boundConds = "Periodic";  % default boundary conditions
    fprintf("Setting boundary conditions to %s.\n\n", boundConds)
else
    boundConds = expParams.boundConds;
    fprintf("Boundary conditions: %s.\n\n", boundConds)
end

%----------------------------------USER-------------------------------------

% Step plot mod
if ~isfield(expParams,"stepPlotMod") || ~IsInteger(expParams.stepPlotMod) || ...
        ((expParams.stepPlotMod <= 0) && expParams.stepPlotMod ~= -1 && expParams.stepPlotMod ~= -2)
    fprintf("Either no or wrong value for the step plot mod 'stepPlotMod'.\n")
    stepPlotMod = 3;  % default step plot mod
    fprintf("Setting step plot mod to %i.\n\n", stepPlotMod)
else
    stepPlotMod = expParams.stepPlotMod;
    fprintf("Step plot mod: %i.\n\n", stepPlotMod)
end

% Step record mod
if ~isfield(expParams,"stepRecMod") || ~IsInteger(expParams.stepRecMod) || ... 
        (expParams.stepRecMod <= 0 && expParams.stepRecMod ~= 0 && expParams.stepRecMod ~= -1)
    fprintf("Either no or wrong value for the record mod 'stepRecMod'.\n")
    stepRecMod = -1;  % default step record mod
    fprintf("Setting step record mod to %i.\n\n", stepRecMod)
else
    stepRecMod = expParams.stepRecMod;
    fprintf("Step record mod: %i.\n\n", stepRecMod)
end

% Theta record mod
if ~isfield(expParams,"thetaRecMod") || ~IsInteger(expParams.thetaRecMod) || ... 
        (expParams.thetaRecMod <= 0 && expParams.thetaRecMod ~= 0 && expParams.thetaRecMod ~= -1)
    thetaRecMod = stepRecMod;  % default step record mod
else
    thetaRecMod = expParams.thetaRecMod;
    fprintf("Theta record mod: %i.\n\n", thetaRecMod)
end


%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------


if ~isfield(expParams, "waitForConf") || expParams.waitForConf == true
    fprintf("----------------------------------\n\n")
    fprintf("Press space to start the simulation.\n\n")
    pause
end


fprintf("----------------------------------\n\n")
fprintf("Starting the simulation")
if isfield(expParams,"expTitle") && isstring(expParams.expTitle)
    fprintf(", title: %", expParams.expTitle)
end
fprintf(".\n\n")


% Set auxiliary variablesrecCount
histCoeff = stepDelay;

% Set output variables
% Auxiliary function to determine the count of to be recorded states
function count = getRecCount(module)
    count = 0;
    if module > 1
        count = floor(T / module);
        % To not record last state twice
        if mod(T, module) == 0
            count = count - 1;
        end
    else
        if module ~= -1
            count = T - 1;
        end
    end
end

% Decide the size of xRec
xRecCount = getRecCount(stepRecMod);

% Setup to record initial state
if stepRecMod == 0
    xRec = zeros([N,d,xRecCount + 2]);
    xRec(:,:,1) = x;
    xRecIndex = 2;
    stepRecMod = 1;
else 
    xRec = zeros([N,d,xRecCount + 1]);
    xRecIndex = 1;
end

% Decide the size of thetaRec
thetaRecCount = getRecCount(thetaRecMod);

% Setup to record initial theta
if thetaRecMod == 0
    thetaRec = zeros([N,2,thetaRecCount + 2]);
    thetaRec(:,1,1) = getTheta(getDelayedDists(x,xHist,histCoeff,"None"));
    thetaRec(:,2,1) = getTheta(getDelayedDists(x,xHist,histCoeff,delayType));
    thetaRecIndex = 2;
    thetaRecMod = 1;
else 
    thetaRec = zeros([N,d,thetaRecCount + 1]);
    thetaRecIndex = 1;
end

    

% Saving maximal values of W to normalize 1D graphs
W_max = W(0);

volume = prod(dims);

% Simulate for t=1:T
for t=1:T
    % Distance matrix
    D = getDelayedDists(x,xHist,histCoeff,delayType);

    % Update history of x before changing x
    if delayType ~= "None"
        xHist(:,:,histCoeff) = x;
        histCoeff = histCoeff - 1;
        if histCoeff <= 0
            histCoeff = stepDelay;
        end
    end
    
    % Local density
    theta = getTheta(D);
        
    % Diffusivity
    G_vals = G(theta);
            
    % Calculate the update
    updt = repmat(G_vals,[1,d]).*randn(N,d);
    
    % Make one time step
    x = x + sqrt(dt)*updt;

    switch boundConds
        case "Periodic"
            % Periodic BCs
            x = mod(x,dims);
        case "Ricochet"
            % Ricochet BCs
            x = abs(x);
            x = dims - abs(dims - x);
    end

    % Plot - to make correct 1D plot, we need current theta
    if stepPlotMod > 0 && mod(t-1,stepPlotMod) == 0
        D = getDists(x,x);
        theta = getTheta(D);
        plotSimState(theta)
    end

    % Record simulation state
    if stepRecMod ~= -1 && mod(t,stepRecMod) == 0
        xRec(:,:,xRecIndex) = x;
        xRecIndex = xRecIndex + 1;
    end

    % Record theta
    if thetaRecMod ~= -1 && mod(t,thetaRecMod) == 0
        thetaRec(:,1,thetaRecIndex) = getTheta(getDelayedDists(x,xHist,histCoeff,"None"));
        thetaRec(:,2,thetaRecIndex) = getTheta(getDelayedDists(x,xHist,histCoeff,delayType));
        thetaRecIndex = thetaRecIndex + 1;
    end
end

function plotSimState(theta)
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


% Calculates the distance matrix of the agents with (possible) delay
% Delay is included by taking the last x from x_history using histCoeff
% Takes into account boundary conditions
function D = getDelayedDists(x, xHist, histCoeff, delayType)
    switch delayType
        case "Reaction"
            D = getDists(xHist(:,:,histCoeff),xHist(:,:,histCoeff));
        case "Transmission"
            D = getDists(x,xHist(:,:,histCoeff));
        case "Memory"
            D = getDists(xHist(:,:,histCoeff),x);
        case "None"
            D = getDists(x,x);
        otherwise
            error('Invalid delay type.');
    end
end

% Calculates the distance matrix from the given position matrices
% Takes into account boundary conditions
function D = getDists(x_1,x_2)
    switch boundConds
        case "Periodic"
            % Distances on torus
            D = torusDistances(x_1,x_2,dims);
        case "Ricochet"
            % Normal distances
            D = distances(x_1,x_2);
        otherwise
            error('Invalid boundary conditions.');
    end
end

% Calculates local density from the ((partially) delayed) distance matrix D
function theta = getTheta(D)
    % Subtract 1 from W to exclude the current individual, since its
    sumMatrix = W(D);
    sumMatrix(1:size(D,1)+1:size(D,1)*size(D,2)) = 0;
    theta = (sum(sumMatrix,2)) / (N-1) * volume; 
    % This means the model behaves similarly for constanc volume density of
    % agents: N/volume, and at the same time behaves the same for different
    % amount of agents in the same volume
end

% Record final simulation state
xRec(:,:,end) = x;
thetaRec(:,1,end) = getTheta(getDelayedDists(x,xHist,histCoeff,"None"));
thetaRec(:,2,end) = getTheta(getDelayedDists(x,xHist,histCoeff,delayType));

% Record last history of x - just permute the history
lastIndex = histCoeff + 1;
if lastIndex >= stepDelay
    lastIndex = 1;
end
permutation = [lastIndex:stepDelay, 1:lastIndex-1];
xHist = xHist(:,:,permutation);

% Return random generator settings
rngSetts = rng;

fprintf("----------------------------------\n\n")
fprintf("Simulation finished.\n\n")

if stepPlotMod ~= -2
    % We need to update theta for the last plot
    D = getDists(x,x);
    theta = getTheta(D);

    fprintf("----------------------------------\n\n")
    fprintf("Plotting agregation groups.\n\n")
    plotAgg(x,theta ./ (W_max * volume), dims, boundConds)
end

end