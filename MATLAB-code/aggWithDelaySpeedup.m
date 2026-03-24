function [xRec, thetaRec, thetaOccur, xInitHist, xHist, rngSetts] = aggWithDelaySpeedup(expParams)

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
%       rngSeed (nonnegative integer) - rng seed to replicate experiments.
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
%
%       MODEL:
%       intRad (positive float) - radius of interactions between agents.
%       boundConds (string) - type of the boundary conditions.
%           Must be one of the following strings:
%               "Periodic"
%               "Reflective"
%
%       TIME & DELAY:
%       T (positive integer) - number of time steps.
%       dt (positive float) - time step length.
%       delayType (string) - type of the delay.
%           Must be one of the following strings:
%               "Reaction"
%               "Transmission"
%               "Inner"
%               "None"
%       stepDelay (nonnegative integer) - number of steps used to delay the simulation.
%       xInitHist (float matrix) - initial history of the matrix of positions used in
%           calculation of the first few iterations.
%           x(:,:,i) - position matrix i steps into the past, 
%           where 1 <= i <= stepDelay.
%
%       USER:
%       waitForConf (logical) - if true, then wait for user to start the simulation.
%       stepPlotMod (positive integer or -1 or -2) - simulation plots the 
%           t-th step if the reminder of t divided by stepRecMod is zero. 
%           The plots are shown in movie-like sense. 
%           If stepPlotMod = -1, then only the last step is plotted.
%           If stepPlotMod = -2, then no steps are plotted.
%       markAgents (positive integer vector) - vector containing indeces
%           of the agents to be marked. See 'markColor' for mark color.
%       markColor (nonegative float matrix) - matrix with 3 columns of the
%           same size as 'markAgents'. In each row must be an eligible RGB
%           triplet, which indicates what color to mark the agent with 
%           (specified by 'markAgents' index in the same row).
%       stepRecMod (positive integer or 0 or -1) - simulation records the t-th matrix 
%           of positions if the reminder of t divided by 'stepRecMod' is zero. 
%           Last matrix is always recorded.
%           If stepRecMod = -1, then only the last matrix is recorded.
%       recInitStep (logical) - if true, then the initial positions matrix
%           x0 is recorded to xRec(:,:,1).
%       thetaRecMod (positive integer or 0 or -1) - simulation records the t-th average
%           density vector (theta) and its delayed value if the reminder of t 
%           divided by 'thetaRecMod' is zero.e 
%           Last vector is always recorded.
%           If thetaRecMod = -1, then only the last vector is recorded.
%           If no value is specified, then the value of 'stepRecMod' is
%           used instead.
%       recInitTheta (logical) - if true, then the initial density vector (theta)
%           and its delayed value are recorded to thetaRec(:,:,1).
%           If no value is specified, then the value of 'recInitStep' is
%           used instead.
%       thetaOccurMod (positive integer or 0 or -1) - simulation records the 
%           occurances of an agents with specific numbers of actual neighbours 
%           and delayed neighbours in the t-th step, if the reminder of t 
%           divided by 'thetaOccurMod' is zero.
%           The first step is always counted.
%           The last step is counted only if mod(T,thetaOccurMod) == 0.
%           If thetaOccurMod = -1, then no in-between steps are counted.
%           If no value is specified, then the value of 'thetaRecMod' is
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
%   thetaRec (float matrix) - 3 dimensional matrix of all recorded average 
%       density vectors (theta). Its dimensions are [N,2,count], 
%       where count is the final count of recorded matrices.
%       thetaRec(:,1,i) is the density vector after i-th step of the simulation
%       calculated with NO delay.
%       thetaRec(:,2,i) is the density vector after i-th step of the simulation
%       calculated WITH delay, based on the delay type.
%       thetaRec(:,:,end) is always the density vector at the end of the simulation. 
%       If input thetaRecMod = 0, then thetaRec(:,:,1) is the initial
%       density vector.
%   thetaOccur (float matrix) - 2 dimensional matrix. 
%       Rows represent all the possible amounts of agent within 
%       the interaction radius of one agent. 
%       Columns represent all the possible amounts of agents within 
%       the interaction radius of one agent, shifted by one to right 
%       (i.e. (1,1) represents (0,0)), while taking into account used delay. 
%       Element recTheta(i,j) is then the number of ocurrences of an
%       agent, who had i-1 actual neighbours and j-1 delayed neighbours, in
%       any step which was recorded based on 'thetaRecMod'.
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

fprintf("----------------------------------\n\n")

function result =  IsInteger(x)
    result = isnumeric(x) && all(x == floor(x));
end

% Set or initialize experiment parameters

%---------------------------------RNG---------------------------------------

% Random generator settings for continuation of experiment
if ~isfield(expParams,"rngSetts") || ~isstruct(expParams.rngSetts)
    % We do not have all settings, so we check just the seed
    % Rng seed for replicable experiment
    if ~isfield(expParams,"rngSeed") || ~IsInteger(expParams.rngSeed) || expParams.rngSeed < 0 || ...
            ~isequal(size(expParams.rngSeed),[1,1])
        fprintf("Either no or wrong value for the rng seed 'rngSeed'.\n")
        fprintf("Randomness is uncontrolled.\n\n")
    else
        rng(expParams.rngSeed)
        fprintf("rngSeed: %i.\n\n", expParams.rngSeed)
    end
else
    rng(expParams.rngSetts)
    fprintf("rngSetts accepted.\n\n")
end

%-----------------------SPACE-&-INITIAL-CONDITIONS------------------------

% Initial positions
if ~isfield(expParams,"x0") || ~isfloat(expParams.x0) || isempty(expParams.x0)
    fprintf("Either no or wrong value for the matrix of initial positions 'x0'.\n")
    fprintf("Looking for the input values for 'N' and 'd'.\n")
    fprintf('   |\n')
    if ~isfield(expParams,"N") || ~IsInteger(expParams.N) || expParams.N <= 0 || ...
            ~isequal(size(expParams.N),[1,1])
        fprintf("   Either no or wrong value for the number of agents 'N'.\n")
        N = 400;                % default number of particles
        fprintf("   Setting N = %i.\n", N);
    else
        N = single(expParams.N);
        fprintf("   N = %i.\n", N)
    end
    if ~isfield(expParams,"d") || ~IsInteger(expParams.d) || expParams.d <= 1 || ...
            ~isequal(size(expParams.d),[1,1])
        fprintf("   Either no or wrong value for the dimension 'd'.\n")
        d = 2;                  % default dimension of the space
        fprintf("   Setting d = %i.\n", d);
    else
        d = single(expParams.d);
        fprintf("   d = %i.\n", d)
    end
    fprintf('   |\n')
    fprintf("Initializing random experiment with N = %i, d = %i.\n\n", N, d);
    x = []; % we will set the initial position later as a random matrix
else
    x = expParams.x0;
    N = size(x,1);
    d = size(x,2);
    fprintf("x0 accepted, N = %i, d = %i.\n\n", N, d)
end

% Setting random initial positions
if isempty(x)
    x = rand(N, d);
end

%----------------------------------MODEL------------------------------------

% Interaction radius
kappa = pi^(d/2) / gamma(d / 2 + 1);  % volume of a unit d-ball
if ~isfield(expParams,"intRad") || ~isnumeric(expParams.intRad) || expParams.intRad < 0
    fprintf("Either no or wrong value for the interaction radius 'intRad'.\n")
    kappa_1 = 2;    % "volume" of a unit 1-ball
    kappa_2 = pi;   % "volume" of a unit 2-ball
    intRad_1 = 0.05^2 / kappa_1 * kappa_2;              % interaction radius in 1D
    intRad = (kappa_1 / kappa * intRad_1)^(1/d);        % interaction radius - it will be 0.05 in 2D
    fprintf("Setting intRad = %.3d\n\n", intRad)
else
    intRad = expParams.intRad;
    fprintf("intRad = %.3d\n\n", intRad)
end
W_norm = kappa * intRad^d;      % W_norm to normalize W (done after the summation for efficiency)
intRadSqrd = intRad*intRad;     % Squared raduius for faster calculations

% Boundary conditions
if ~isfield(expParams,"boundConds") || ~isstring(expParams.boundConds) || ...
        ~isequal(size(expParams.boundConds),[1,1])
    fprintf("Either no or wrong value for the boundary conditions 'boundConds'.\n")
    boundConds = "Periodic";    % default boundary conditions
    fprintf("Setting boundary conditions to '%s'.\n\n", boundConds)
else
    boundConds = expParams.boundConds;
    fprintf("boundConds: '%s'.\n\n", boundConds)
end

%------------------------------TIME-&-DELAY---------------------------------

% Step count
if ~isfield(expParams,"T") || ~IsInteger(expParams.T) || expParams.T < 0 || ...
        ~isequal(size(expParams.T),[1,1])
    fprintf("Either no or wrong value for the number of time steps 'T'.\n")
    T = 1000;                    % default number of time steps
    fprintf("Setting T = %i.\n\n", T)
else
    T = expParams.T;
    fprintf("T = %i.\n\n", T)
end

% Time step length
if ~isfield(expParams,"dt") || ~isfloat(expParams.dt) || expParams.dt <= 0 || ...
        ~isequal(size(expParams.dt),[1,1])
    fprintf("Either no or wrong value for the time step length 'dt'.\n")
    dt = 1e-3;                  % default time step length
    fprintf("Setting dt = %.3d\n\n", dt)
else
    dt = expParams.dt;
    fprintf("dt = %.3d\n\n", dt)
end

% Delay type
if ~isfield(expParams,"delayType") || ~isstring(expParams.delayType) || ...
        ~isequal(size(expParams.delayType),[1,1])
    fprintf("Either no or wrong value for the delay type 'delayType'.\n")
    delayType = "Reaction";     % default delay type
    fprintf("Setting delay type to '%s'.\n\n", delayType)
else
    delayType = expParams.delayType;
    fprintf("delayType: '%s'.\n\n", delayType)
end

% Step delay
if ~isfield(expParams,"stepDelay") || ~IsInteger(expParams.stepDelay) || expParams.stepDelay < 0 || ...
        ~isequal(size(expParams.stepDelay),[1,1])
    fprintf("Either no or wrong value for the step delay 'stepDelay'.\n")
    stepDelay = 5;
    fprintf("Setting step delay to %i.\n\n", stepDelay)
else
    stepDelay = expParams.stepDelay;
    fprintf("stepDelay: %i.\n\n", stepDelay)
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
    xHist = genInitHist(x,dt,stepDelay,boundConds,ones(1,d));  % default initial history
    fprintf("Initializing experiment with 'blind motion' initial history.\n\n")
else
    xHist = expParams.xInitHist;
    fprintf("xInitHist accepted.\n\n")
end

% Return initial history (for example if randomly generated)
xInitHist = xHist;

%----------------------------------USER-------------------------------------

% Step plot mod
if ~isfield(expParams,"stepPlotMod") || ~IsInteger(expParams.stepPlotMod) || ...
        ((expParams.stepPlotMod <= 0) && expParams.stepPlotMod ~= -1 && expParams.stepPlotMod ~= -2) || ...
        ~isequal(size(expParams.stepPlotMod),[1,1])
    fprintf("Either no or wrong value for the step plot mod 'stepPlotMod'.\n")
    stepPlotMod = 3;  % default step plot mod
    fprintf("Setting step plot mod to %i.\n\n", stepPlotMod)
else
    stepPlotMod = expParams.stepPlotMod;
    fprintf("stepPlotMod: %i.\n\n", stepPlotMod)
end

if stepPlotMod > 0
    % Marked agents indeces
    if ~isfield(expParams,"markAgents") || ~IsInteger(expParams.markAgents) || any(expParams.markAgents <= 0) || ...
            ~isequal(size(expParams.markAgents),[size(expParams.markAgents,1),1])
        fprintf("Either no or wrong value for the marked agents indeces 'markAgents'.\n")
        markAgents = [];  % default marks
        fprintf("None of the agents will be marked.\n\n")
    else
        markAgents = expParams.markAgents;
        fprintf("markAgents: ")
        disp(markAgents.')
    end
    
    % Marked agents colors
    if ~isfield(expParams,"markColors") || ~isfloat(expParams.markColors) || ...
            any(expParams.markColors < 0) || any(expParams.markColors > 1) || ...
            ~isequal(size(expParams.markColors),[N,3])
        fprintf("Either no or wrong value for the marked agents colors 'markAgents'.\n")
        markColors = repmat([1,0,0],[length(markAgents),1]);
        fprintf("Setting marked agents colors to red.\n\n")
    else
        markColors = expParams.markColors;
        fprintf("markColors accepted.\n\n")
    end
else
    markColors = [];
    fprintf("No agent marking, since plotting is disabled.\n\n")
end

% Step record mod
if ~isfield(expParams,"stepRecMod") || ~IsInteger(expParams.stepRecMod) || ... 
        (expParams.stepRecMod <= 0 && expParams.stepRecMod ~= -1) || ...
        ~isequal(size(expParams.stepRecMod),[1,1])
    fprintf("Either no or wrong value for the step record mod 'stepRecMod'.\n")
    stepRecMod = -1;  % default step record mod
    fprintf("Setting step record mod to %i.\n\n", stepRecMod)
else
    stepRecMod = expParams.stepRecMod;
    fprintf("stepRecMod: %i.\n\n", stepRecMod)
end

% Theta record mod
if ~isfield(expParams,"thetaRecMod") || ~IsInteger(expParams.thetaRecMod) || ... 
        (expParams.thetaRecMod <= 0 && expParams.thetaRecMod ~= -1) || ...
        ~isequal(size(expParams.thetaRecMod),[1,1])
    fprintf("Either no or wrong value for the theta record mod 'thetaRecMod'.\n")
    thetaRecMod = stepRecMod;  % default theta record mod
    fprintf("Setting theta record mod to %i.\n\n", thetaRecMod)
else
    thetaRecMod = expParams.thetaRecMod;
    fprintf("thetaRecMod: %i.\n\n", thetaRecMod)
end

% Theta occurance mod
if ~isfield(expParams,"thetaOccurMod") || ~IsInteger(expParams.thetaOccurMod) || ... 
        (expParams.thetaOccurMod <= 0 && expParams.thetaOccurMod ~= -1) || ...
        ~isequal(size(expParams.thetaOccurMod),[1,1])
    fprintf("Either no or wrong value for the theta occurance mod 'thetaOccurMod'.\n")
    thetaOccurMod = thetaRecMod;  % default theta record mod
    fprintf("Setting theta occurance mod to %i.\n\n", thetaOccurMod)
else
    thetaOccurMod = expParams.thetaOccurMod;
    fprintf("thetaOccurMod: %i.\n\n", thetaOccurMod)
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
if isfield(expParams,"expTitle") && isstring(expParams.expTitle) && isequal(size(expParams.expTitle),[1,1])
    fprintf(", title: %s", expParams.expTitle)
end
fprintf(".\n\n")


% Set auxiliary variables

% Coefficient to access history
histCoeff = stepDelay;

% Colors of agents in plots (to visualize marked agents)
scatterColors = repmat([0.1,0.6,1],[N,1]); % Unmarked color
if ~isempty(markColors)
    scatterColors(markAgents,:) = markColors;
end

% Set output variables

% Auxiliary function to determine the count of to be recorded steps
function count = getRecCount(module)
    % We do not want to record
    if module > 0
        count = ceil(T / module);
        if count == 0
            count = 1; % We always record the last step
        end
    else
        count = 1; % Record just the last step
    end
end

% Decide the size of xRec
xRecCount = getRecCount(stepRecMod);

% Setup to record initial step
if isfield(expParams,"recInitStep") && expParams.recInitStep == true
    xRec = zeros([N,d,xRecCount + 1]);
    xRec(:,:,1) = x;
    xRecIndex = 2;
else 
    xRec = zeros([N,d,xRecCount]);
    xRecIndex = 1;
end

% Decide the size of thetaRec
thetaRecCount = getRecCount(thetaRecMod);

% Setup to record initial theta
if (isfield(expParams,"recInitTheta") && expParams.recInitTheta == true) || ...
        isfield(expParams,"recInitStep") && expParams.recInitStep == true
    thetaRec = zeros([N,2,thetaRecCount + 1]);
    thetaRec(:,1,1) = getTheta(getDistsSqrd(x,x));
    thetaRec(:,2,1) = getTheta(getDelayedDistsSqrd(x,xHist,histCoeff,delayType));
    thetaRecIndex = 2;
else 
    thetaRec = zeros([N,d,thetaRecCount]);
    thetaRecIndex = 1;
end

% Setup to count occurances
thetaOccur = zeros(N+1,N+1);    % 1 <= number of neighbours <= N

% Count initial occurances (so the edge case T = 0 works properly)
if thetaOccurMod ~= -1
    intCountsRealTime = getIntCountsFromDSqrd(getDistsSqrd(x,x));
    intCounts = getIntCountsFromDSqrd(getDelayedDistsSqrd(x,xHist,histCoeff,delayType));
    indexes = [intCountsRealTime(:),intCounts(:)] + 1; % Shift by one to include case where intCount is 0
    thetaOccur = accumarray(indexes, 1, size(thetaOccur));
end

memorizable = (delayType == "Reaction");
% If the delayType is 'Reaction', we can save some calculations
% on recording theta and incrementing thetaOccur
if memorizable
    intCountsHist = zeros(N,stepDelay);
    % Initialize the intCounts history from the initial history of x
    for i = stepDelay:-1:1
        oldX = xHist(:,:,i);
        intCountsHist(:,i) = getIntCountsFromDSqrd(getDistsSqrd(oldX,oldX));
    end
end

% Simulate for t=1:T, this loop calculates the  step t from the step t-1
for t=1:T
    % Interaction counts:
    % Vector, where i-th element is the count of agents 
    % in the interaction radius of the i-th agent (including the i-th)
    if memorizable % We have already calculated it in the past
        intCounts = intCountsHist(:,histCoeff);
        % Update the history
        DSqrd = getDistsSqrd(x,x); % We calculate the current one to be used int the future
        intCountsRealTime = getIntCountsFromDSqrd(DSqrd);
        intCountsHist(:,histCoeff) = intCountsRealTime;
        % The hist coeff will be updated later
    else % We must calculate it anew
        DSqrd = getDelayedDistsSqrd(x,xHist,histCoeff,delayType);
        intCounts = getIntCountsFromDSqrd(DSqrd);
    end

    if delayType == "None"
        intCountsRealTime = intCounts;
    end

    % Update history of x and histCoeff before changing x
    if delayType ~= "None"
        xHist(:,:,histCoeff) = x;
        histCoeff = histCoeff - 1;
        if histCoeff <= 0
            histCoeff = stepDelay;
        end
    end   
    
    % Count theta occurrences from the previous (t-1) step 
    % (to reduce calls of getDelayedDistsSqrd())
    if t > 1 % Initial step is already recorded
        if thetaOccurMod ~= -1 && mod(t-1,thetaOccurMod) == 0
            % If it is memorizable or no delay, intCountsRealTime were already
            % calculated, if not, we need to calculate them now
            if ~memorizable && delayType ~= "None"
                intCountsRealTime = getIntCountsFromDSqrd(getDistsSqrd(x,x));
            end
            
            % Count the occurrences of the theta couples and add them to
            % the whole count
            indexes = [intCountsRealTime(:),intCounts(:)] + 1; % Shift by one to include case where intCount is 0
            prevStepThetaOccur = accumarray(indexes, 1, size(thetaOccur));
            thetaOccur = thetaOccur + prevStepThetaOccur;
        end
    end
    
    % Local density
    theta = getThetaFromIntCounts(intCounts);

    % Record theta from the previous (t-1) step 
    % (to reduce calls of getDelayedDistsSqrd())
    if t > 1 % Initial step is already recorded
        if thetaRecMod ~= -1 && mod(t-1,thetaRecMod) == 0
            % If it is memorizable or no delay, intCountsRealTime were already
            % calculated, if not, we need to calculate them now
            if ~memorizable && delayType ~= "None"
                intCountsRealTime = getIntCountsFromDSqrd(getDistsSqrd(x,x));
            end

            thetaRealTime = getThetaFromIntCounts(intCountsRealTime);
            thetaRec(:,1,thetaRecIndex) = thetaRealTime;
            thetaRec(:,2,thetaRecIndex) = theta;
            thetaRecIndex = thetaRecIndex + 1;
        end
    end
        
    % Diffusivity
    G_vals = exp(-theta);
            
    % Calculate the update
    updt = G_vals .* randn(N,d);
    
    % Make one time step
    x = x + sqrt(dt) * updt;

    % Apply BCs
    switch boundConds
        case "Periodic"
            % Periodic BCs
            x = mod(x,1);
        case "Reflective"
            % Reflective BCs
            x = abs(x);
            x = 1 - abs(1 - x);
    end

    % Plot - to make correct 1D plot, we need current theta
    if stepPlotMod > 0 && mod(t,stepPlotMod) == 0
        DSqrd = getDistsSqrd(x,x);
        theta = getTheta(DSqrd);
        plotSimStep(theta)
    end

    % Record simulation step
    if t < T % Last step is recorded after this loop (since we could get mod(T,stepRecMod) ~= 0, so we record it specially after the loop)
        if stepRecMod ~= -1 && mod(t,stepRecMod) == 0
            xRec(:,:,xRecIndex) = x;
            xRecIndex = xRecIndex + 1;
        end
    end
end

function plotSimStep(theta)
    switch d
        case 1
            scatter(x(:,1),theta,[],scatterColors); 
            axis([0 1 0 1]);
            getframe;
        case 2
            scatter(x(:,1),x(:,2),[],scatterColors); 
            axis([0 1 0 1]);
            getframe;
        case 3
            scatter3(x(:,1),x(:,2),x(:,3),[],scatterColors);
            axis([0 1 0 1 0 1]);
            getframe;
    end
end

% Calculates the distance matrix of the agents with (possible) delay.
% Delay is included by taking the oldest x from x_history using histCoeff.
% Takes into account boundary conditions.
function DSqrd = getDelayedDistsSqrd(x, xHist, histCoeff, delayType)
    switch delayType
        case "Reaction"
            DSqrd = getDistsSqrd(xHist(:,:,histCoeff),xHist(:,:,histCoeff));
        case "Transmission"
            DSqrd = getDistsSqrd(x,xHist(:,:,histCoeff));
        case "Inner"
            DSqrd = getDistsSqrd(xHist(:,:,histCoeff),x);
        case "None"
            DSqrd = getDistsSqrd(x,x);
        otherwise
            error('Invalid delay type.');
    end
end

% Calculates the distance matrix from the given position matrices
% Takes into account boundary conditions
function DSqrd = getDistsSqrd(x_1,x_2)
    switch boundConds
        case "Periodic"
            % Distances on torus
            DSqrd = torusDistancesSqrd(x_1,x_2);
        case "Reflective"
            % Normal distances
            DSqrd = distancesSqrd(x_1,x_2);
        otherwise
            error('Invalid boundary conditions.');
    end
end

% Calculates the local density from the ((partially) delayed) distance matrix DSqrd
function theta = getTheta(DSqrd)
    theta = getThetaFromIntCounts(getIntCountsFromDSqrd(DSqrd));
end

% Calculates the local density of any i-th agent from the given numbers of
% agents that are in the interaction radius of that i-th agent (returns the vector theta)
function theta = getThetaFromIntCounts(intCounts)
    theta = intCounts / W_norm / N;
end

% Calculates the number of agents in the interaction radius of any i-th
% agent (returns a vector of these numbers)
function intCounts = getIntCountsFromDSqrd(DSqrd)
    intCounts = sum((DSqrd < intRadSqrd),2);
end

% Record final simulation step (the final step might not have been possible to record in the loop)
xRec(:,:,end) = x;
thetaRec(:,1,end) = getTheta(getDistsSqrd(x,x));
thetaRec(:,2,end) = getTheta(getDelayedDistsSqrd(x,xHist,histCoeff,delayType));

% Record last history of x - just permute the history
% The coefficient was decreased in the final step, so we need to increase it back
lastIndex = histCoeff + 1; 
if lastIndex >= stepDelay
    lastIndex = 1;
end
permutation = [lastIndex:stepDelay, 1:lastIndex-1];
xHist = xHist(:,:,permutation);

% Return random generator settings
rngSetts = rng;

fprintf("----------------------------------\n\n")
fprintf("Simulation")
if isfield(expParams,"expTitle") && isstring(expParams.expTitle) && isequal(size(expParams.expTitle),[1,1])
    fprintf(" titled %s", expParams.expTitle)
end
fprintf(" finished.\n\n")

if stepPlotMod ~= -2
    % We need to update theta for the last plot
    DSqrd = getDistsSqrd(x,x);
    theta = getTheta(DSqrd);

    fprintf("----------------------------------\n\n")
    fprintf("Plotting agregation groups.\n\n")
    plotAgg(x,theta,ones(1,d),boundConds)
end

end