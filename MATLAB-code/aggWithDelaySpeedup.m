function [xRec, thetaRec, thetaOccur, xHist, rngSetts] = aggWithDelaySpeedup(expParams)

% Runs the discrete simulation of agregation with a constant delay.
%
%--------------------------------------------------------------------------
%
% INPUT:
%   expParams (structure) - struct, which can contain following fields. 
%       If an important field is missing, default value is used.
%
%       POSSIBLE FIELDS:
%
%       RNG:
%       rngSeed (nonnegative integer) - rng seed to replicate experiments.
%       rngSetts (struct) - struct returned by the rng function, containing
%           random generator settings. This struct is also an output of
%           this function, and is created after the simulation is done.
%           This field's main goal is to enable user to continue in an
%           experiment which finished "in the middle", by using the values 
%           it returned as initial conditions and random generator settings.
%
%       SPACE & INITIAL CONDITIONS:
%       domainDims (positive float ROW vector) - dimensions of the domain 
%           of the simulation, i.e., dimensions of the box, in which the 
%           agents move. This box "starts" at vector 0 and "ends" at vector
%           'domainDims'.
%           Alternatively, set 'd'. In that case 'ones(1,d)' will be set as
%           dimensions of the domain
%       d (positive integer) - dimension.
%           Alternative to 'domainDims'.
%       x0 (float matrix) - matrix of initial positions.
%           x(i,:) - position (float vector) of the i-th agent. Positions
%           should be contained within the box defined by 'domainDims'.
%           Alternatively, set the dimension 'd' and the number of
%           agents 'N'. In that case, uniformly random initial positions
%           will be selected.
%       xInitHist (float matrix) - initial history of the matrix of 
%           positions used in the calculation of the first few iterations.
%           x(:,:,i) - position matrix i steps into the past
%           Alternatively, set the dimension 'd', the number of
%           agents 'N' and the step delay 'stepDelay'. In that case, random
%           'blind motion' will be generated as the initial history.
%       N (positive integer) - number of agents.
%           Alternative to 'x0' and 'xInitHist'.
%
%       MODEL:
%       intRad (positive float) - radius of interactions between agents.
%       rd (float) - response decay. The 'a' in 'G(theta) = exp(-a*theta)'.
%       boundConds (string) - type of the boundary conditions.
%           Must be one of the following strings:
%               "NoBoundary"
%               "Periodic"
%               "Reflective"
%
%       TIME & DELAY:
%       dt (positive float) - time step length.
%       stepCount (positive integer) - number of time steps.
%           Alternatively, set the total time 'T'. In that case,
%           'round(T / dt)' will be set as step count.
%       T (positive float) - total time of the simulation.
%           Alternative to 'stepCount'.
%       delayType (string) - type of the delay.
%           Must be one of the following strings:
%               "Reaction"
%               "Transmission"
%               "Inner"
%               "None"
%       stepDelay (nonnegative integer) - number of steps used to delay 
%           the simulation.
%           Alternative to 'xInitHist'.
%           Alternatively, set the (real time) delay 'tau'. In that case,
%           'round(tau / dt)' will be set as step delay.
%       tau (nonnegative float) - delay in the simulation.
%           Alternative to 'stepDelay'.
%
%       USER:
%       waitForConf (logical) - if true, then wait for the user to start 
%           the simulation.
%       stepPlotMod (positive integer or -1) - simulation plots the t-th 
%           step if the reminder of 't' divided by stepPlotMod is 0. 
%           The plots are shown in movie-like sense. 
%           If 'stepPlotMod > 0', then the initial conditions are plotted
%           as well.
%       plotLastAgg (logical) - if true, then the last step of the simulation
%           is plotted together with detected clusters.
%           If no value is specified and 'stepPlotMod > 0', then it will
%           be set 'plotLastAgg = true'.
%       markAgents (positive integer vector) - vector containing indices
%           of agents to be marked. See 'markColors' for the mark color.
%       markColors (nonegative float matrix) - matrix with 3 columns of the
%           same length as 'markAgents'. In each row must be an eligible 
%           RGB triplet, which indicates what color to mark the agent with 
%           (specified by corresponding 'markAgents' index).
%       stepRecMod (positive integer or -1) - simulation records the t-th 
%           matrix of positions if the reminder of 't' divided by 
%           'stepRecMod' is 0. 
%           Last matrix is always recorded.
%           If 'stepRecMod = -1', then only the last matrix is recorded.
%       recInitStep (logical) - if true, then the initial positions matrix
%           'x0' is recorded to 'xRec(:,:,1)'.
%           If no value is specified and 'stepRecMod > 0', then it will
%           be set 'recInitStep = true'.
%       thetaRecMod (positive integer or -1) - simulation records the t-th 
%           density vector (theta) and its delayed value if the reminder 
%           of 't' divided by 'thetaRecMod' is 0.
%           Last vector is always recorded.
%           If 'thetaRecMod = -1', then only the last vector is recorded.
%           If no value is specified, then the value of 'stepRecMod' is
%           used instead.
%       recInitTheta (logical) - if true, then the initial density vector 
%           (theta) and its delayed value are recorded to 
%           'thetaRec(:,:,1)'. 
%           If no value is specified and 'recInitTheta > 0', then it will
%           be set 'recInitTheta = true'.
%       thetaOccurMod (positive integer or -1) - simulation records the 
%           occurances of agents with specific numbers of actual neighbours 
%           and delayed neighbours in the t-th step, if the reminder of 't' 
%           divided by 'thetaOccurMod' is 0.
%           The last step is always counted.
%           If 'thetaOccurMod = -1', then only the last step is counted.
%           If no value is specified, then the value of 'thetaRecMod' is
%           used instead.
%       countInitThetaOccur (logical) - if true, then the initial theta
%           occurences (as in 'thetaOccurMod') are counted as well.
%           If no value is specified and 'thetaOccurMod > 0', then it will
%           be set 'countInitThetaOccur = true'.
%       expTitle (string) - title to be printed before the experiment 
%           begins, and after the experiment ends.
%       videoPath (string) - if provided, script records the plotted
%           positions as video onto this given path. This path should
%           include the filename, but NOT the file extension. It will be
%           saved as '.avi'.
%
% OUTPUT:
%   xRec (float matrix) - 3 dimensional matrix of all recorded matrices of
%       positions. Its dimensions are '[N,d,count]', where 'count' is the 
%       final number of recorded matrices. 
%       'xRec(:,:,end)' is always the matrix of positions at the end of the 
%       simulation. 
%       If 'recInitStep = true', then 'xRec(:,:,1)' is the matrix of 
%       initial positions.
%   thetaRec (float matrix) - 3 dimensional matrix of all recorded density 
%       vectors (theta). Its dimensions are '[N,2,count]', where 'count' is 
%       the final number of recorded matrices.
%       'thetaRec(:,1,i)' is the actual (without delay) density vector 
%       after the i-th step of the simulation.
%       'thetaRec(:,2,i)' is the delayed density vector after the i-th step
%       of the simulation.
%       'thetaRec(:,:,end)' are always (non/delayed) density vectors at the
%       end of the simulation. 
%       If 'recInitTheta = true', then 'thetaRec(:,:,1)' are the initial
%       density vectors.
%   thetaOccur (float matrix) - 2 dimensional matrix. 
%       Row indices represent all the possible amounts of agents within 
%       the interaction radius of one agent, shifted by one to right.
%       Column indices represent all the possible amounts of agents within 
%       the interaction radius of one agent, shifted by one to right, while
%       taking into account used delay. 
%       The shift means that, e.g., '(1,1)' represents '(0,0)'.
%       Element 'thetaOccur(i,j)' is then the number of ocurrences of any
%       agent, who had 'i-1' actual neighbours and 'j-1' delayed 
%       neighbours, in any step which was recorded based on 'thetaRecMod'.
%   xHist (float matrix) - 3 dimensional matrix of the last history of
%       positions, which still affects the following simulation steps.
%       Its dimensions are '[N,d,stepDelay]'. It can be directly plugged in
%       as the input to this function to continue in the experiment.
%   rngSetts (struct) - struct returned by the rng function, containing
%       random generator settings right after the simulation have finished.
%       It can be directly plugged in as the input to this function to 
%       continue in the experiment.

arguments
    expParams (1,1) struct = struct()
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

fprintf("----------------------------------\n\n")

fprintf("Initializing the experiment: Agregation with delay.\n\n")

fprintf("----------------------------------\n\n")

function result = IsAllModInteger(x)
    result = IsAllInteger(x) && all((x > 0) | (x == -1));
end

function result = IsAllNonnegInteger(x)
    result = IsAllInteger(x) && all(x >= 0);
end

function result = IsAllPosInteger(x)
    result = IsAllInteger(x) && all(x > 0);
end

function result = IsAllInteger(x)
    result = IsAllReal(x) && all(x == floor(x));
end

function result = IsAllNonnegReal(x)
    result = IsAllReal(x) && all(x >= 0);
end

function result = IsAllPosReal(x)
    result = IsAllReal(x) && all(x > 0);
end

function result = IsAllReal(x)
    result = isnumeric(x) && isreal(x) && all(isfinite(x));
end

% Set or initialize experiment parameters ---------------------------------

%---------------------------------RNG--------------------------------------

% Random generator settings for continuation of experiment
if ~isfield(expParams,"rngSetts") || ~isscalar(expParams.rngSetts) || ...
        ~isstruct(expParams.rngSetts)
    % We do not have all settings, so we check just the seed
    % Rng seed for replicable experiment
    if ~isfield(expParams,"rngSeed") || ~isscalar(expParams.rngSeed) || ...
            ~IsAllNonnegInteger(expParams.rngSeed)
        fprintf("rngSeed - Either no or wrong value for the rng seed.\n")
        fprintf("Randomness is uncontrolled.\n\n")
    else
        rng(expParams.rngSeed)
        fprintf("rngSeed = %i.\n\n", expParams.rngSeed)
    end
else
    rng(expParams.rngSetts)
    fprintf("rngSetts accepted.\n\n")
end

%-----------------------SPACE-&-INITIAL-CONDITIONS-------------------------

% Dimension & dimensions
if ~isfield(expParams, "domainDims") || ~isvector(expParams.domainDims) || ...
        ~IsAllPosReal(expParams.domainDims)
    fprintf("domainDims - Either no or wrong value for the vector of dimensions.\n")
    fprintf("Looking for the input value for 'd'.\n")
    if ~isfield(expParams,"d") || ~isscalar(expParams.d) || ~IsAllPosInteger(expParams.d)
        fprintf("   d - Either no or wrong value for the dimension.\n")
        d = 2;                  % default dimension of the space
        fprintf("   Setting d = %i.\n", d);
    else
        d = double(expParams.d);
        fprintf("   d = %i.\n", d)
    end
    domainDims = ones(1,d);     % default dimensions
    fprintf("Setting domainDims:")
    disp(domainDims)
else 
    domainDims = double(expParams.domainDims(:).');
    d = size(domainDims, 2);
    fprintf("domainDims:")
    disp(domainDims)
end
volume = prod(domainDims);

% Initial positions & number of agents
if ~isfield(expParams,"x0") || ~IsAllReal(expParams.x0) || isempty(expParams.x0) || ...
        ~ismatrix(expParams.x0) || size(expParams.x0, 2) ~= d
    fprintf("x0 - Either no or wrong value for the matrix of initial positions.\n")
    fprintf("Looking for the input value for 'N'.\n")
    fprintf('   |\n')
    if ~isfield(expParams,"N") || ~isscalar(expParams.N) || ~IsAllPosInteger(expParams.N)
        fprintf("   N - Either no or wrong value for the number of agents.\n")
        N = 400;                % default number of particles
        fprintf("   Setting N = %i.\n", N);
    else
        N = double(expParams.N);
        fprintf("   N = %i.\n", N)
    end
    fprintf('   |\n')
    fprintf("Initializing random experiment with N = %i, d = %i.\n\n", N, d);
    x = []; % we will set the initial position later as a random matrix
else
    x = double(expParams.x0);
    N = size(x,1);
    fprintf("x0 accepted, N = %i, d = %i.\n\n", N, d)
end

% Initial history of positions
if ~isempty(x)
    if ~isfield(expParams, "xInitHist") || ~IsAllReal(expParams.xInitHist) || ...
            size(expParams.xInitHist, 1) ~= N || size(expParams.xInitHist, 2) ~= d || ...
            (isempty(expParams.xInitHist) && ndims(expParams.xInitHist) ~= 3) || ... % check that it is empty only when it is 3d = no delay defined => can accept
            ndims(expParams.xInitHist) > 3
        fprintf("xInitHist - Either no or wrong value for the matrix of initial history of positions.\n")
        xHist = [];  % default initial history will be made after stepDelay
        fprintf("Initializing experiment with 'blind motion' initial history.\n\n")
    else
        xHist = double(expParams.xInitHist);
        fprintf("xInitHist accepted.\n\n")
    end
else % It does not make sense to have a history for unknown initial conditions
    xHist = [];
    fprintf("xInitHist ignored, since 'x0' was not specified")
end

% Setting random initial positions
if isempty(x)
    x = rand(N, d) .* (domainDims);
end

%----------------------------------MODEL-----------------------------------

% Interaction radius
if ~isfield(expParams, "intRad") || ~isscalar(expParams.intRad) || ...
        ~IsAllPosReal(expParams.intRad)
    fprintf("intRad - Either no or wrong value for the interaction radius.\n")
    intRad = getIntRad(d);
    fprintf("Setting intRad = %.3e\n\n", intRad)
else
    intRad = double(expParams.intRad);
    fprintf("intRad = %.3e\n\n", intRad)
end
kappa = pi^(d/2) / gamma(d / 2 + 1);    % volume of a unit d-ball
W_norm = kappa * intRad^d;              % W_norm to normalize W (done after the summation for efficiency)
intRadSqrd = intRad*intRad;             % Squared raduius for faster calculations

% Response decay
if ~isfield(expParams, "rd") || ~isscalar(expParams.rd) || ~IsAllReal(expParams.rd)
    fprintf("rd - Either no or wrong value for the response decay.\n")
    rd = 1;
    fprintf("Setting rd = %.3e\n\n", rd)
else
    rd = double(expParams.rd);
    fprintf("rd = %.3e\n\n", rd)
end

definedBoundConds = ["NoBoundary", "Periodic", "Reflective"];

% Boundary conditions
if ~isfield(expParams, "boundConds") || ~isscalar(expParams.boundConds) || ...
        ~isstring(expParams.boundConds) || ~ismember(expParams.boundConds, definedBoundConds)     
    fprintf("boundConds - Either no or wrong value for the boundary conditions.\n")
    boundConds = "Periodic";    % default boundary conditions
    fprintf("Setting boundConds = '%s'.\n\n", boundConds)
else
    boundConds = expParams.boundConds;
    fprintf("boundConds: '%s'.\n\n", boundConds)
end

%------------------------------TIME-&-DELAY--------------------------------

% Time step length
if ~isfield(expParams, "dt") || ~isscalar(expParams.dt) || ~IsAllPosReal(expParams.dt)
    fprintf("dt - Either no or wrong value for the time step length.\n")
    dt = 1e-3;                  % default time step length
    fprintf("Setting dt = %.3e\n\n", dt)
else
    dt = double(expParams.dt);
    fprintf("dt = %.3e\n\n", dt)
end

% Step count
if ~isfield(expParams, "stepCount") || ~isscalar(expParams.stepCount) || ...
        ~IsAllPosInteger(expParams.stepCount)
    if ~isfield(expParams,"T") || ~isscalar(expParams.T) || ~IsAllPosReal(expParams.T)
        fprintf("stepCount (T) - Either no or wrong value for both the number of time steps" + ...
            " and the total simulation time.\n")
        stepCount = 1000;                    % default number of time steps
        fprintf("Setting stepCount = %i (T = %.3e).\n\n", stepCount, stepCount * dt)
    else
        stepCount = round(double(expParams.T) / dt);
        if stepCount == 0
            stepCount = 1;
        end
        T = stepCount * dt;
        fprintf("T = %.3e (stepCount = %i).\n\n", T, stepCount)
    end
else
    stepCount = double(expParams.stepCount);
    fprintf("stepCount = %i (T = %.3e).\n\n", stepCount, stepCount * dt)
end

% Step delay
if isempty(xHist) && ndims(xHist) ~= 3 % = wrong or no xInitHist input

    if ~isfield(expParams, "stepDelay") || ~isscalar(expParams.stepDelay) || ...
            ~IsAllNonnegInteger(expParams.stepDelay)

        if ~isfield(expParams, "tau") || ~isscalar(expParams.tau) || ...
                ~IsAllNonnegReal(expParams.tau)
            fprintf("stepDelay (tau) - Either no or wrong value for both the step delay" + ...
                " and the real time delay.\n")
            stepDelay = getStepDelay(d);                    % default stepDelay
            fprintf("Setting stepDelay = %i (tau = %.3e).\n\n", stepDelay, stepDelay * dt)
        else
            stepDelay = round(double(expParams.tau) / dt);
            tau = stepDelay * dt;
            fprintf("tau = %.3e (stepDelay = %i).\n\n", tau, stepDelay)
        end

    else
        stepDelay = double(expParams.stepDelay);
        fprintf("stepDelay: %i (tau = %.3e).\n\n", stepDelay, stepDelay * dt)
    end
% We got correct xInitHist
else
    stepDelay = size(xHist, 3);
    fprintf("Inferring step delay from 'xInitHist': stepDelay = %i (tau = %.3e).\n\n", ...
        stepDelay, stepDelay * dt)
end

% Make default xInihHist (if wrong or no input)
if isempty(xHist) && ndims(xHist) ~= 3
    xHist = genInitHist(x, dt, stepDelay, boundConds, domainDims);
end

definedDelayTypes = ["Reaction", "Transmission", "Inner", "None"];

% Delay type
if stepDelay ~= 0
    if ~isfield(expParams, "delayType") || ~isscalar(expParams.delayType) || ...
            ~isstring(expParams.delayType) || ~ismember(expParams.delayType, definedDelayTypes)
        fprintf("delayType - Either no or wrong value for the delay type.\n")
        delayType = "Reaction";     % default delay type
        fprintf("Setting delayType = '%s'.\n\n", delayType)
    else
        delayType = expParams.delayType;
        fprintf("delayType = '%s'.\n\n", delayType)
    end
else
    delayType = "None";
    fprintf("Zero step delay detected. Forcing delayType = %s.\n\n", delayType)
end

if stepDelay > 0 && delayType == "None"
    fprintf("There is nonzero 'stepDelay' but 'delayType' is 'None' - simulation will have no delay.\n\n")
end

%----------------------------------USER------------------------------------

% Step plot mod
if ~isfield(expParams,"stepPlotMod") || ~isscalar(expParams.stepPlotMod) || ...
        ~IsAllModInteger(expParams.stepPlotMod)
    fprintf("stepPlotMod - Either no or wrong value for the step plot mod.\n")
    stepPlotMod = getStepPlotMod(d);
    fprintf("Setting stepPlotMod = %i.\n\n", stepPlotMod)
else
    stepPlotMod = expParams.stepPlotMod;
    fprintf("stepPlotMod: %i.\n\n", stepPlotMod)
end

% Plot aggregation
if ~isfield(expParams, "plotLastAgg") || ~isscalar(expParams.plotLastAgg) || ...
        ~islogical(expParams.plotLastAgg)
    if stepPlotMod > 0
        plotLastAgg = true;
    else
        plotLastAgg = false;
    end
else
    plotLastAgg = expParams.plotLastAgg;
end

if stepPlotMod > 0
    % Marked agents indeces
    if ~isfield(expParams,"markAgents") || ~isvector(expParams.markAgents) || ...
            ~IsAllPosInteger(expParams.markAgents) || ~all(expParams.markAgents <= N)
        fprintf("markAgents - Either no or wrong value for the marked agents indeces.\n")
        markAgents = [];  % default marks
        fprintf("None of the agents will be marked.\n\n")
    else
        markAgents = expParams.markAgents;
        if isrow(markAgents)
            markAgents = markAgents.';
        end
        fprintf("markAgents: ")
        disp(markAgents.')
    end
    
    % Marked agents colors
    if ~isfield(expParams,"markColors") || ~IsAllPosReal(expParams.markColors) || ...
            ~all(expParams.markColors <= 1) || ...
            ~isequal(size(expParams.markColors), [length(markAgents), 3])
        fprintf("markAgents - Either no or wrong value for the marked agents colors.\n")
        markColors = repmat([1,0,0], [length(markAgents),1]);
        fprintf("Setting all marked agents colors to red.\n\n")
    else
        markColors = expParams.markColors;
        fprintf("markColors accepted.\n\n")
    end
else
    markAgents = [];
    markColors = [];
    fprintf("No agent marking, since plotting is disabled.\n\n")
end

% Step record mod
if ~isfield(expParams,"stepRecMod") || ~isscalar(expParams.stepRecMod) || ...
        ~IsAllModInteger(expParams.stepRecMod) 
    fprintf("stepRecMod - Either no or wrong value for the step record mod.\n")
    stepRecMod = -1;  % default step record mod
    fprintf("Setting stepRecMod = %i.\n\n", stepRecMod)
else
    stepRecMod = expParams.stepRecMod;
    fprintf("stepRecMod: %i.\n\n", stepRecMod)
end

% Record initial step
if ~isfield(expParams, "recInitStep") || ~isscalar(expParams.recInitStep) || ...
        ~islogical(expParams.recInitStep)
    if stepRecMod > 0
        recInitStep = true;
    else
        recInitStep = false;
    end
else
    recInitStep = expParams.recInitStep;
end

% Theta record mod
if ~isfield(expParams,"thetaRecMod") || ~isscalar(expParams.thetaRecMod) || ...
        ~IsAllModInteger(expParams.thetaRecMod)        
    fprintf("thetaRecMod - Either no or wrong value for the theta record mod.\n")
    thetaRecMod = stepRecMod;  % default theta record mod
    fprintf("Setting thetaRecMod = %i.\n\n", thetaRecMod)
else
    thetaRecMod = expParams.thetaRecMod;
    fprintf("thetaRecMod: %i.\n\n", thetaRecMod)
end

% Record initial theta
if ~isfield(expParams, "recInitTheta") || ~isscalar(expParams.recInitTheta) || ...
        ~islogical(expParams.recInitTheta)
    if thetaRecMod > 0
        recInitTheta = true;
    else
        recInitTheta = false;
    end
else
    recInitTheta = expParams.recInitTheta;
end

% Theta occurance mod
if ~isfield(expParams, "thetaOccurMod") || ~isscalar(expParams.thetaOccurMod) || ...
        ~IsAllModInteger(expParams.thetaOccurMod)
    fprintf("thetaOccurMod - Either no or wrong value for the theta occurance mod.\n")
    thetaOccurMod = thetaRecMod;  % default theta record mod
    fprintf("Setting thetaOccurMod = %i.\n\n", thetaOccurMod)
else
    thetaOccurMod = expParams.thetaOccurMod;
    fprintf("thetaOccurMod: %i.\n\n", thetaOccurMod)
end

% Record initial theta
if ~isfield(expParams, "countInitThetaOccur") || ~isscalar(expParams.countInitThetaOccur) || ...
        ~islogical(expParams.countInitThetaOccur)
    if thetaOccurMod > 0
        countInitThetaOccur = true;
    else
        countInitThetaOccur = false;
    end
else
    countInitThetaOccur = expParams.countInitThetaOccur;
end

% Video path
if ~isfield(expParams, "videoPath") || ~isscalar(expParams.videoPath) || ...
        ~isstring(expParams.videoPath)
    videoPath = [];
else
    videoPath = expParams.videoPath;
end

% Wait for conformation
if ~isfield(expParams, "waitForConf") || ~isscalar(expParams.waitForConf) || ...
        ~islogical(expParams.waitForConf)
    waitForConf = true;
else
    waitForConf = expParams.waitForConf;
end


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


% Set auxiliary variables -------------------------------------------------


fprintf("----------------------------------\n\n")

fprintf("Preparing inner variables.\n\n")


% Coefficient to access history
histCoeff = stepDelay;

% Colors of agents in plots (to visualize marked agents)
scatterColors = repmat([0.1,0.6,1],[N,1]); % Unmarked color
if ~isempty(markColors)
    scatterColors(markAgents,:) = markColors;
end

% Setup VideoWriter if the simulation is recorded
recordVideo = false;
if ~isempty(videoPath)
    if stepPlotMod > 0
        fprintf("Making the video of the plot at %s.\n\n", expParams.videoPath)
        recordVideo = true;
        vidWrit = VideoWriter(sprintf("%s.avi",expParams.videoPath));
        frameTime = stepPlotMod * dt;
        vidWrit.FrameRate = 1 / frameTime;
        open(vidWrit)
    else
        fprintf("NO video recording since there is NO plotting.\n\n")
    end
else
    fprintf("NO video recording.\n\n")
end

% Make new figure if simulation is plotted and plot initial positions
if stepPlotMod > 0
    figure
    hold off
    DSqrd = getDistsSqrd(x,x);
    theta = getTheta(DSqrd);
    plotSimStep(theta)
end

maxPlot = 0;


% Set output variables ----------------------------------------------------


fprintf("----------------------------------\n\n")

fprintf("Preparing output variables.\n\n")


% Auxiliary function to determine the count of to be recorded steps
function count = getRecCount(module)
    % We want to record
    if module > 0
        count = ceil(stepCount / module);
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
if recInitStep
    fprintf("Recording initial step.\n\n")
    xRec = zeros([N,d,xRecCount + 1]);
    xRec(:,:,1) = x;
    xRecIndex = 2;
else 
    fprintf("NOT recording initial step.\n\n")
    xRec = zeros([N,d,xRecCount]);
    xRecIndex = 1;
end

% Decide the size of thetaRec
thetaRecCount = getRecCount(thetaRecMod);

% Setup to record initial theta
if recInitTheta
    fprintf("Recording initial theta.\n\n")
    thetaRec = zeros([N, 2, thetaRecCount + 1]);
    thetaRec(:,1,1) = getTheta(getDistsSqrd(x, x));
    thetaRec(:,2,1) = getTheta(getDelayedDistsSqrd(x, xHist, histCoeff, delayType));
    thetaRecIndex = 2;
else
    fprintf("NOT recording initial theta.\n\n")
    thetaRec = zeros([N, 2, thetaRecCount]);
    thetaRecIndex = 1;
end

% Setup to count occurances
thetaOccur = zeros(N+1,N+1);    % 0 <= number of neighbours <= N (N + 1 due to possible zero neighbours in Transmission type delay)


% Count initial occurances (so the edge case stepCount = 0 works properly)
if countInitThetaOccur
    fprintf("Recording initial theta occurances.\n\n")
    intCountsRealTime = getIntCountsFromDSqrd(getDistsSqrd(x, x));
    intCountsDelayed = getIntCountsFromDSqrd(getDelayedDistsSqrd(x, xHist, histCoeff, delayType));
    indexes = [intCountsRealTime(:), intCountsDelayed(:)] + 1; % Shift by one to include case where intCount is 0 (this can happen in Transmission type delay)
    thetaOccur = accumarray(indexes, 1, size(thetaOccur));
else
    fprintf("NOT recording initial theta occurances.\n\n")
end


% Decide on & prepare memorization ----------------------------------------


memorizable = (delayType == "Reaction");
% If the delayType is 'Reaction', we can save some calculations
% on recording theta and incrementing thetaOccur
if memorizable
    fprintf("----------------------------------\n\n")
    fprintf("Using memorization for speedup.\n\n")
    intCountsHist = zeros(N, stepDelay);
    % Initialize the intCounts history from the initial history of x
    for i = stepDelay:-1:1
        oldX = xHist(:,:,i);
        intCountsHist(:,i) = getIntCountsFromDSqrd(getDistsSqrd(oldX, oldX));
    end
end


%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------


if waitForConf
    fprintf("----------------------------------\n\n")
    fprintf("Press any key to start the simulation.\n\n")
    pause
end


fprintf("----------------------------------\n\n")
fprintf("Starting the simulation")
printTitle()
fprintf(".\n\n")


% Simulate for t = 1:stepCount --------------------------------------------
% This loop calculates the step t from the step t-1 -----------------------

for t = 1:stepCount
    % Interaction counts:
    % Vector, where i-th element is the count of agents 
    % in the interaction radius of the i-th agent (including the i-th)
    if memorizable % We have already calculated it in the past
        intCountsDelayed = intCountsHist(:,histCoeff);
        % Update the history
        DSqrd = getDistsSqrd(x, x); % We calculate the current one to be used int the future
        intCountsRealTime = getIntCountsFromDSqrd(DSqrd);
        intCountsHist(:,histCoeff) = intCountsRealTime;
        % The hist coeff will be updated later
    else % We must calculate it anew
        DSqrd = getDelayedDistsSqrd(x, xHist, histCoeff, delayType);
        intCountsDelayed = getIntCountsFromDSqrd(DSqrd);
    end

    if delayType == "None"
        intCountsRealTime = intCountsDelayed;
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
                intCountsRealTime = getIntCountsFromDSqrd(getDistsSqrd(x, x));
            end
            
            % Count the occurrences of the theta couples and add them to
            % the whole count
            indexes = [intCountsRealTime(:), intCountsDelayed(:)] + 1; % Shift by one to include case where intCount is 0
            prevStepThetaOccur = accumarray(indexes, 1, size(thetaOccur));
            thetaOccur = thetaOccur + prevStepThetaOccur;
        end
    end
    
    % Local density
    thetaDelayed = getThetaFromIntCounts(intCountsDelayed);

    % Record theta from the previous (t-1) step
    % (to reduce calls of getDelayedDistsSqrd())
    if t > 1 % Initial step is already recorded
        if thetaRecMod ~= -1 && mod(t-1,thetaRecMod) == 0
            % If it is memorizable or no delay, intCountsRealTime were already
            % calculated, if not, we need to calculate them now
            if ~memorizable && delayType ~= "None"
                intCountsRealTime = getIntCountsFromDSqrd(getDistsSqrd(x, x));
            end

            thetaRealTime = getThetaFromIntCounts(intCountsRealTime);
            thetaRec(:,1,thetaRecIndex) = thetaRealTime;
            thetaRec(:,2,thetaRecIndex) = thetaDelayed;
            thetaRecIndex = thetaRecIndex + 1;
        end
    end
        
    % Diffusivity
    G_vals = exp(-rd * thetaDelayed);
            
    % Calculate the update
    updt = G_vals .* randn(N,d);
    
    % Make one time step
    x = x + sqrt(dt) * updt;

    % Apply BCs
    x = applyBCs(x, boundConds, domainDims);

    % Plot - to make correct 1D plot, we need current theta
    if stepPlotMod > 0 && mod(t, stepPlotMod) == 0
        DSqrd = getDistsSqrd(x,x);
        theta = getTheta(DSqrd);
        plotSimStep(theta)
    end

    % Record simulation step
    if t < stepCount % Last step is recorded after this loop (since we could get mod(stepCount,stepRecMod) ~= 0, so we record it specially after the loop)
        if stepRecMod ~= -1 && mod(t,stepRecMod) == 0
            xRec(:,:,xRecIndex) = x;
            xRecIndex = xRecIndex + 1;
        end
    end
end

function plotSimStep(theta)
    f = [];
    switch d
        case 1
            maxPlot = max(maxPlot, max(theta) * getMaxPlotMultip());
            scatter(x(:,1), theta, getPointSize()^2, scatterColors, '.');
            axis([0 domainDims(1) 0 1]);
            f = getframe;
        case 2
            scatter(x(:,1), x(:,2), getPointSize()^2, scatterColors, '.'); 
            axis([0 domainDims(1) 0 domainDims(2)]);
            f = getframe;
        case 3
            scatter3(x(:,1), x(:,2), x(:,3), getPointSize()^2, scatterColors, '.');
            axis([0 domainDims(1) 0 domainDims(2) 0 domainDims(3)]);
            f = getframe;
    end

    if recordVideo && ~isempty(f)
        writeVideo(vidWrit,f)
    end
end

% Calculates the distance matrix of the agents with (possible) delay.
% Delay is included by taking the oldest x from x_history using histCoeff.
% Takes into account boundary conditions.
function DSqrd = getDelayedDistsSqrd(x, xHist, histCoeff, delayType)
    switch delayType
        case "Reaction"
            DSqrd = getDistsSqrd(xHist(:,:,histCoeff), xHist(:,:,histCoeff));
        case "Transmission"
            DSqrd = getDistsSqrd(x, xHist(:,:,histCoeff));
        case "Inner"
            DSqrd = getDistsSqrd(xHist(:,:,histCoeff), x);
        case "None"
            DSqrd = getDistsSqrd(x, x);
        otherwise
            error('Invalid delay type.');
    end
end

% Calculates the distance matrix from the given position matrices
% Takes into account boundary conditions
function DSqrd = getDistsSqrd(x_1, x_2)
    switch boundConds
        case "NoBoundary"
            % No BCs
            DSqrd = distancesSqrd(x_1, x_2);
        case "Periodic"
            % Distances on torus
            DSqrd = torusDistancesSqrd(x_1, x_2, domainDims);
        case "Reflective"
            % Normal distances
            DSqrd = distancesSqrd(x_1, x_2);
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
    theta = intCounts / W_norm / N * volume;
end

% Calculates the number of agents in the interaction radius of any i-th
% agent (returns a vector of these numbers)
function intCounts = getIntCountsFromDSqrd(DSqrd)
    intCounts = sum((DSqrd < intRadSqrd), 2);
end


% Finish

% Close video writer
if recordVideo
    close(vidWrit)
end

% Record final simulation step (the final step might not have been possible to record in the loop)
% xRec
xRec(:,:,end) = x;
% thetaRec
thetaRec(:,1,end) = getTheta(getDistsSqrd(x, x));
thetaRec(:,2,end) = getTheta(getDelayedDistsSqrd(x, xHist, histCoeff, delayType));
% thetaOccur
intCountsRealTime = getIntCountsFromDSqrd(getDistsSqrd(x, x));
intCountsDelayed = getIntCountsFromDSqrd(getDelayedDistsSqrd(x, xHist, histCoeff, delayType));
indexes = [intCountsRealTime(:), intCountsDelayed(:)] + 1; % Shift by one to include case where intCount is 0 (this can happen in Transmission type delay)
thetaOccur = thetaOccur + accumarray(indexes, 1, size(thetaOccur));

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
printTitle()
fprintf(" finished.\n\n")

if plotLastAgg
    % We need to update theta for the last plot
    DSqrd = getDistsSqrd(x,x);
    theta = getTheta(DSqrd);

    fprintf("----------------------------------\n\n")
    fprintf("Plotting agregation groups.\n\n")
    plotAgg(x, theta * W_norm / volume, boundConds, domainDims)
end

function printTitle()
    if isfield(expParams,"expTitle") && isstring(expParams.expTitle) && isscalar(expParams.expTitle)
    fprintf(" titled %s", expParams.expTitle)
    end
end
        

end