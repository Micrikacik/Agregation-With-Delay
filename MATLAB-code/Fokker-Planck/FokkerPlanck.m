function [rhoRec, rhoHist] = FokkerPlanck(expParams)

% Runs the discretized continuous simulation of agregation with a constant delay.
%
%---------------------------------------------------------------------------
%
% INPUT:
%   expParams (structure) - struct, which can contain following fields, if an important
%       field is missing, default value is used. If no input is required,
%       set to '{}'.
%
%       POSSIBLE FIELDS:
%
%       RNG:
%       rngSeed (nonnegative integer) - rng seed to replicate experiments.
%       rngSetts (struct) - struct returned by the rng function, containing
%           random generator settings.
%
%       SPACE & INITIAL CONDITIONS:
%       rho0 (positive float COLUMN vector) - vector of initial discretized density. 
%           Alternatively, set number of grid points 'gridPointCount'.
%       gridPointCount (integer scalar) - number of points at which we
%           disretize the differential equation using uniform step.
%       L (positive float) - length of the simulation interval [0, L].
%       maxInitFlucts (positive float scalar) - if either 'rho0'
%           'rhoInitHist' is not provided, this value is used to generate
%           them with random fluctuations, which are at most 'maxInitFlucts'.
%
%       MODEL:
%       intRad (positive float) - radius of interactions between agents.
%       boundConds (string) - type of the boundary conditions.
%           Must be one of the following strings:
%               "Periodic"
%               "Reflective"
%       rd (positive float scalar) - it is a multiplier 'a' in
%           the response function G(-a*theta)
%
%       TIME & DELAY:
%       dt (positive float) - time step length.
%       stepCount (positive integer) - number of time steps.
%       T (positive float) - total time of the simulation, alternative to 
%           'stepCount', which will be 'stepCount = round(T / dt)', if not provided.
%       stepDelay (nonnegative integer) - number of steps used to delay the simulation.
%       tau (nonnegative float) - delay in the simulation, alternative to 
%           'stepDelay', which will be 'stepDelay = round(tau / dt)', if not provided.
%       rhoInitHist (float matrix) - initial history of the matrix of positions used in
%           calculation of the first few iterations.
%           rhoInitHist(:,i) - density vector i steps into the past, 
%           where 1 <= i <= stepDelay.
%
%       USER:
%       waitForConf (logical scalar) - if true, then wait for user to start the simulation.
%       stepPlotMod (positive integer or -1 or -2) - simulation plots the 
%           t-th step if the reminder of t divided by stepRecMod is zero. 
%           The plots are shown in movie-like sense. 
%           If stepPlotMod = -1, then only the last step is plotted.
%           If stepPlotMod = -2, then no steps are plotted.
%       stepRecMod (positive integer or -1) - simulation records the t-th 
%           density vector if the reminder of t divided by 'stepRecMod' is zero. 
%           Last density vector is always recorded.
%           If stepRecMod = -1, then only the last vector is recorded.
%       recInitStep (logical scalar) - if true, then the initial density
%           vector rho0 is recorded to rhoRec(:,1).
%       expTitle (string) - title to be printed before the experiment begins
%       recordVideoPath (string) - if provided, script records the plotted
%           positions as video onto this given path. This path should
%           include the filename, but NOT the file extension. It will be
%           saved as '.avi'.
%
% OUTPUT:
%   rhoRec (float matrix) - 2 dimensional matrix of all recorded density 
%       vectors. Its dimensions are [gridPointCount, count], where count 
%       is the final count of recorded matrices. 
%       rhoRec(:,end) is always the density vector at the end of the simulation. 
%       If input stepRecMod = 0, then rhoRec(:,1) is the vector of initial
%       density.
%   rhoHist (float matrix) - 2 dimensional matrix of the last history of
%       densities, which still affect the following simulation steps.
%       Its dimensions are [gridPointCount, stepDelay]. 
%       It can be directly plugged in as an input to this function to 
%       continue in the experiment.

arguments
    expParams (1,1) struct = struct()
end

%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------


fprintf(outputDivider())

fprintf("Initializing the experiment: Agregation with delay - Fokker-Planck.\n\n")

fprintf(outputDivider())


function result = IsAllModInteger(x)
    result = IsAllInteger(x) && all((x > 0) | (x == -1), "all");
end

function result = IsAllNonnegInteger(x)
    result = IsAllInteger(x) && all(x >= 0, "all");
end

function result = IsAllPosInteger(x)
    result = IsAllInteger(x) && all(x > 0, "all");
end

function result = IsAllInteger(x)
    result = IsAllReal(x) && all(x == floor(x), "all");
end

function result = IsAllNonnegReal(x)
    result = IsAllReal(x) && all(x >= 0, "all");
end

function result = IsAllPosReal(x)
    result = IsAllReal(x) && all(x > 0, "all");
end

function result = IsAllReal(x)
    result = isnumeric(x) && isreal(x) && all(isfinite(x), "all");
end

function result = ArrayDim(x)
    result = ndims(x);
    if isvector(x)
        result = 1;
    end
end


% Set or initialize experiment parameters

%---------------------------------RNG---------------------------------------

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

%-----------------------SPACE-&-INITIAL-CONDITIONS------------------------

gridDims = [];
d = [];

% Domain dimensions
if ~isfield(expParams, "domainDims") || ~isvector(expParams.domainDims) || ...
        ~IsAllPosReal(expParams.domainDims)
    fprintf("domainDims - Either no or wrong value for the vector of dimensions.\n")
    domainDims = [];     % default dimensions will be set later
    fprintf("Setting all domain dimensions to 1.\n\n")
else    
    domainDims = double(expParams.domainDims(:).');
    d = length(domainDims);
    fprintf("domainDims:")
    disp(domainDims)
end

% Initial density, grid dimensions and dimension
if ~isfield(expParams,"rho0") || ~IsAllReal(expParams.rho0) || isempty(expParams.rho0) || ...
        (~isempty(d) && ~isequal(ArrayDim(expParams.rho0), d))
    fprintf("rho0 - Either no or wrong value for the matrix of initial density.\n")
    fprintf("Initializing random experiment.\n\n");
    rho = []; % we will set the initial density later as a random vector
% Found rho0
else
    rho = double(expParams.rho0);
    if isrow(rho)
        rho = rho.'; % Accept any vector input but work with column vectors
    end
    if isvector(rho)
        gridDims = numel(rho);
    else
        gridDims = size(rho);
    end
    d = numel(gridDims);
    rho = reshape(rho, numel(rho), 1); % make it into vector
    fprintf("rho0 accepted, d = %i, gridDims:", d);
    disp(gridDims)
end

% Grid dimensions (if not inferred)
if isempty(gridDims)
    if ~isfield(expParams, "gridDims") || ~isvector(expParams.gridDims) || ...
            ~IsAllPosInteger(expParams.gridDims)
        fprintf("gridDims - Either no or wrong value for the grid dimensions.\n") 
        fprintf("   |\n")
        gridDims = [];
        % Set later after d
    else
        gridDims = double(expParams.gridDims(:).'); % Make it row vector
        d = numel(gridDims);
        fprintf("d = %i, gridDims:", d)
        disp(gridDims)
    end
end

% Dimension (if not inferred)
if isempty(d)
    if ~isfield(expParams,"d") || ~isscalar(expParams.d) || ~IsAllPosInteger(expParams.d)
        fprintf("   d - Either no or wrong value for the dimension.\n")
        d = 1;                  % default dimension of the space
        fprintf("   Setting d = %i.\n", d);
    else
        d = double(expParams.d);
        fprintf("   d = %i.\n", d)
    end
end

% Set defaul gridDims
if isempty(gridDims)
    switch d
        case 1 
            gridDims = 10000;
        case 2 
            gridDims = [100, 100];
        case 3 
            gridDims = [20, 20, 20];
        otherwise 
            gridDims = ones(1, d) * 10;
    end
    fprintf("   |")
    fprintf("Setting gridDims:")
    disp(gridDims)
end

% Set default domainDims
if isempty(domainDims)
    domainDims = ones(1, d);    % default domainDims
end
volume = prod(domainDims);

% Grid point count
N = prod(gridDims);

% Initial fluctuations
if ~isfield(expParams, "maxInitFlucts") || ~isscalar(expParams.maxInitFlucts) || ...
        ~IsAllPosReal(expParams.maxInitFlucts)
    fprintf("maxInitFlucts - Either no or wrong value for the maximal initial fluctuations.\n")
    maxInitFlucts = 1e-3;           % default maximal initial fluctuations
    fprintf("Setting maxInitFlucts = %.3e.\n\n", maxInitFlucts)
else 
    maxInitFlucts = double(expParams.maxInitFlucts);
    fprintf("maxInitFlucts = %.3e\n\n", maxInitFlucts)
end

% Initial history of density
if ~isempty(rho)
    if ~isfield(expParams,"rhoInitHist") || ~IsAllReal(expParams.rhoInitHist) || ...
            ~isequal(size(expParams.rhoInitHist, 1:d), gridDims) || ...
            ~isequal(size(expParams.rhoInitHist, d+1), 0)
        fprintf("rhoInitHist - Either no or wrong value for the matrix of initial densities.\n")
        rhoHist = [];  % default initial history
        %rhoHist = genRhoInitHist(maxInitFlucts, rho, stepDelay);  % default initial history
        fprintf("Initializing experiment with random initial density history.\n\n")
    else
        rhoHist = double(expParams.rhoInitHist);
        rhoHist = reshape(rhoHist, [N, size(expParams.rhoInitHist, d+1)]); % make it into columns of a matrix
        fprintf("rhoInitHist accepted.\n\n")
    end
else % It does not make sense to have a history for unknown initial conditions
    rhoHist = [];
    fprintf("rhoHist ignored, since 'rho0' was not specified.\n\n")
end

%----------------------------------MODEL------------------------------------

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

% Response decay
if ~isfield(expParams, "rd") || ~isscalar(expParams.rd) || ~IsAllReal(expParams.rd)
    fprintf("rd - Either no or wrong value for the response decay.\n")
    rd = 1;
    fprintf("Setting rd = %.3e\n\n", rd)
else
    rd = double(expParams.rd);
    fprintf("rd = %.3e\n\n", rd)
end

% Boundary conditions
definedBoundConds = ["NoBoundary", "Periodic", "Reflective"];
if ~isfield(expParams, "boundConds") || ~isscalar(expParams.boundConds) || ...
        ~isstring(expParams.boundConds) || ~ismember(expParams.boundConds, definedBoundConds)     
    fprintf("boundConds - Either no or wrong value for the boundary conditions.\n")
    boundConds = "Periodic";    % default boundary conditions
    fprintf("Setting boundConds = '%s'.\n\n", boundConds)
else
    boundConds = expParams.boundConds;
    fprintf("boundConds: '%s'.\n\n", boundConds)
end

%------------------------------TIME-&-DELAY---------------------------------

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
        stepCount = 2000;                    % default number of time steps
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
if isempty(rhoHist) % = wrong or no rhoInitHist input

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
% We got correct rhoInitHist
else
    stepDelay = size(rhoHist, 2); % it got reshaped earlier
    fprintf("Inferring step delay from 'rhoInitHist': stepDelay = %i (tau = %.3e).\n\n", ...
        stepDelay, stepDelay * dt)
end

% Make default rhoInitHist
if isempty(rhoHist)
    % Make default rho
    if isempty(rho)
        rho = ones(N, 1);
        rhoHist = genRhoInitHistSmooth(maxInitFlucts, rho, gridDims, stepDelay + 1);
        rhoHist = reshape(rhoHist, [N, stepDelay + 1]);
        rho = rhoHist(:,1);
        rhoHist = rhoHist(:,2:end);
    else
        rhoHist = genRhoInitHistSmooth(maxInitFlucts, rho, gridDims, stepDelay);
        rhoHist = reshape(rhoHist, [N, stepDelay]);
    end
end

%----------------------------------USER-------------------------------------

% Plot only for supported dimensions
if d <= 3

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

else
    stepPlotMod = -1;
    plotLastAgg = false;
    fprintf("Dimension larger than 3 detected. Forcing stepPlotMod = %i, plotLastAgg = %i.\n\n", stepPlotMod, plotLastAgg)
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

% Construct matrices
if ~isfield(expParams, "A") || ~IsAllReal(expParams.A) || ...
        ~isequal(size(expParams.A), [N,N]) || ...
        ~isfield(expParams,"distM") || ~IsAllNonnegReal(expParams.distM) || ...
        ~isequal(size(expParams.distM), [N,N])
    constrMats = true;
else
    constrMats = false;
end

% Video path
if ~isfield(expParams, "videoPath") || ~isscalar(expParams.videoPath) || ...
        ~isstring(expParams.videoPath) || strlength(expParams.videoPath) == 0
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


%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------


if waitForConf
    fprintf(outputDivider())
    fprintf("Press any key to confirm input.\n\n")
    pause
end


% Set auxiliary variables -------------------------------------------------

fprintf(outputDivider())

fprintf("Preparing inner variables.\n\n")

% Coefficient to access history
histCoeff = stepDelay;

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

% Step length
dDim = domainDims ./ gridDims;

maxPlot = 0;

% Make new figure if simulation is plotted and plot initial densisy
if stepPlotMod > 0
    figure
    plotSimStep(rho, rhoHist(:,histCoeff), 0)
end


% Set output variables ----------------------------------------------------


fprintf(outputDivider())

fprintf("Preparing output variables.\n\n")


% Auxiliary function to determine the count of to be recorded steps
function count = getRecCount(module)
    if module > 0
        count = ceil(stepCount / module);
        if count == 0
            count = 1; % We always record the last step
        end
    else
        count = 1; % Record just the last step
    end
end

% Decide the size of rhoRec
rhoRecCount = getRecCount(stepRecMod);

% Setup to record initial step
if recInitStep
    fprintf("Recording initial step.\n\n")
    rhoRecCount = rhoRecCount + 1; % this is used later, so we update it to match how many steps we record
    rhoRec = zeros([N, rhoRecCount]);
    rhoRec(:,1) = rho;
    rhoRecIndex = 2;
else 
    fprintf("NOT recording initial step.\n\n")
    rhoRec = zeros([N, rhoRecCount]);
    rhoRecIndex = 1;
end


% Pre-calculations --------------------------------------------------------

fprintf(outputDivider())

if constrMats
    fprintf("A, distM - Either no or wrong value for either the discretization matrix or the distance matrix.\n")
    fprintf("Constructing default discretization and distance matrices.\n")
    fprintf("(This can take some time for larger systems and higher dimensions, so it is advised to precalculate the matrices.)\n\n")
    [A, distM] = makeFPMatricesOpt(gridDims, boundConds, domainDims, makeStencil(d, gridDims));
else
    A = double(expParams.A);
    distM = double(expParams.distM);
    fprintf("'A' & 'distM' accepted.\n\n");
end

% Tol for roundoffs
tol = 100 * sum(eps(domainDims)) / d;

intRadTol = intRad - tol;

% Normalization multiplier
dbv = discretizedBallVolume(intRadTol, dDim);

% Pre-calculate the kernel W, normalized such that its integral is 1
% Remark: Don't use simple 'W = (dist <  intRad)' since this gives issues
% with roundoff errors (staircasing) when intRad is a multiple of dx
W = sparse(double(distM < intRadTol));   % strict cutoff, robust near multiples of dx
W = W / dbv;


%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------


fprintf(outputDivider())
fprintf("Starting the simulation")
printTitle()
fprintf(".\n\n")


% Solve for t = 1:stepCount -----------------------------------------------
% This loop calculates the step t from the step t-1 -----------------------

for t = 1:stepCount
    
    rhoDelayed = rhoHist(:,histCoeff);

    % Store present rho to the rho history
    rhoHist(:,histCoeff) = rho;

    % Update hist coeff by -1
    histCoeff = mod(histCoeff - 2, stepDelay) + 1;

    % Calculate G
    WConvRhoDelayed = W * rhoDelayed;
    G_sqrd = volume / 2 * exp(-rd * WConvRhoDelayed) .^ 2;

    % Create the matrix for semi-implicit finite differences
    M = speye(N,N) + dt * A .* (G_sqrd.');
    
    % Make the step
    rho = M \ rho;
    
    % Plot
    if stepPlotMod > 0 && mod(t, stepPlotMod) == 0
        plotSimStep(rho, rhoDelayed, t)
    end

    % Record simulation step
    if t < stepCount % Last step is recorded after this loop (since we could get mod(stepCount,stepRecMod) ~= 0, so we record it specially after the loop)
        if stepRecMod ~= -1 && mod(t, stepRecMod) == 0
            rhoRec(:,rhoRecIndex) = rho;
            rhoRecIndex = rhoRecIndex + 1;
        end
    end
end

function plotSimStep(rho, rhoDelayed, t)
    maxPlot = max(maxPlot, max(rho) * getMaxPlotMultip());
    title(sprintf("T = %d", t))
    switch d
        case 1
            grid1D = ((1:gridDims).' - 0.5) * dDim;
            clf
            hold on
            plot(grid1D, rhoDelayed)
            plot(grid1D, rho)
            hold off
            axis([0, domainDims, 0, maxPlot])
            legend("rho delayed", "rho")
        case 2
            imagesc(reshape(rho, gridDims))
            set(gca, 'YDir', 'normal');
            clim([0, maxPlot])
            colorbar("Ticks", linspace(0, maxPlot, 5))
        case 3
            gridX = ((1:gridDims(1)) - 0.5) * dDim(1);
            gridY = ((1:gridDims(2)) - 0.5) * dDim(2);
            gridZ = ((1:gridDims(3)) - 0.5) * dDim(3);
            rho50prc = mean(rho);
            rho95prc = prctile(rho,95);
            clf
            colormap("turbo")
            isosurface(gridX, gridY, gridZ, reshape(rho, gridDims), rho50prc);
            isosurface(gridX, gridY, gridZ, reshape(rho, gridDims), 1);
            isosurface(gridX, gridY, gridZ, reshape(rho, gridDims), rho95prc);
            clim([0, maxPlot])
            colorbar("Ticks", linspace(0, maxPlot, 5))
        otherwise
            error('Dimension d = %.i not implemented.', d);
    end
    
    title(sprintf("T = %d", t * dt))

    f = getframe;

    if recordVideo
        writeVideo(vidWriter, f)
    end
end


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


% Finish ------------------------------------------------------------------

% Close video writer
if recordVideo
    close(vidWriter)
end

% Record final simulation step (the final step might not have been possible to record in the loop)
rhoRec(:,end) = rho;
rhoRec = reshape(rhoRec, [gridDims, rhoRecCount]);

% Record last history of rho - just permute the history
% The coefficient was decreased in the final step, so we need to increase it back
lastIndex = histCoeff + 1; 
if lastIndex >= stepDelay
    lastIndex = 1;
end
permutation = [lastIndex:stepDelay, 1:lastIndex-1];
rhoHist = rhoHist(:,permutation);

fprintf(outputDivider())
fprintf("Simulation")
printTitle()
fprintf(" finished.\n\n")

function printTitle()
    if isfield(expParams,"expTitle") && isstring(expParams.expTitle) && isscalar(expParams.expTitle)
    fprintf(" titled %s", expParams.expTitle)
    end
end

% Plot aggregation
if stepPlotMod ~= -2
    fprintf(outputDivider())
    fprintf("Plotting agregation groups.\n\n")
    plotSimStep(rho, rhoHist(:, stepDelay), stepCount) % TODO - agregation detection
end


end