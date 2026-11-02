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


fprintf("----------------------------------\n\n")

fprintf("Initializing the experiment: Agregation with delay - Fokker-Planck.\n\n")

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
            ~isscalar(expParams.rngSeed)
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

% Initial density, grid dimensions and dimension
if ~isfield(expParams,"rho0") || ~isfloat(expParams.rho0) || isempty(expParams.rho0)
    fprintf("Either no or wrong value for the matrix of initial density 'rho0'.\n")
    fprintf("Looking for the input value for 'gridDims'.\n")
    fprintf('   |\n')
    % Look for gridDims
    if ~isfield(expParams,"gridDims") || ~IsInteger(expParams.gridDims) || ...
            any(expParams.gridDims <= 0) || ~isrow(expParams.gridDims)
        fprintf("   Either no or wrong value for the grid dimensions 'gridDims'.\n")
        fprintf("   Looking for the input value for 'd'.\n")
        fprintf('       |\n')
        % Look for d
        if ~isfield(expParams,"d") || ~isfloat(expParams.d) || ...
            expParams.d <= 0 || ~isscalar(expParams.d) || ~IsInteger(expParams.d)
            fprintf("       Either no or wrong value for the dimension 'd'.\n")
            d = 1;                          % default d
            gridDims = 10000;                 % default gridDims
            fprintf("       Setting 1D simulation with gridDims = %i.\n", gridDims);
        % Found d
        else
            d = expParams.d;
            switch d
                case 1 
                    gridDims = 10000;
                case 2 
                    gridDims = [100, 100];
                case 3 
                    gridDims = [30, 30, 30];
                otherwise 
                    gridDims = ones(1, d) * 10;
            end
            fprintf("       d = %i.\n", d)
        end
        fprintf('       |\n')
        fprintf('   -----\n')
    % Found gridDims
    else
        gridDims = expParams.gridDims;
        d = numel(gridDims);
        fprintf("   gridDims:")
        disp(gridDims)
    end
    fprintf('   |\n')
    fprintf("Initializing random experiment with gridDims:");
    disp(gridDims)
    rho = []; % we will set the initial density later as a random vector
% Found rho0
else
    rho = expParams.rho0;
    if isvector(rho)
        gridDims = numel(rho);
    else
        gridDims = size(rho);
    end
    d = numel(gridDims);
    rho = reshape(rho, numel(rho), 1); % make it into vector
    fprintf("rho0 accepted, gridDims:");
    disp(gridDims)
end

% Grid point count
N = prod(gridDims);

% Domain dimensions
if ~isfield(expParams,"domainDims") || ~isfloat(expParams.domainDims) || ...
        any(expParams.domainDims <= 0) || ~isequal(size(expParams.domainDims), size(gridDims))
    fprintf("Either no or wrong value for the domain dimensions 'domainDims'.\n")
    domainDims = ones(size(gridDims));           % default domain domainDims
    fprintf("Setting domainDims:")
    disp(domainDims)
else 
    domainDims = expParams.domainDims;
    fprintf("domainDims:")
    disp(domainDims)
end

% Domain volume
volume = prod(domainDims);

% Initial fluctuations
if ~isfield(expParams,"maxInitFlucts") || ~isfloat(expParams.maxInitFlucts) || ...
        expParams.maxInitFlucts <= 0 || ~isscalar(expParams.L)
    fprintf("Either no or wrong value for the maximal initial fluctuations 'maxInitFlucts'.\n")
    maxInitFlucts = 1e-3;           % default maximal initial fluctuations
    fprintf("Setting maxInitFlucts = %.3d.\n\n", maxInitFlucts)
else 
    maxInitFlucts = expParams.maxInitFlucts;
    fprintf("maxInitFlucts = %.3d\n\n", maxInitFlucts)
end

% Setting random initial density
if isempty(rho)
    flucts = makeFlucts(maxInitFlucts, N);
    rho = 1 / volume + flucts; % it is density, so integral must be 1
end

%----------------------------------MODEL------------------------------------

% Interaction radius
if ~isfield(expParams,"intRad") || ~isnumeric(expParams.intRad) || expParams.intRad < 0
    fprintf("Either no or wrong value for the interaction radius 'intRad'.\n")
    intRad = getIntRad(d);
    fprintf("Setting intRad = %.3d\n\n", intRad)
else
    intRad = expParams.intRad;
    fprintf("intRad = %.3d\n\n", intRad)
end

% Response decay
if ~isfield(expParams,"rd") || ~isnumeric(expParams.rd) || ...
        ~isscalar(expParams.rd)
    fprintf("Either no or wrong value for the response decay 'rd'.\n")
    rd = 1;
    fprintf("Setting rd = %.3d\n\n", rd)
else
    rd = expParams.rd;
    fprintf("rd = %.3d\n\n", rd)
end

% Boundary conditions
if ~isfield(expParams,"boundConds") || ~isstring(expParams.boundConds) || ...
        ~isscalar(expParams.boundConds)
    fprintf("Either no or wrong value for the boundary conditions 'boundConds'.\n")
    boundConds = "Periodic";    % default boundary conditions
    fprintf("Setting boundary conditions to '%s'.\n\n", boundConds)
else
    boundConds = expParams.boundConds;
    fprintf("boundConds: '%s'.\n\n", boundConds)
end

%------------------------------TIME-&-DELAY---------------------------------

% Time step length
if ~isfield(expParams,"dt") || ~isfloat(expParams.dt) || expParams.dt <= 0 || ...
        ~isscalar(expParams.dt)
    fprintf("Either no or wrong value for the time step length 'dt'.\n")
    dt = 1e-3;                  % default time step length
    fprintf("Setting dt = %.3d\n\n", dt)
else
    dt = expParams.dt;
    fprintf("dt = %.3d\n\n", dt)
end

% Step count
if ~isfield(expParams,"stepCount") || ~IsInteger(expParams.stepCount) || ...
        expParams.stepCount < 0 || ~isscalar(expParams.stepCount)
    fprintf("Either no or wrong value for the number of time steps 'stepCount'.\n")
    if ~isfield(expParams,"T") || ~isfloat(expParams.T) || expParams.T <= 0 || ...
        ~isscalar(expParams.T)
        fprintf("Either no or wrong value for the total simulation time 'T'.\n")
        stepCount = 1000;                    % default number of time steps
        fprintf("Setting stepCount = %i (T = %.3d).\n\n", stepCount, stepCount * dt)
    else
        stepCount = round(expParams.T / dt);
        T = stepCount * dt;
        fprintf("T = %.3d (stepCount = %i).\n\n", T, stepCount)
    end
else
    stepCount = expParams.stepCount;
    fprintf("stepCount = %i (T = %.3d).\n\n", stepCount, stepCount * dt)
end

% Step delay
if ~isfield(expParams,"stepDelay") || ~IsInteger(expParams.stepDelay) || ...
        expParams.stepDelay < 0 || ~isscalar(expParams.stepDelay)
    fprintf("Either no or wrong value for the step delay 'stepDelay'.\n")
    if ~isfield(expParams,"tau") || ~isfloat(expParams.tau) || expParams.tau < 0 || ...
        ~isscalar(expParams.tau)
        fprintf("Either no or wrong value for the delay 'tau'.\n")
        stepDelay = 60;                    % default step delay
        fprintf("Setting stepDelay = %i (tau = %.3d).\n\n", stepDelay, stepDelay * dt)
    else
        stepDelay = round(expParams.tau / dt);
        tau = stepDelay * dt;
        fprintf("tau = %.3d (stepDelay = %i).\n\n", tau, stepDelay)
    end
else
    stepDelay = expParams.stepDelay;
    fprintf("stepDelay: %i (tau = %.3d).\n\n", stepDelay, stepDelay * dt)
end

% Initial history of density
if ~isfield(expParams,"rhoInitHist") || ~isfloat(expParams.rhoInitHist) || ...
        ~isequal(size(expParams.rhoInitHist), [gridDims, stepDelay])
    fprintf("Either no or wrong value for the matrix of initial densities 'rhoInitHist'.\n")
    rhoHist = genRhoInitHist(maxInitFlucts, rho, stepDelay);  % default initial history
    fprintf("Initializing experiment with random initial density history.\n\n")
else
    rhoHist = reshape(expParams.rhoInitHist, [N, stepDelay]); % make it into columns of a matrix
    fprintf("rhoInitHist accepted.\n\n")
end

%----------------------------------USER-------------------------------------

% Step plot mod
if ~isfield(expParams,"stepPlotMod") || ~IsInteger(expParams.stepPlotMod) || ...
        ((expParams.stepPlotMod <= 0) && expParams.stepPlotMod ~= -1 && expParams.stepPlotMod ~= -2) || ...
        ~isscalar(expParams.stepPlotMod)
    fprintf("Either no or wrong value for the step plot mod 'stepPlotMod'.\n")
    stepPlotMod = 5;  % default step plot mod
    fprintf("Setting step plot mod to %i.\n\n", stepPlotMod)
else
    stepPlotMod = expParams.stepPlotMod;
    fprintf("stepPlotMod: %i.\n\n", stepPlotMod)
end

% Step record mod
if ~isfield(expParams,"stepRecMod") || ~IsInteger(expParams.stepRecMod) || ... 
        (expParams.stepRecMod <= 0 && expParams.stepRecMod ~= -1) || ...
        ~isscalar(expParams.stepRecMod)
    fprintf("Either no or wrong value for the step record mod 'stepRecMod'.\n")
    stepRecMod = -1;  % default step record mod
    fprintf("Setting step record mod to %i.\n\n", stepRecMod)
else
    stepRecMod = expParams.stepRecMod;
    fprintf("stepRecMod: %i.\n\n", stepRecMod)
end


%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------


if ~isfield(expParams, "waitForConf") || expParams.waitForConf == true
    fprintf("----------------------------------\n\n")
    fprintf("Press space to start the simulation.\n\n")
    pause
end


% Set auxiliary variables -------------------------------------------------

fprintf("----------------------------------\n\n")
fprintf("Preparing variables.\n\n")

% Coefficient to access history
histCoeff = stepDelay;

% Setup VideoWriter if the simulation is recorded
recordVideo = false;
if isfield(expParams,"recordVideoPath") && isstring(expParams.recordVideoPath)
    recordVideo = true;
    vidWriter = VideoWriter(sprintf("%s.avi", expParams.recordVideoPath));
    frameTime = stepPlotMod * dt;
    vidWriter.FrameRate = 1 / frameTime;
    open(vidWriter)
end

% Step length
dDim = domainDims ./ gridDims;


% Set output variables ----------------------------------------------------

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
if isfield(expParams,"recInitStep") && expParams.recInitStep == true
    rhoRecCount = rhoRecCount + 1; % this is used later
    rhoRec = zeros([N, rhoRecCount]);
    rhoRec(:,1) = rho;
    rhoRecIndex = 2;
else 
    rhoRec = zeros([N, rhoRecCount]);
    rhoRecIndex = 1;
end


% Pre-calculations --------------------------------------------------------

fprintf("----------------------------------\n\n")

if ~isfield(expParams,"A") || ~isfloat(expParams.A) || ~isequal(size(expParams.A), [N,N]) || ...
        ~isfield(expParams,"distM") || ~isfloat(expParams.distM) || ~isequal(size(expParams.distM), [N,N])
    fprintf("Either no or wrong value for either the discretization matrix 'A' or the distance matrix 'distM'.\n")
    fprintf("Constructing default discretization and distance matrices.\n")
    fprintf("(This can take some time for larger systems and higher dimensions, so it is advised to precalculate the matrices.)\n\n")
    [A, distM] = makeFPMatricesAIOpt(gridDims, boundConds, domainDims, makeStencil(d, gridDims));
else
    A = expParams.A;
    distM = expParams.distM;
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

%sum(W,"all") / N
%spy(W)

% Make new figure if simulation is plotted and plot initial density
if stepPlotMod > 0
    figure
    plotSimStep(rho, rhoHist(:,histCoeff), 0)
end


% Solve for t = 1:stepCount -----------------------------------------------

fprintf("----------------------------------\n\n")
fprintf("Starting the simulation")
if isfield(expParams,"expTitle") && isstring(expParams.expTitle) && isscalar(expParams.expTitle)
    fprintf(", title: %s", expParams.expTitle)
end
fprintf(".\n\n")

for t = 1:stepCount
    
    if stepDelay > 0
        rhoDelayed = rhoHist(:,histCoeff);

        % Store present rho to the rho history
        rhoHist(:,histCoeff) = rho;

        % Update hist coeff by -1
        histCoeff = mod(histCoeff - 2, stepDelay) + 1;
    else
        rhoDelayed = rho;
    end

    % Calculate G
    WConvRhoDelayed = W * rhoDelayed;
    G_sqrd = 1 / 2 * exp(-rd * WConvRhoDelayed) .^ 2;

    % Create the matrix for semi-implicit finite differences
    M = speye(N,N) + dt * A .* (G_sqrd.');
    
    %calculate the L^2 norm
    %E(t)=sum((rho.^2).*G);
    
    % Make the step
    rho = M \ rho;
    
    %print the total mass to check mass conservation
    %sum(rho)
        
    
    % Plot
    if stepPlotMod > 0 && mod(t,stepPlotMod) == 0
        plotSimStep(rho, rhoDelayed, t)
    end

    % Record simulation step
    if t < stepCount % Last step is recorded after this loop (since we could get mod(stepCount,stepRecMod) ~= 0, so we record it specially after the loop)
        if stepRecMod ~= -1 && mod(t,stepRecMod) == 0
            rhoRec(:,rhoRecIndex) = rho;
            rhoRecIndex = rhoRecIndex + 1;
        end
    end
end

function plotSimStep(rho, rhoDelayed, t)
    maxPlot = max(1.2 * max(rho));
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
            rho80prc = prctile(rho,80);
            rho95prc = prctile(rho,95);
            clf
            colormap("turbo")
            isosurface(gridX, gridY, gridZ, reshape(rho, gridDims), rho50prc);
            isosurface(gridX, gridY, gridZ, reshape(rho, gridDims), rho80prc);
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


% Finish

% Close video writer
if recordVideo
    close(vidWriter)
end

% Record final simulation step (the final step might not have been possible to record in the loop)
rhoRec(:,end) = rho;
indices = num2cell(gridDims);
rhoRec = reshape(rhoRec, indices{:}, rhoRecCount);

% Record last history of rho - just permute the history
% The coefficient was decreased in the final step, so we need to increase it back
lastIndex = histCoeff + 1; 
if lastIndex >= stepDelay
    lastIndex = 1;
end
permutation = [lastIndex:stepDelay, 1:lastIndex-1];
rhoHist = rhoHist(:,permutation);

fprintf("----------------------------------\n\n")
fprintf("Simulation")
if isfield(expParams,"expTitle") && isstring(expParams.expTitle) && isscalar(expParams.expTitle)
    fprintf(" titled %s", expParams.expTitle)
end
fprintf(" finished.\n\n")

if stepPlotMod ~= -2
    fprintf("----------------------------------\n\n")
    fprintf("Plotting agregation groups.\n\n")
    plotSimStep(rho, rhoHist(:, stepDelay), stepCount) % TODO - agregation detection
end


end