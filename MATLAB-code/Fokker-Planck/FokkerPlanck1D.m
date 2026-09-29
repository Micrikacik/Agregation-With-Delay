function [rhoRec, rhoHist, rngSetts] = FokkerPlanck1D(expParams)

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
%           random generator settings. This struct is also an output of
%           this function, and is created after the simulation is done.
%           This field's main goal is to enable user to continue in an
%           experiment which already finished, by using the values it
%           returned as initial conditions and random generator settings.
%
%       SPACE & INITIAL CONDITIONS:
%       rho0 (positive float COLUMN vector) - vector of initial discretized density. 
%           Alternatively, set number of grid points 'gridPointCount'.
%       gridPointCount (integer scalar) - number of points at which we
%           disretize the differential equation using uniform step.
%       L (positive float) - length of the simulation interval [0, L].
%       maxInitFluc (positive float scalar) - if either 'rho0'
%           'rhoInitHist' is not provided, this value is used to generate
%           them with random fluctuations, which are at most 'maxInitFluc'.
%
%       MODEL:
%       intRad (positive float) - radius of interactions between agents.
%       boundConds (string) - type of the boundary conditions.
%           Must be one of the following strings:
%               "Periodic"
%               "Reflective"
%       respDecay (positive float scalar) - it is a multiplier 'a' in
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
%   rngSetts (struct) - struct returned by the rng function, containing
%       random generator settings right after the simulation have finished.
%       It can be directly plugged in as an input to this function to 
%       continue in the experiment.


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

% Initial density
if ~isfield(expParams,"rho0") || ~isfloat(expParams.rho0) || isempty(expParams.rho0)
    fprintf("Either no or wrong value for the matrix of initial density 'rho0'.\n")
    fprintf("Looking for the input value for 'gridPointCount'.\n")
    fprintf('   |\n')
    if ~isfield(expParams,"gridPointCount") || ~IsInteger(expParams.gridPointCount) || ...
            expParams.gridPointCount <= 0 || ~isequal(size(expParams.gridPointCount),[1,1])
        fprintf("   Either no or wrong value for the grid point count 'gridPointCount'.\n")
        gridPointCount = 400;                % default grid point count
        fprintf("   Setting gridPointCount = %i.\n", gridPointCount);
    else
        gridPointCount = single(expParams.gridPointCount);
        fprintf("   gridPointCount = %i.\n", gridPointCount)
    end
    fprintf('   |\n')
    fprintf("Initializing random experiment with gridPointCount = %i.\n\n", gridPointCount);
    rho = []; % we will set the initial density later as a random vector
else
    rho = expParams.rho0;
    gridPointCount = size(rho, 1);
    fprintf("rho0 accepted, gridPointCount = %i.\n\n", gridPointCount)
end

% Dimension
d = 1;

% Interval length
if ~isfield(expParams,"L") || ~isfloat(expParams.L) || expParams.L <= 0 || ...
        ~isequal(size(expParams.L),[1,1])
    fprintf("Either no or wrong value for the interval length 'L'.\n")
    L = 1;           % default length
    fprintf("Setting interval length to 1.\n\n")
else 
    L = expParams.L;
    fprintf("L = %.3d\n\n", L)
end

% Initial fluctuations
if ~isfield(expParams,"maxInitFluc") || ~isfloat(expParams.maxInitFluc) || ...
        expParams.maxInitFluc <= 0 || ~isequal(size(expParams.L),[1,1])
    fprintf("Either no or wrong value for the maximal initial fluctuations 'maxInitFluc'.\n")
    maxInitFluc = 1e-2;           % default maximal initial fluctuations
    fprintf("Setting maxInitFluc = %.3d.\n\n", maxInitFluc)
else 
    maxInitFluc = expParams.maxInitFluc;
    fprintf("maxInitFluc = %.3d\n\n", maxInitFluc)
end

% Setting random initial density
if isempty(rho)
    fluctuations = makeFluctuations(maxInitFluc, gridPointCount, d);
    rho = 1 + fluctuations;
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

% Response decay
if ~isfield(expParams,"respDecay") || ~isnumeric(expParams.respDecay) || expParams.respDecay < 0
    fprintf("Either no or wrong value for the response decay 'respDecay'.\n")
    respDecay = 1;
    fprintf("Setting respDecay = %.3d\n\n", respDecay)
else
    respDecay = expParams.respDecay;
    fprintf("respDecay = %.3d\n\n", respDecay)
end

%------------------------------TIME-&-DELAY---------------------------------

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

% Step count
if ~isfield(expParams,"stepCount") || ~IsInteger(expParams.stepCount) || expParams.stepCount < 0 || ...
        ~isequal(size(expParams.stepCount),[1,1])
    fprintf("Either no or wrong value for the number of time steps 'stepCount'.\n")
    if ~isfield(expParams,"T") || ~isfloat(expParams.T) || expParams.T <= 0 || ...
        ~isequal(size(expParams.T), [1,1])
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
if ~isfield(expParams,"stepDelay") || ~IsInteger(expParams.stepDelay) || expParams.stepDelay < 0 || ...
        ~isequal(size(expParams.stepDelay),[1,1])
    fprintf("Either no or wrong value for the step delay 'stepDelay'.\n")
    if ~isfield(expParams,"tau") || ~isfloat(expParams.tau) || expParams.tau < 0 || ...
        ~isequal(size(expParams.tau),[1,1])
        fprintf("Either no or wrong value for the delay 'tau'.\n")
        stepDelay = 5;                    % default number of time steps
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
        ~isequal(size(expParams.rhoInitHist), [gridPointCount,stepDelay])
    fprintf("Either no or wrong value for the matrix of initial densities 'rhoInitHist'.\n")
    rhoHist = genRhoInitHist(maxInitFluc, gridPointCount, stepDelay, L, d);  % default initial history
    fprintf("Initializing experiment with random initial density history.\n\n")
else
    rhoHist = expParams.rhoInitHist;
    fprintf("rhoInitHist accepted.\n\n")
end

%----------------------------------USER-------------------------------------

% Step plot mod
if ~isfield(expParams,"stepPlotMod") || ~IsInteger(expParams.stepPlotMod) || ...
        ((expParams.stepPlotMod <= 0) && expParams.stepPlotMod ~= -1 && expParams.stepPlotMod ~= -2) || ...
        ~isequal(size(expParams.stepPlotMod),[1,1])
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
        ~isequal(size(expParams.stepRecMod),[1,1])
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


fprintf("----------------------------------\n\n")
fprintf("Starting the simulation")
if isfield(expParams,"expTitle") && isstring(expParams.expTitle) && isequal(size(expParams.expTitle),[1,1])
    fprintf(", title: %s", expParams.expTitle)
end
fprintf(".\n\n")


% Set auxiliary variables

% Coefficient to access history
histCoeff = stepDelay;

% Setup VideoWriter if the simulation is recorded
recordVideo = false;
if isfield(expParams,"recordVideoPath") && isstring(expParams.recordVideoPath)
    recordVideo = true;
    vidWrit = VideoWriter(sprintf("%s.avi", expParams.recordVideoPath));
    frameTime = stepPlotMod * dt;
    vidWrit.FrameRate = 1 / frameTime;
    open(vidWrit)
end

% Step length
dx = L / gridPointCount;

% Coefficient in finite differences
dtxx = dt / (dx^2);

% Equidistant grid
plotGrid = (1:gridPointCount)*dx;


% Set output variables

% Auxiliary function to determine the count of to be recorded steps
function count = getRecCount(module)
    % We do not want to record
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
    rhoRec = zeros([gridPointCount, rhoRecCount + 1]);
    rhoRec(:,1) = rho;
    rhoRecIndex = 2;
else 
    rhoRec = zeros([gridPointCount, rhoRecCount]);
    rhoRecIndex = 1;
end


% Pre-calculations

% Tol for roundoffs
tol = 100 * eps(L);

% Normalization multiplier
multip = 1 / (1 + 2 * floor(((intRad - tol) - dx) / 2 / dx));

sum(rho)

% Pre-calculate the distances over the periodic/reflected domain
switch boundConds 

    case "Periodic"
        k = 0:gridPointCount-1;
        dist = min(k, gridPointCount-k) * dx;

        % Pre-calculate the kernel W, normalized such that its integral is 1
        % Remark: Don't use simple 'W = (dist <  intRad)' since this gives issues
        % with roundoff errors (staircasing) when intRad is a multiple of dx
        w = double(dist < intRad - tol);   % strict cutoff, robust near multiples of dx
        w = w * multip;
        
        W = zeros(gridPointCount, gridPointCount);
        for i = 1:gridPointCount
            W(i,:) = circshift(w, i-1); % w is 1 by gridPointCount
        end

    case "Reflective"
        % Pre-calculate the kernel W, normalized such that its integral is 1
        % Remark: Don't use simple 'W = (dist <  intRad)' since this gives issues
        % with roundoff errors (staircasing) when intRad is a multiple of dx      
        W = zeros(gridPointCount, gridPointCount);
        for i = 1:gridPointCount
            k = 1-i:gridPointCount-i;
            dist = abs(k) * dx;   
            
            w = double(dist < intRad - tol);   % strict cutoff, robust near multiples of dx
            multip = 1 / round((intRad - tol) / dx);
            w = w * multip;

            W(i,:) = w;
        end

    otherwise
        error("Undefined boundary conditions: '%.i'", boundConds);
end


% Make new figure if simulation is plotted and plot initial density
if stepPlotMod > 0
    figure
    plotSimStep(rho, 0)
end


% Solve for t = 1:stepCount
for t = 1:stepCount
    
    if stepDelay > 0
        rhoDelayed = rhoHist(:,histCoeff);

        % Store present rho to the rho history
        rhoHist(:,histCoeff) = rho;

        % Update hist coeff
        histCoeff = histCoeff - 1;
        histCoeff = mod(histCoeff - 1, stepDelay) + 1;
    else
        rhoDelayed = rho;
    end

    % Calculate G
    WConvRhoDelayed = W * rhoDelayed;
    G_sqrd = exp(-respDecay * WConvRhoDelayed); % midstep
    G_sqrd = G_sqrd.^2;
    
    %calculate the L^2 norm
    %E(t)=sum((rho.^2).*G);
        
    % Compose the matrix for semi-implicit finite differences
    A = diag(1+dtxx*G_sqrd) - 0.5*dtxx*diag(G_sqrd(2:gridPointCount),1) - 0.5*dtxx*diag(G_sqrd(1:gridPointCount-1),-1);
    
    % Adjust to BCs
    switch boundConds
        case "Periodic"
            A(1,gridPointCount) = -0.5*dtxx*G_sqrd(gridPointCount);
            A(gridPointCount,1) = -0.5*dtxx*G_sqrd(1);
        case "Reflective"
            A(1,2) = -0.5*dtxx*G_sqrd(2);
            A(gridPointCount,gridPointCount-1) = -0.5*dtxx*G_sqrd(gridPointCount-1);
        otherwise
            error("Undefined boundary conditions: '%.i'", boundConds);
    end

    % Make the step
    rho = A \ rho;
    
    %print the total mass to check mass conservation
    %sum(rho)
        
    
    % Plot - to make correct 1D plot, we need current theta
    if stepPlotMod > 0 && mod(t,stepPlotMod) == 0
        plotSimStep(rho, t)
    end

    % Record simulation step
    if t < stepCount % Last step is recorded after this loop (since we could get mod(stepCount,stepRecMod) ~= 0, so we record it specially after the loop)
        if stepRecMod ~= -1 && mod(t,stepRecMod) == 0
            rhoRec(:,:,rhoRecIndex) = rho;
            rhoRecIndex = rhoRecIndex + 1;
        end
    end
end

function plotSimStep(rho, t)
    f = [];
    switch d
        case 1
            % Plot rho
            %subplot(2,1,1);
            plot(plotGrid, rho); axis([0 L 0 max(4, max(rho + 0.5))]);
            ttl = sprintf('rho; physical time = %e', t * dt);
            title(ttl);

            % Plot G_sqrd and convolution
            %subplot(2,1,2);
            %plot(plotGrid, WConvRhoDelayed, plotGrid, G_sqrd,'--');
            %title('Value of G^2 & value of W convolution with delayed rho')
            %legend('convolution', 'G^2')

            f = getframe;
            
            %pause to see the IC
            %if (t==1) pause; end
            
            %fname = sprintf('num_mf1/num_mf1_%f.eps',phys_time);
            %saveas(gcf, fname);
            %fname = sprintf('num_mf1/num_mf1_%f.fig',phys_time);
            %saveas(gcf, fname); 
        otherwise
            error('Dimension d = %.i not implemented.', d);
    end

    if recordVideo && ~isempty(f)
        writeVideo(vidWrit,f)
    end
end


% Finish

% Close video writer
if recordVideo
    close(vidWrit)
end

% Record final simulation step (the final step might not have been possible to record in the loop)
rhoRec(:,:,end) = rho;

% Record last history of rho - just permute the history
% The coefficient was decreased in the final step, so we need to increase it back
lastIndex = histCoeff + 1; 
if lastIndex >= stepDelay
    lastIndex = 1;
end
permutation = [lastIndex:stepDelay, 1:lastIndex-1];
rhoHist = rhoHist(:,permutation);

% Return random generator settings
rngSetts = rng;

fprintf("----------------------------------\n\n")
fprintf("Simulation")
if isfield(expParams,"expTitle") && isstring(expParams.expTitle) && isequal(size(expParams.expTitle),[1,1])
    fprintf(" titled %s", expParams.expTitle)
end
fprintf(" finished.\n\n")

if stepPlotMod ~= -2
    fprintf("----------------------------------\n\n")
    fprintf("Plotting agregation groups.\n\n")
    plotSimStep(rho, stepCount) % TODO - agregation detection
end


end