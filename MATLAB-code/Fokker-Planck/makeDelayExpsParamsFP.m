function [expsParams, folderPathFunc, filePostfixFunc] = makeDelayExpsParamsFP(d, stepDelays)

% d - dimension
% stepDelays - vector of step delays
% delayType - delay type
% nMC - number of Monte Carlo simulations
% fileName - file name

arguments
    d
    stepDelays
end

stepCount = 1e+6;

baseParams = struct( ...
        "rngSeed", 1, ...
        ...% No random generator settings controll
        ...% No initial density controll
        ...% Default gridDims
        "d", d, ...% FIXED
        ...% Default intRad
        ...% Default boundConds
        "stepCount", stepCount, ...% FIXED
        "dt", 1e-3, ...% FIXED
        ...% VARIABLE step delay - will be set later to make multiple expParams
        ...% No initial history controll
        ...% Default max fluctuations
        "waitForConf", false, ...% FIXED
        "stepPlotMod", -1, ...% FIXED
        "stepRecMod", ceil(stepCount/min(stepCount, 100)), ...% FIXED - record either all or up to 100 steps
        "recInitStep", true ...% FIXED
        ...% No experiment title
    );

expsParams = struct( ...
        "rngSeed", [], ...
        "d", cell(length(stepDelays), 1), ... % Set the length of the structure array
        "stepCount", [], ...
        "dt", [], ...
        "waitForConf", [], ...
        "stepPlotMod", [], ...
        "stepRecMod", [], ...
        "recInitStep", [] ...
    );

% Set all base expParams
for i = 1:length(stepDelays)
    expsParams(i) = baseParams;
end

% Set step delays
for i = 1:length(stepDelays)
    stepDelay = stepDelays(i);
    expsParams(i).stepDelay = stepDelay;
end

folderPathFunc = @(expParams, i_exp) FPFolderPath(expParams.d, expParams.stepDelay);
filePostfixFunc = @(expParams, i_exp) FPFilePostfix(expParams.d, expParams.stepDelay);