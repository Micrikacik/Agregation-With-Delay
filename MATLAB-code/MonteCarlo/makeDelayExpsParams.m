function [expsParams, folderPathFunc, filePostfixFunc] = makeDelayExpsParams(d, stepDelays, delayType)

% d - dimension
% stepDelays - vector of step delays
% delayType - delay type
% nMC - number of Monte Carlo simulations
% fileName - file name

arguments
    d
    stepDelays
    delayType
end

stepCount = 1e+6;

baseParams = struct( ...
        ...% No rng controll
        ...% No random generator settings controll
        ...% No initial position controll
        "N", 400, ...% FIXED
        "d", d, ...% FIXED
        ...% Default intRad
        ...% Default boundConds
        "stepCount", stepCount, ...% FIXED
        "dt", 1e-3, ...% FIXED
        "delayType", delayType, ...% FIXED
        ...% VARIABLE step delay - will be set later to make multiple expParams
        ...% No initial history controll
        "waitForConf", false, ...% FIXED
        "stepPlotMod", -1, ...% FIXED
        ...% No agent marking
        ...% No color for agent marking
        "stepRecMod", ceil(stepCount/min(stepCount, 100)), ...% FIXED - record either all or up to 100 steps
        "recInitStep", true, ...% FIXED
        ...% Same rec mod for theta as for x
        ...% Record initial theta if recording initial x
        "thetaOccurMod", 1 ...
        ...% No experiment title
    );

expsParams = struct( ...
        "N", cell(length(stepDelays), 1), ... % Set the length of the structure array
        "d", [], ...
        "stepCount", [], ...
        "dt", [], ...
        "delayType", [], ...
        "waitForConf", [], ...
        "stepPlotMod", [], ...
        "stepRecMod", [], ...
        "recInitStep", [], ...
        "thetaOccurMod", [] ...
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

folderPathFunc = @(expParams,i_exp) MCFolderPath(expParams.delayType, expParams.d, expParams.stepDelay);
filePostfixFunc = @(expParams,i_exp) MCFilePostfix(expParams.delayType, expParams.d, expParams.stepDelay);