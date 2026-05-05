function [] = MonteCarloManager(d, stepDelays, delayType, nMC, startGroup, endGroup, overwrite, useBaseSeed, fileName)

% d - dimension
% stepDelays - vector of step delays
% delayType - delay type
% nMC - number of Monte Carlo simulations
% fileName - file name

arguments
    d
    stepDelays
    delayType
    nMC
    startGroup = []
    endGroup = []
    overwrite = false
    useBaseSeed = false
    fileName = "MCData"
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
        ...% VARIABLE step delay - will be set later to make multiple params
        ...% No initial history controll
        "waitForConf", false, ...% FIXED
        "stepPlotMod", -2, ...% FIXED
        ...% No agent marking
        ...% No color for agent marking
        "stepRecMod", ceil(stepCount/min(stepCount,100)), ...% FIXED - record either all or up to 100 steps
        "recInitStep", true, ...% FIXED
        ...% Same rec mod for theta as for x
        ...% Record initial theta if recording initial x
        "thetaOccurMod", 1 ...
        ...% No experiment title
    );

rng('default')
rng('shuffle')

p = gcp("nocreate");
if isempty(p)
    p = parpool; % There is no parpool, so we start the default one
    poolsize = p.NumWorkers;
else
    poolsize = p.NumWorkers;
end

groupCount = ceil(nMC / poolsize);
if  isempty(endGroup)
    endGroup = groupCount;
end
if isempty(startGroup)
    startGroup = 1;
end

fprintf("Base parameters:\n")
disp(baseParams)
fprintf("Step delays:")
disp(stepDelays)
fprintf("Number of Monte Carlo simulations: %i\n\n", nMC)
fprintf("Group count: %i\n\n", groupCount)
fprintf("Pool size: %i\n\n", poolsize)
fprintf("Start group: %i\n\n", startGroup)
fprintf("End group: %i\n\n", endGroup)
fprintf("Overwrite: %i\n\n", overwrite)
fprintf("Use base seed: %i\n\n", useBaseSeed)
fprintf("File name: %s\n\n", fileName)

digitCountOfNMC = floor(log10(nMC + poolsize - 1)) + 1;
baseSeed = [];

% Divide the run into groups, which have the same size as there is workers,
% to be able to save a part of the results as soon as possible
for group = startGroup:endGroup
    % Run through all of the given stepDelay values
    for stepDelay = stepDelays
        postfix = MCFilePostfix(delayType,d,stepDelay,group,groupCount);

        params = baseParams;
        params.stepDelay = stepDelay; % VARIABLE step delay
        params.expTitle = postfix;
        
        if useBaseSeed
            baseSeed = stepDelay * 10^digitCountOfNMC + (group-1) * poolsize;
        end
        [results, time] = MonteCarlo(params, poolsize, baseSeed);
        
        folderPath = MCFolderPath(delayType,d,stepDelay,true);
        path = MCFilePath(folderPath,fileName,postfix);
        if isfile(path) && (overwrite == false)
            fprintf("File on path \n%s \nalready exists, and overwriting is not allowed. Saving as 'DUPLICATE'.\n\n", path)
            path = MCFilePath(folderPath,"DUPLICATE_" + fileName,postfix);
        end
        save(path, 'results', 'params', "poolsize", "time");
    end
end

% Save info files
for stepDelay = stepDelays
    folderPath = MCFolderPath(delayType,d,stepDelay);
    infoPath = sprintf("%s/info.mat", folderPath);
    if isfile(infoPath) && (overwrite == false)
        fprintf("File on path \n%s \nalready exists, and overwriting is not allowed. Saving as 'DUPLICATE'.\n\n", infoPath)
        infoPath = sprintf("%s/DUPLICATE_info.mat", folderPath);
    end 
    save(infoPath, "poolsize", "groupCount");
end