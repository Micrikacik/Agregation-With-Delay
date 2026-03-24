function [] = MonteCarloManager(d, stepDelays, delayType, nMC, startGroup, endGroup, overwrite, fileName)

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
    startGroup = 1
    endGroup = 0
    overwrite = false
    fileName = "MCData"
end

T = 1e+6;

baseParams = struct( ...
        ...% No rng controll
        ...% No random generator settings controll
        ...% No initial position controll
        "N", 400, ...% FIXED
        "d", d, ...% FIXED
        ...% Default intRad
        ...% Default boundConds
        "T", T, ...% FIXED
        "dt", 1e-3, ...% FIXED
        "delayType", delayType, ...% FIXED
        ...% VARIABLE step delay - will be set later to make multiple params
        ...% No initial history controll
        "waitForConf", false, ...% FIXED
        "stepPlotMod", -2, ...% FIXED
        ...% No agent marking
        ...% No color for agent marking
        "stepRecMod", ceil(T/min(T,100)), ...% FIXED - record either all or up to 100 steps
        "recInitStep", true, ...% FIXED
        ...% Same rec mod for theta as for x
        ...% Record initial theta if recording initial x
        "thetaOccurMod", 1 ...
        ...% No experiment title
    );

p = gcp("nocreate");
if isempty(p)
    p = parpool; % There is no parpool, so we start the default one
    poolsize = p.NumWorkers;
else
    poolsize = p.NumWorkers;
end

groupCount = ceil(nMC / poolsize);
if endGroup < 1
    endGroup = groupCount;
end

% Divide the run into groups, which have the same size as there is workers,
% to be able to save a part of the results as soon as possible
for group = startGroup:endGroup
    % Run through all of the given stepDelay values
    for stepDelay = stepDelays
        postfix = MCFilePostfix(delayType,d,stepDelay,group,groupCount);

        params = baseParams;
        params.stepDelay = stepDelay; % VARIABLE step delay
        params.expTitle = postfix;

        [results, time] = MonteCarlo(params,poolsize);
        
        folderPath = MCFolderPath(delayType,d,stepDelay,true);
        path = MCFilePath(folderPath,fileName,postfix);
        if isfile(path) && (overwrite == false)
            fprintf("File on path \n%s \nalready exists, and overwriting is not allowed.", path)
        else 
            save(path, 'results', 'params', "poolsize", "time");
        end

        infoPath = sprinf("%s/info", folderPath);
        if isfile(infoPath) && (overwrite == false)
            fprintf("File on path \n%s \nalready exists, and overwriting is not allowed.", infoPath)
        else 
            save(infoPath, "poolsize", "groupCount", "params");
        end
    end
end