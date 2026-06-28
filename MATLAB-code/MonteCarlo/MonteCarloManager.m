function [] = MonteCarloManager(expsParams, nMC, folderPathFunc, filePostfixFunc, startGroup, endGroup, overwrite, useBaseSeed, fileName)

% expsParams (structure array) - contains parameters to the experiments
% nMC (positive integer) - number of Monte Carlo simulations
% folderPathFunc & filePostfixFunc (function handle) - functions with 
%   experiment parameters and experiment index (i_exp = 1:length(expsParams)) as
%   input, which returns folderPath to where the results will be saved and filePostfix
% TODOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
% startGroup (positive integer)
% endGroup (positive integer)

arguments
    expsParams (:,1) struct
    nMC (1,1) double {mustBeInteger, mustBePositive}
    folderPathFunc function_handle
    filePostfixFunc function_handle
    startGroup double {mustBeInteger, mustBePositive} = []
    endGroup double {mustBeInteger, mustBePositive} = []
    overwrite (1,1) logical = false
    useBaseSeed (1,1) logical = false
    fileName (1,1) string = "MCData"
end

for i_exp = 1:length(expsParams)
    % Disable plotting, if not specified
    if ~isfield(expsParams(i_exp), "stepPlotMod") || isempty(expsParams(i_exp).stepPlotMod)
        expsParams(i_exp).stepPlotMod = -2;
    end
    % Disable recording, if not specified
    if ~isfield(expsParams(i_exp), "stepRecMod") || isempty(expsParams(i_exp).stepRecMod)
        expsParams(i_exp).stepRecMod = -1;
    end
    % Disable recording of initial step, if not specified
    if ~isfield(expsParams(i_exp), "recInitStep") || isempty(expsParams(i_exp).recInitStep)
        expsParams(i_exp).recInitStep = false;
    end
    % Disable waiting for confirmation, if not specified
    if ~isfield(expsParams(i_exp), "waitForConf") || isempty(expsParams(i_exp).waitForConf)
        expsParams(i_exp).waitForConf = false;
    end
end

p = gcp("nocreate"); % Get active parpool
if isempty(p)
    p = parpool; % There is no parpool, so we start the default one
    poolsize = p.NumWorkers;
else
    poolsize = p.NumWorkers;
end

groupCount = ceil(nMC / poolsize);
% Set default values for 'endGroup' and 'startGroup'
if  isempty(endGroup)
    endGroup = groupCount;
end
if isempty(startGroup)
    startGroup = 1;
end

fprintf("Experiment parameters (values will be displayed when simulations start):\n")
disp(expsParams)
fprintf("Number of Monte Carlo simulations: %i\n\n", nMC)
fprintf("Folder path function:\n")
disp(folderPathFunc)
fprintf("File postfix function:\n")
disp(filePostfixFunc)
fprintf("Group count: %i\n\n", groupCount)
fprintf("Pool size: %i\n\n", poolsize)
fprintf("Start group: %i\n\n", startGroup)
fprintf("End group: %i\n\n", endGroup)
fprintf("Overwrite: %i\n\n", overwrite)
fprintf("Use base seed: %i\n\n", useBaseSeed)
fprintf("File name: %s\n\n", fileName)

% Save info files with poolsize and groupcount for later processing
for i_exp = 1:length(expsParams)
    params = expsParams(i_exp);
    folderPath = folderPathFunc(params, i_exp);

    % Create the directory if it does not exist
    if ~exist(folderPath,"dir")
        mkdir(folderPath);
    end

    infoPath = sprintf("%s/info.mat", folderPath);
    if isfile(infoPath) && (overwrite == false)
        fprintf("File on path \n%s \nalready exists, and overwriting is not allowed. Saving as 'DUPLICATE'.\n\n", infoPath)
        infoPath = sprintf("%s/DUPLICATE_info.mat", folderPath);
    end 
    save(infoPath, "poolsize", "groupCount");
end

% Divide the run into groups, which have the same size as there is workers,
% to be able to save a part of the results as soon as possible
for group = startGroup:endGroup
    % Run through all of the given experiment parameters
    for i_exp = 1:length(expsParams)
        params = expsParams(i_exp);
        postfix = filePostfixFunc(params,i_exp) + "_" + MCGroupFilePostfix(group,groupCount);
        params.expTitle = postfix;
        
        if useBaseSeed
            baseSeed = (group-1) * poolsize;
            % Run parallel simulations with baseSeed
            [results, time] = MonteCarlo(params, poolsize, baseSeed);
        else
            % Run parallel simulations without baseSeed
            [results, time] = MonteCarlo(params, poolsize);
        end

        folderPath = folderPathFunc(params,i_exp);
        filePath = MCFilePath(folderPath,fileName,postfix);

        if isfile(filePath) && (overwrite == false)
            fprintf("File on path \n%s \nalready exists, and overwriting is not allowed. Saving as 'DUPLICATE'.\n\n", filePath)
            filePath = MCFilePath(folderPath,"DUPLICATE_" + fileName,postfix);
        end
        save(filePath, 'results', 'params', "time");
    end
end