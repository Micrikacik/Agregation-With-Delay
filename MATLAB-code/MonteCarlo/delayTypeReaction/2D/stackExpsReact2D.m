clearvars
clc

delayType = "Reaction";
d = 2;
% stepDelays = [0,60:30:420];
stepDelays = [0,90:30:300];
% fileName = "MCData";
fileName = "RepairData";
% infoFile = "info";
infoFile = "DUPLICATE_info";

for stepDelay = stepDelays
    folderPath = MCFolderPath(delayType, 2, stepDelay);
    info = load(sprintf("%s/%s.mat", folderPath, infoFile));
    poolsize = info.poolsize;
    groupCount = info.groupCount;
    stackedResults = {}; % start as an empty cell
    time = 0;
    params = 0;
    for group = 1:groupCount
        postfix = MCFilePostfix(delayType,d,stepDelay,group,groupCount);
        filePath = MCFilePath(folderPath,fileName,postfix);
        experiment = load(filePath);
        if ~isequal(params, 0) && ~isequal(rmfield(params, "expTitle"), rmfield(experiment.params, "expTitle"))
            fprintf("\nThere are different experiment parameters between groups!\n")
        end
        params = experiment.params;
        time = time + experiment.time;
        results = experiment.results;
        stackedResults = [stackedResults; results];
        fprintf("%i",group)
    end
    fprintf(" saving file ...")
    results = stackedResults; % to make the name in file correct
    newPostfix = MCFilePostfix(delayType,d,stepDelay);
    newFilePath = MCFilePath(folderPath,fileName,newPostfix);
    save(newFilePath, "params", "results", "time")
    fprintf(" stepDelay %i done\n", stepDelay)
end
