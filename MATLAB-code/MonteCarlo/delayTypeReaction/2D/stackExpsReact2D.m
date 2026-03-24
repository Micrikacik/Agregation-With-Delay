clearvars
clc

delayType = "Reaction";
d = 2;
stepDelays = [0,60:30:300];
fileName = "MCData";

for stepDelay = stepDelays
    folderPath = MCFolderPath(delayType, 2, stepDelay);
    info = load(sprintf("%s/info", folderPath));
    poolsize = info.poolsize;
    groupCount = info.groupCount;
    stackedResults = {}; % start as an empty cell
    time = 0;
    for group = 1:groupCount
        postfix = MCFilePostfix(delayType,d,stepDelay,group,groupCount);
        filePath = MCFilePath(folderPath,fileName,postfix);
        experiment = load(filePath);
        params = experiment.params;
        time = time + experiment.time;
        results = experiment.results;
        stackedResults = [stackedResults, results];
        fprintf("%i",group)
    end
    results = stackedResults; % to make the name in file correct
    newPostfix = MCFilePostfix(delayType,d,stepDelay);
    newFilePath = MCFilePath(folderPath,fileName,newPostfix);
    save(newFilePath, "params", "results", "time")
    fprintf(" %i done\n", stepDelay)
end
