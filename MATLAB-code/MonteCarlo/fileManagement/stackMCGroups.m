function [] = stackMCGroups(d, stepDelays, delayType, exceptFields, infoFile, fileName)

arguments
    d 
    stepDelays 
    delayType 
    exceptFields (:,1) string = []
    infoFile (1,1) string = "info"
    fileName (1,1) string = "MCData"
end

for stepDelay = stepDelays
    folderPath = MCFolderPath(delayType, d, stepDelay);
    info = load(sprintf("%s/%s.mat", folderPath, infoFile));
    poolsize = info.poolsize;
    groupCount = info.groupCount;
    stackedResults = struct([]); % start as an empty struct array
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
        for field = exceptFields
            results = rmfield(results,field);
        end
        stackedResults = [  stackedResults;
                            results         ];
        fprintf("%i",group)
    end
    fprintf(" saving file ...")
    results = stackedResults; % to make the name in file correct
    newPostfix = MCFilePostfix(delayType,d,stepDelay);
    newFilePath = MCFilePath(folderPath,fileName,newPostfix);
    save(newFilePath, "params", "results", "time")
    fprintf(" stepDelay %i done\n", stepDelay)
end
