d = 2;
delayType = "Reaction";
fileName = "MCData";



for stepDelay = 30:30:420
    path = MCFilePath(MCFolderPath(delayType,d,stepDelay),fileName,MCFilePostfix(delayType,d,stepDelay));
    load(path);
    fprintf("Started %i\n", stepDelay)
    res = struct("xRec",{},"thetaRec",{},"thetaOccur",{});
    for i = 1:length(results)
        res(i).xRec = results{i}.xRec;
        res(i).thetaRec = results{i}.thetaRec;
        res(i).thetaOccur = results{i}.thetaOccur;
    end
    fprintf("\nSaving %i\n", stepDelay)
    results = res';
    save(path, "params", "results", "time")
    fprintf("\nfinished %i\n", stepDelay)
end