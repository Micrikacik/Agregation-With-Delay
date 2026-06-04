%% Backup original files
clearvars

d = 2;
delayType = "Reaction";
fileName = "MCData";
backupFileName = "BackupMCData";
stepDelays = [0,90:30:300];

for stepDelay = stepDelays
    path = MCFilePath(MCFolderPath(delayType,d,stepDelay),fileName,MCFilePostfix(delayType,d,stepDelay));
    backupPath = MCFilePath(MCFolderPath(delayType,d,stepDelay),backupFileName,MCFilePostfix(delayType,d,stepDelay));
    copyfile(path,backupPath)
    fprintf("Backuped %i\n", stepDelay)
end

%% Load backup
clearvars

d = 2;
delayType = "Reaction";
fileName = "MCData";
backupFileName = "BackupMCData";
stepDelays = 120;

for stepDelay = stepDelays
    path = MCFilePath(MCFolderPath(delayType,d,stepDelay),fileName,MCFilePostfix(delayType,d,stepDelay));
    backupPath = MCFilePath(MCFolderPath(delayType,d,stepDelay),backupFileName,MCFilePostfix(delayType,d,stepDelay));
    copyfile(backupPath,path)
    fprintf("Loaded %i\n", stepDelay)
end

%% Get bad indices
clearvars

d = 2;
delayType = "Reaction";
fileName = "MCData";
repairFileName = "RepairData";
stepDelays = [120];
badIndices = -ones(length(stepDelays),52);

for s = 1:length(stepDelays)
    stepDelay = stepDelays(s);
    path = MCFilePath(MCFolderPath(delayType,d,stepDelay),fileName,MCFilePostfix(delayType,d,stepDelay));
    results = load(path).results;
    fprintf("Checking %i\n", stepDelay)
    same = 0;
    for i = 1:length(results)
        for j = i+1:length(results)
            if norm(results{i}.xRec(:,:,1) - results{j}.xRec(:,:,1), "fro") < 0.001
                fprintf("(%i,%i) ", i, j)
                same = same + 1;
                badIndices(s,same) = i;
                break
            end
        end
    end
    fprintf("\nfinished %i, same: %i\n", stepDelay, same)
end

% Replace bad data with repaired data

for s = 1:length(stepDelays)
    stepDelay = stepDelays(s);
    path = MCFilePath(MCFolderPath(delayType,d,stepDelay),fileName,MCFilePostfix(delayType,d,stepDelay));
    load(path)
    repairPath = MCFilePath(MCFolderPath(delayType,d,stepDelay),repairFileName,MCFilePostfix(delayType,d,stepDelay));
    repairedResults = load(repairPath).results;
    fprintf("Step delay %i\n", stepDelay)
    for b = 1:length(repairedResults)
        index = badIndices(s,b);
        if index > 0
            fprintf("%i,", index)
            results{index} = repairedResults{b};
        else
            break % there are no more bad indices for this step delay
        end
    end
    fprintf(" saving repaired data ... ")
    save(path,"results","params","time")
    fprintf("done\n\n")
end