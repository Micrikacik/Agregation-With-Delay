clearvars
clc

N = 400;
d = 2;
nMC = 100;
stepDelays = 0:30:420;
count = numel(stepDelays);

xData = zeros(N,d,nMC,count);

for stepDelay = stepDelays
    folder = MCFolderPath("Reaction",d,stepDelay);
    postfix = MCFilePostfix("Reaction",d,stepDelay);
    file = MCFilePath(folder,"MCData",postfix);
    results = load(file).results;
    for i_Sim = 1:nMC
        xData(:, :, i_Sim, find(stepDelays == stepDelay)) = results{i_Sim}.xRec(:,:,end);
    end
end

processMCData(xData)
