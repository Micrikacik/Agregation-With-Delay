function xData = gatherMCData(N, d, nMC, stepDelays, delayType, iRec, fileName)

arguments
    N 
    d 
    nMC 
    stepDelays 
    delayType
    iRec = -1
    fileName = "MCData"
end

count = numel(stepDelays);

xData = zeros(N,d,nMC,count);

wBar = waitbar(0,"Gathering data...");

for i = 1:count
    stepDelay = stepDelays(i);
    folder = MCFolderPath(delayType,d,stepDelay);
    postfix = MCFilePostfix(delayType,d,stepDelay);
    file = MCFilePath(folder,fileName,postfix);
    results = load(file).results;
    xRec = results.xRec;
    len = size(xRec,3);
    if iRec == -1
        index = len;
    else
        index = mod(iRec - 1, len) + 1;
    end
    for i_Sim = 1:nMC
        xData(:, :, i_Sim, i) = results(i_Sim).xRec(:,:,index);
    end
    waitbar(i / count,wBar,"Gathering data...")
end

close(wBar)
