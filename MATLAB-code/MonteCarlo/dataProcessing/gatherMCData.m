function xData = gatherMCData(expsParams, folderPathFunc, filePostfixFunc, nMC, iRec, fileName)

% TODO
% iRec (positive integer or -1) - index of the record to gather

arguments
    expsParams (:,1) struct
    folderPathFunc function_handle
    filePostfixFunc function_handle
    nMC (1,1) double {mustBeInteger, mustBePositive}
    iRec (1,1) double {mustBeInteger} = -1
    fileName (1,1) string = "MCData"
end

count = length(expsParams);
% There might be different values of N or d
maxN = max([expsParams.N]);
maxd = max([expsParams.d]);

% Preallocate matrix
xData = zeros(maxN, maxd, nMC, count);

wString = "Gathering data ...";
wBar = waitbar(0, wString);

for i_exp = 1:count
    params = expsParams(i_exp);
    folder = folderPathFunc(params, i_exp);
    postfix = filePostfixFunc(params, i_exp);
    file = MCFilePath(folder, fileName, postfix);

    results = load(file).results;

    if nMC > length(results)
        fprintf("The desired value of the input 'nMC' is larger than the actual amount of simulations for experiment %i.\n", i_exp)
    end 

    for i_sim = 1:min(nMC, length(results))
        xRec = results(i_sim).xRec;

        len = size(xRec,3);
        if iRec == -1
            index = len;
        else
            index = mod(iRec - 1, len) + 1;
        end

        xData(1:params.N, 1:params.d, i_sim, i_exp) = results(i_sim).xRec(:,:,index);
    end

    waitbar(i_exp / count, wBar, wString)
end

close(wBar)
