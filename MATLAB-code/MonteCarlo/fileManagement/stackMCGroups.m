function [] = stackMCGroups(expsParams, folderPathFunc, filePostfixFunc, exceptFields, infoFile, fileName)

arguments
    expsParams (:,1) struct
    folderPathFunc function_handle
    filePostfixFunc function_handle
    exceptFields (:,1) string = []
    infoFile (1,1) string = "info"
    fileName (1,1) string = "MCData"
end

count = length(expsParams);
wBar = waitbar(0);

for i_exp = 1:count
    wString = sprintf("Stacking experiment %i of %i ...", i_exp, count);
    waitbar(0, wBar, wString)

    params = expsParams(i_exp);
    folderPath = folderPathFunc(params, i_exp);
    info = load(sprintf("%s/%s.mat", folderPath, infoFile));
    % poolsize = info.poolsize; % Not needed
    groupCount = info.groupCount;
    stackedResults = struct([]); % start as an empty struct array
    time = 0;
    group_params = 0;

    for group = 1:groupCount
        postfix = filePostfixFunc(params, i_exp) + "_" + MCGroupFilePostfix(group, groupCount);
        filePath = MCFilePath(folderPath, fileName, postfix);
        experiment = load(filePath);

        % Check if parameters between groups are the same
        if ~isequal(group_params, 0) && ~isequal(rmfield(group_params, "expTitle"), rmfield(experiment.params, "expTitle"))
            fprintf("\nThere are different experiment parameters between groups!\n")
        end

        % Check if parameters of the group are the same as the ones from
        % the input
        group_params = experiment.params;
        if ~isequal(rmfield(group_params, "expTitle"), params)
            fprintf("\nA group has different experiment parameters from the input!\n")
        end

        % Add results together
        time = time + experiment.time;
        results = experiment.results;
        for field = exceptFields
            results = rmfield(results,field);
        end
        stackedResults = [  stackedResults;
                            results         ];

        fprintf("Finished group %i.\n",group)
        waitbar(group / groupCount, wBar, wString)
    end

    wString = sprintf("Saving file %i of %i ...", i_exp, count);
    waitbar(1, wBar, wString)

    results = stackedResults; % to make the name in file correct
    newPostfix = filePostfixFunc(params, i_exp);
    newFilePath = MCFilePath(folderPath,fileName,newPostfix);
    save(newFilePath, "params", "results", "time")
    fprintf("Saved file %i of %i.\n", i_exp, count);
end

close(wBar)
