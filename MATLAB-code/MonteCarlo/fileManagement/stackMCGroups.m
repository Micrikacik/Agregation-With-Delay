function [] = stackMCGroups(expsParams, folderPathFunc, filePostfixFunc, startGroup, endGroup, exceptFields, infoFile, fileName)

arguments
    expsParams (:,1) struct
    folderPathFunc function_handle
    filePostfixFunc function_handle
    startGroup double {mustBeInteger, mustBePositive} = []
    endGroup double {mustBeInteger, mustBePositive} = []
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

    if  isempty(endGroup)
        endGroup = groupCount;
    end
    if isempty(startGroup)
        startGroup = 1;
    end

    stackedResults = struct([]); % start as an empty struct array
    time = 0;
    group_params = [];

    for group = startGroup:endGroup
        postfix = filePostfixFunc(params, i_exp) + "_" + MCGroupFilePostfix(group, groupCount);
        filePath = MCFilePath(folderPath, fileName, postfix);
        experiment = load(filePath);

        % Check if parameters between groups are the same
        new_group_params = experiment.params;
        new_group_params = rmfield(new_group_params, "expTitle");
        if ~isempty(group_params) && ~isequal(group_params, new_group_params)
            fprintf("There are different experiment parameters between groups!\n")
            compareStructs(new_group_params, group_params, "group_" + group + "_i_exp_" + i_exp, "group_" + group-1 + "_i_exp_" + i_exp)
        end

        % Check if parameters of the group are the same as the ones from
        % the input
        group_params = new_group_params;
        if ~isequal(group_params, params)
            fprintf("Group %i of i_exp %i might have different experiment parameters from the input!\n", group, i_exp)
            compareStructs(group_params, params, "group_" + group + "_i_exp_" + i_exp, "input")
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
    fprintf("Saved file %i of %i.\n\n", i_exp, count);
end

close(wBar)
