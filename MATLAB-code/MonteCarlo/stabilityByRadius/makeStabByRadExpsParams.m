function [expsParams, folderPathFuncHandle, filePostfixFuncHandle] = makeStabByRadExpsParams(d, intRads)

    stepCount = 1e+6;
    
    baseParams = struct( ...
            ...% No rng controll
            ...% No random generator settings controll
            ...% No initial position controll
            "N", 400, ...% FIXED
            "d", d, ...% FIXED
            ...% Variable intRad
            ...% Default boundConds
            "stepCount", stepCount, ...% FIXED
            "dt", 1e-3, ...% FIXED
            "delayType", "None", ...% FIXED
            ...% No stepDelay
            ...% No initial history controll
            "waitForConf", false, ...% FIXED
            "stepPlotMod", -2, ...% FIXED
            ...% No agent marking
            ...% No color for agent marking
            "stepRecMod", ceil(stepCount/min(stepCount, 100)), ...% FIXED - record either all or up to 100 steps
            "recInitStep", true, ...% FIXED
            ...% Same rec mod for theta as for x
            ...% Record initial theta if recording initial x
            "thetaOccurMod", 1 ...
            ...% No experiment title
        );
    
    expsParams = struct( ...
            "N", cell(length(intRads), 1), ... % Set the length of the structure array
            "d", [], ...
            "stepCount", [], ...
            "dt", [], ...
            "delayType", [], ...
            "waitForConf", [], ...
            "stepPlotMod", [], ...
            "stepRecMod", [], ...
            "recInitStep", [], ...
            "thetaOccurMod", [] ...
        );
    
    % Set all base params
    for i = 1:length(intRads)
        expsParams(i) = baseParams;
    end
    
    % Set interaction radii
    for i = 1:length(intRads)
        intRad = intRads(i);
        expsParams(i).intRad = intRad;
    end
    
    function path = folderPathFunc(d, intRad)
        path = "MonteCarlo/stabilityByRadius";
        path = sprintf("%s/%iD/intRad0p%.0f", path, d, intRad * 100);
    end
    
    function postfix = filePostfixFunc(d, intRad)
        postfix = sprintf("%iD_0p%.0f", d, intRad * 100);
    end
    
    folderPathFuncHandle = @(params, i_exp) folderPathFunc(params.d, params.intRad);
    filePostfixFuncHandle = @(params, i_exp) filePostfixFunc(params.d, params.intRad);
end