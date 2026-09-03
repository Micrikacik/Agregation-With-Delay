function [expsParams, folderPathFunc, filePostfixFunc] = makeDelayExpsParamsSmallstep(d, tau)

arguments
    d (1,1) {mustBeInteger}
    tau (1,1)
end

expsParams = struct( ...
        ...% No rng controll
        ...% No random generator settings controll
        ...% No initial position controll
        "N", 400, ...% FIXED
        "d", d, ...% FIXED
        ...% Default intRad
        ...% Default boundConds
        "T", 1e+3, ...% FIXED
        "dt", 1e-5, ...% FIXED
        "delayType", "Reaction", ...% FIXED
        "tau", tau, ...% FIXED
        ...% No initial history controll
        "waitForConf", false, ...% FIXED
        "stepPlotMod", -2 ...% FIXED
        ...% No agent marking
        ...% No color for agent marking
        ...% No position recording
        ...% No initial position recording
        ...% Same rec mod for theta as for x
        ...% Record initial theta if recording initial x
        ...% No position recording
        ...% No experiment title
    );

folderPathFunc = @(params, i_exp) smallstepMCFolderPath(params.delayType, params.d, round(params.tau / params.dt));
filePostfixFunc = @(params, i_exp) smallstepMCFilePostfix(params.delayType, params.d, round(params.tau / params.dt));

function folderPath =  smallstepMCFolderPath(delayType, d, stepDelay)
    folderPath = sprintf("MonteCarlo/delayType%s/%iD/smallstep/stepDelay%i", delayType, d, stepDelay);
end

function filePostfix =  smallstepMCFilePostfix(delayType, d, stepDelay)
    filePostfix = sprintf("%s_%iD_smallstep_%i", delayType, d, stepDelay);
end

end