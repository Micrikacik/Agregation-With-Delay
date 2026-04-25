function [results, time] = MonteCarlo(expParams, nMC, baseSeed)

% expParams - Parameters used as input to aggWithDelaySpeedup.m
% nMC - Number of Monte Carlo simulations.
% baseSeed - Number to which the worker index is added (baseSeed + i) and
%   the result is used as the seed for rng in the worker as rng(baseSeed + i). 
%   If 'baseSeed = []', then rng(baseSeed + i) is not used in the worker.
%
% Calls rng('shuffle') before starting the parallel simulations.

arguments
    expParams 
    nMC 
    baseSeed = []
end

if isfield(expParams,"rngSetts") || isfield(expParams,"rngSeed")
    fprintf("There are rng settings present in experiment parameters! We shall remove them.")
    expParams = rmfield(expParams, {'rngSetts', 'rngSeed'});
end

results = struct( ...
    "xRec", cell(nMC,1), ...
    "thetaRec", cell(nMC,1), ...
    "thetaOccur", cell(nMC,1), ...
    "xInitHist", cell(nMC,1), ...
    "xHist", cell(nMC,1), ...
    "rngSetts", cell(nMC,1), ...
    "initRngState", cell(nMC,1) ...
    );

tic
parfor i = 1:nMC
    params = expParams;
    params.expTitle = params.expTitle + sprintf("_worker%i", i);
    if ~isempty(baseSeed)
        params.rngSeed = baseSeed + i; % Set the random seed for the worker
    end
    initRngState = rng;
    [xRec, thetaRec, thetaOccur, xInitHist, xHist, rngSetts] = aggWithDelaySpeedup(params)
    results(i) = struct("xRec", xRec, "thetaRec", thetaRec, "thetaOccur", thetaOccur, "xInitHist", xInitHist, "xHist", xHist, "rngSetts", rngSetts, "initRngState", initRngState)
end
time = toc;