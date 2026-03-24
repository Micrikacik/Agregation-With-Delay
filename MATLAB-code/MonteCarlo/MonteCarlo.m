function [results, time] = MonteCarlo(expParams, nMC)

% params - parameters used as input to aggWithDelay.m
% nMC - number of Monte Carlo simulations

% Randomize randomness
rng('shuffle')

results = cell(nMC,1);

tic
parfor i = 1:nMC
    params = expParams;
    params.expTitle = params.expTitle + sprintf("_worker%i", i);
    [xRec, thetaRec, thetaOccur] = aggWithDelaySpeedup(params) % rng can not be controlled in parallel calculations
    results{i} = struct("xRec", xRec, "thetaRec", thetaRec, "thetaOccur", thetaOccur)
end
time = toc;