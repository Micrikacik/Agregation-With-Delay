% Starts default parpool and based on poolsize runs MC simulations to make
% simCount simulations of 2D smallstep experiment.

clearvars

simCount = 74; % I have 26 already out of 100

p = gcp("nocreate"); % Get active parpool or start a new one
poolsize = p.NumWorkers;

startGroup = 3; % I already have groups 1 & 2, so we start at group 3
groupCount = ceil(simCount / poolsize);
endGroup = startGroup + groupcount - 1; % endGroup is inclusive, so we subtract 1

d = 2;
tau = 0.3;

[expParams, folderFun, fileFun] = makeDelayExpsParamsSmallstep(d, tau);
MonteCarloManager(expParams, folderFun, fileFun, startGroup, endGroup);