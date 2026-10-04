% Starts default parpool and based on poolsize runs MC simulations to make
% simCount simulations of 2D smallstep experiment.

clearvars

simCount = 22; % I have 6x13=78 already out of 100

p = gcp("nocreate"); % Get active parpool or start a new one
poolsize = p.NumWorkers;

startGroup = 7; % I already have groups 1 to 6, so we start at group 7
groupCount = ceil(simCount / poolsize);
endGroup = startGroup + groupCount - 1; % endGroup is inclusive, so we subtract 1

d = 1;
tau = 0.09;

[expParams, folderFun, fileFun] = makeDelayExpsParamsSmallstep(d, tau);
MonteCarloManager(expParams, folderFun, fileFun, startGroup, endGroup);