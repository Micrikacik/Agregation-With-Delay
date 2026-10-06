% Starts default parpool and based on poolsize runs MC simulations to make
% simCount simulations of 2D smallstep experiment.

clearvars

d = 2;
tau = 0.3;

[expParams, folderFun, fileFun] = makeDelayExpsParamsSmallstep(d, tau);
MonteCarloManager(expParams, folderFun, fileFun);