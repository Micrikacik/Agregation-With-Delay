% Starts default parpool and based on poolsize runs MC simulations to make
% simCount simulations of 2D smallstep experiment.

clearvars

d = 1;
tau = 0.09;

[expParams, folderFun, fileFun] = makeDelayExpsParamsSmallstep(d, tau);
MonteCarloManager(expParams, folderFun, fileFun);