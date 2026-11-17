function [rho, rhoInitHist] = genRhoInitHist(maxFlucts, rho, stepDelay)

% maxFluc - maximal fluctuations added in one step
% rho0 - vectorized density
% stepDelay - step delay

N = numel(rho);

stepDelayPlus1 = stepDelay + 1;

histDims = [N, stepDelayPlus1];
rhoInitHist = makeFlucts(maxFlucts, histDims);
rhoInitHist = rhoInitHist + rho;
rho = rhoInitHist(:,1);
rhoInitHist = rhoInitHist(:,2:end);