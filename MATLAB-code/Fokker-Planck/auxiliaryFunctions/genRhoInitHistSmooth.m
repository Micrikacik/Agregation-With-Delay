function [rho, rhoInitHist] = genRhoInitHistSmooth(maxFlucts, rho, gridDims, stepDelay, centralFreqs)

% maxFluc - maximal fluctuations added in one step
% rho - vectorized density to make fluctuations from
% gridDims - dimension of the grid used to discretize rho (without vectorization)
% stepDelay - step delay

N = numel(rho);

stepDelayPlus1 = stepDelay + 1;

histDims = [gridDims, stepDelayPlus1];
rhoInitHist = reshape(makeFluctsSmooth(maxFlucts, histDims, centralFreqs), [N, stepDelayPlus1]);
rhoInitHist = rhoInitHist + rho;
rho = rhoInitHist(:,1);
rhoInitHist = rhoInitHist(:,2:end);

