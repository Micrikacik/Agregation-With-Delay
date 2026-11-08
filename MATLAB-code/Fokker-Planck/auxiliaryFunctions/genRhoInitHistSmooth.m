function [rhoInitHist] = genRhoInitHistSmooth(maxFlucts, rho, gridDims, stepDelay)

% maxFluc - maximal fluctuations added in one step
% rho - vectorized density
% gridDims - dimension of the grid used to discretize rho (without vectorization)
% stepDelay - step delay

N = numel(rho);

histDims = [gridDims, stepDelay];
rhoInitHist = reshape(makeFluctsSmooth(maxFlucts, histDims), [N, stepDelay]);
rhoInitHist = rhoInitHist + rho;

