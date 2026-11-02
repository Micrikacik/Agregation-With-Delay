function [rhoInitHist] = genRhoInitHist(maxFlucts, rho0, stepDelay)

% maxFluc - maximal fluctuations added in one step
% rho0 - vectorized density
% stepDelay - step delay

N = numel(rho0);

rhoInitHist = zeros(N, stepDelay);

for t=1:stepDelay
    tempFluc = makeFlucts(maxFlucts, N);
    tempRho = rho0 + tempFluc;
    rhoInitHist(:,t) = tempRho;
end