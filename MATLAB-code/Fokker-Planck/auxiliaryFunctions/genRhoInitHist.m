function [rhoInitHist] = genRhoInitHist(maxFluc, rho0, stepDelay)

% maxFluc - maximal fluctuations added in one step
% rho0 - vectorized density
% stepDelay - step delay

V = numel(rho0);

rhoInitHist = zeros(V, stepDelay);

tempFluc = makeFluctuations(maxFluc, V);
tempRho = rho0 + tempFluc;
rhoInitHist(:,1) = tempRho;

for t=2:stepDelay
    tempFluc = makeFluctuations(maxFluc, V);
    tempRho = rhoInitHist(:,t-1) + tempFluc;
    rhoInitHist(:,t) = tempRho;
end