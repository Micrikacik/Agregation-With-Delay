function [rhoInitHist] = genRhoInitHist(maxInitFluc, gridPointCount, stepDelay, L, d)

rhoInitHist = zeros(gridPointCount, stepDelay);

for t=1:stepDelay
    tempFluc = makeFluctuations(maxInitFluc, gridPointCount, d);
    tempRho = 1 / L + tempFluc;
    rhoInitHist(:,t) = tempRho;
end