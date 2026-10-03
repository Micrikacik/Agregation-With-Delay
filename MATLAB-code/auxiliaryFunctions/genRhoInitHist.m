function [rhoInitHist] = genRhoInitHist(maxInitFluc, gridPointCount, stepDelay, volume, d)

rhoInitHist = zeros(gridPointCount, stepDelay);

for t=1:stepDelay
    tempFluc = makeFluctuations(maxInitFluc, gridPointCount, d);
    tempRho = 1 / volume + tempFluc;
    rhoInitHist(:,t) = tempRho;
end