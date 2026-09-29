function [fluctuations] = makeFluctuations(maxFluc, gridPointCount, d)

switch d
    case 1
        fluctuations = maxFluc * 2 * (rand(gridPointCount, 1) - 0.5);
    otherwise
        error('Dimension d = %.i not implemented.', d);
end