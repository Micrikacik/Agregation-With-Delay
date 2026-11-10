function stepDelay = getStepDelay(d)

% Returns the default step delay for the given dimension 'd'.
%
% INPUT:
%   d (positive integer) - dimension.
%
% OUTPUT:
%   stepDelay (nonnegative integer) - default 'stepDelay'.

switch d
    case 1 
        stepDelay = 120;
    case 2 
        stepDelay = 60;
    case 3 
        stepDelay = 40;
    otherwise 
        stepDelay = 0;
end