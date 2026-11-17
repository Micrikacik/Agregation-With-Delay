function stepDelay = getStepDelay(d, nonzero)

% Returns the default step delay for the given dimension 'd'.
%
% INPUT:
%   d (positive integer) - dimension.
%   nonzero (logical) - if true, instead of 0 this function returns 1.
%
% OUTPUT:
%   stepDelay (nonnegative integer) - default 'stepDelay'.

arguments
    d 
    nonzero = false
end

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

if nonzero && stepDelay == 0
    stepDelay = 1;
end