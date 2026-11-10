function stepPlotMod = getStepPlotMod(d)

% Returns the default step plot mod for the given dimension 'd'.
%
% INPUT:
%   d (positive integer) - dimension.
%
% OUTPUT:
%   stepPlotMod (positive integer or -1) - default 'stepPlotMod'.

switch d
    case 1
        stepPlotMod = 2;  
    case 2
        stepPlotMod = 4;  
    case 3
        stepPlotMod = 8;  
    otherwise
        stepPlotMod = -1;
end