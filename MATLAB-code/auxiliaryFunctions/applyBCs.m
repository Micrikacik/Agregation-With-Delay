function [x] = applyBCs(x, boundConds, domainDims)

% Apply the boundary conditions on 'x'.
% 
% INPUT:
%   x (float matrix) - matrix of positions.
%       x(i,:) - position (float vector) in d-dim of the i-th agent.
%   boundConds (string) - type of the boundary conditions.
%       Must be one of the following strings:
%           "NoBoundary"
%           "Periodic"
%           "Reflective"
%   domainDims (positive float vector) - dimensions of the simulation, 
%       i.e., dimensions of the box, in which the agents move.
%       This row vector must have its length the same as is the second
%       dimension of the matrix 'x'.
%
% OUTPUT:
%   x (float matrix) - the input 'x' on which BCs were applied.

arguments
    x (:,:) double {mustBeFinite, mustBeNonempty, mustBeReal}
    boundConds (1,1) string
    domainDims (1,:) double {mustBeFinite, mustBeReal, mustBePositive}
end

if size(x, 2) ~= length(domainDims)
    error("Input 'domainDims' has wrong size, which does not correspond to 'x'.")
end

switch boundConds
    case "NoBoundary"
        % No BCs
    case "Periodic"
        % Periodic BCs
        x = mod(x, domainDims);
    case "Reflective"
        % Reflective BCs (local)
        x = domainDims - abs(domainDims - abs(x));
    otherwise
        error('Invalid boundary conditions.');
end