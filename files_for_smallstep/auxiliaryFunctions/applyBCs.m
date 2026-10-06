function [x] = applyBCs(x, boundConds, dims)

% Apply the boundary conditions on 'x'.
% 
% INPUT:
%   x (float matrix) - matrix of initial positions.
%       x(i,:) - position (float vector) in d-dim torus of the i-th agent.
%   boundConds (string) - type of the boundary conditions.
%       Must be one of the following strings:
%           "NoBoundary"
%           "Periodic"
%           "Reflective"
%   dims (positive float ROW vector) - dimensions of the simulation, i.e., 
%       dimensions of the box, in which the agents move.
%       This row vector must have its length the same as is the second
%       dimension of the matrix 'x'.
%
% OUTPUT:
%   x (float matrix) - the input 'x' on which BCs were applied.

arguments
    x (:,:) double
    boundConds (1,1) string
    dims (1,:) double {mustBePositive}
end

switch boundConds
    case "NoBoundary"
        % No BCs
    case "Periodic"
        % Periodic BCs
        x = mod(x ,dims);
    case "Reflective"
        % Reflective BCs (local)
        x = dims - abs(dims - abs(x));
    otherwise
        error('Invalid boundary conditions.');
end