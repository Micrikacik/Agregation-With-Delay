function xInitHist = genInitHist(x, dt, stepDelay, boundConds, domainDims)

% Generates initial history of 'blind motion'.
%
% INPUT:
%   x (float matrix) - matrix of initial positions.
%       x(i,:) - position (float vector) in d-dim torus of the i-th agent.
%   dt (positive float) - time step length.
%   stepDelay (nonnegative integer) - number of steps used to delay the simulation.
%       Must be > 0 for this function.
%   boundConds (string) - type of the boundary conditions.
%       Must be one of the following strings:
%           "NoBoundary"
%           "Periodic"
%           "Reflective"
%   domainDims (positive float ROW vector) - dimensions of the simulation, i.e., 
%       dimensions of the box, in which the agents move.
%       This row vector must have its length the same as is the second
%       dimension of the matrix 'x'.

arguments
    x (:,:) double
    dt (1,1) double {mustBePositive}
    stepDelay (1,1) double {mustBePositive, mustBeInteger}
    boundConds (1,1) string
    domainDims (1,:) double {mustBePositive}
end

[N, d] = size(x);

xInitHist = zeros([N, d, stepDelay]);

if d ~= length(domainDims)
    error("Input 'domainDims' has wrong size, which does not correspond to 'x'.")
end

xInitHist(:,:,1) = x - sqrt(dt) * randn(N, d);    

% Apply BCs
xInitHist(:,:,1) = applyBCs(xInitHist(:,:,1), boundConds, domainDims);

for i = 2:stepDelay
    xInitHist(:,:,i) = xInitHist(:,:,i-1) - sqrt(dt) * randn(N, d);   

    % Apply BCs
    xInitHist(:,:,i) = applyBCs(xInitHist(:,:,i), boundConds, domainDims);
end
