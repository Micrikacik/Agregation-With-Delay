function D = torusDistancesSqrd(x_1, x_2, dims)

% Calculates euclidian distances SQUARED on torus with dimensions 'dims'
% (as a periodic hypercube) between the positions in 'x_1' and 'x_2', 
% returning a distance matrix 'D'.
% This function is fast implementation of torusDistances.m, but needs 
% special output usage (returns SQUARED distances).
%
% INPUT:
%   x_1, x_2 (float matrices) - N by d matrices, each row represents
%       position vector.
%       Distances are calculated between vectors x_1(i,:) and x_2(j,:).
%   dims (float vector) - vector of length d (second dimension of 
%       'x_1, x_2'), values represent dimensions of the torus (as a periodic 
%       hypercube), i.e. 'dims(i)' is the length of the hypercube in the i-th
%       dimension.
%       
% OUTPUT:
%   D (nonnegative float matrix) - N by N (symmetric) distance matrix,
%       element D(i,j) is SQUARED distance on torus with dimensions 'dims'
%       "||x_1(i,:) - x_2(j,:)||^2"

arguments
    x_1 (:,:) double
    x_2 (:,:) double = x_1 
    dims = ones(1, size(x_1, 2))
end

if size(x_1) ~= size(x_2)
    error("Input matrices do not have equal sizes.")
end

d = size(x_1, 2);

if d > 3
    error("Dimension %i not supported.", d)
end

switch d
    case 1
        % pairwise differences
        dx = abs(x_1(:,1) - x_2(:,1).');
        
        % apply periodic boundary conditions
        dx = min(dx, dims(1) - dx);
        
        % calculate Euclidean distance SQUARED
        D = dx.^2;
    case 2
        % pairwise differences
        dx = abs(x_1(:,1) - x_2(:,1).');
        dy = abs(x_1(:,2) - x_2(:,2).');
        
        % apply periodic boundary conditions
        dx = min(dx, dims(1) - dx);
        dy = min(dy, dims(2) - dy);
        
        % calculate Euclidean distance SQUARED
        D = dx.^2 + dy.^2;
    case 3
        % pairwise differences
        dx = abs(x_1(:,1) - x_2(:,1).');
        dy = abs(x_1(:,2) - x_2(:,2).');
        dz = abs(x_1(:,3) - x_2(:,3).');
        
        % apply periodic boundary conditions
        dx = min(dx, dims(1) - dx);
        dy = min(dy, dims(2) - dy);
        dz = min(dz, dims(3) - dz);
        
        % calculate Euclidean distance SQUARED
        D = dx.^2 + dy.^2 + dz.^2;
    otherwise
        error("Dimension %i not supported.", d)
end