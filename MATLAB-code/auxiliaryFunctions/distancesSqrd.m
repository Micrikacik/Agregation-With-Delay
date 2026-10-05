function D = distancesSqrd(x_1,x_2)

% Calculates euclidian distances SQUARED between the positions 
% in 'x_1' and 'x_2', returning a distance matrix 'D'.
% This function is fast implementation of distances.m, but needs special
% output usage (returns SQUARED distances).
%
% INPUT:
%   x_1, x_2 (float matrices) - N by d matrices, each row represents
%       position vector.
%       Distances are calculated between vectors x_1(i,:) and x_2(j,:).
%
% OUTPUT:
%   D (nonnegative float matrix) - N by N distance matrix,
%       element D(i,j) is SQUARED distance ||x_1(i,:) - x_2(j,:)||^2

arguments
    x_1 (:,:) double
    x_2 (:,:) double = x_1 
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
        dx = x_1(:,1) - x_2(:,1).';
        
        % calculate Euclidean distance SQUARED
        D = dx.^2;
    case 2
        % pairwise differences
        dx = x_1(:,1) - x_2(:,1).';
        dy = x_1(:,2) - x_2(:,2).';
        
        % calculate Euclidean distance SQUARED
        D = dx.^2 + dy.^2;
    case 3
        % pairwise differences
        dx = x_1(:,1) - x_2(:,1).';
        dy = x_1(:,2) - x_2(:,2).';
        dz = x_1(:,3) - x_2(:,3).';
        
        % calculate Euclidean distance SQUARED
        D = dx.^2 + dy.^2 + dz.^2;
    otherwise
        error("Dimension %i not supported.", d)
end