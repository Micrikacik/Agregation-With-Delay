function D = distancesSqrd(x_1,x_2)

% Calculates euclidian distances SQUARED between the positions 
% in 'x_1' and 'x_2', returning a distance matrix 'D'.
% This function is fast implementation of distances.m, but needs special
% output usage (returns SQUARED distances).
%
% INPUT:
%   x_1, x_2 (float matrices) - N by d matrices, each row represents
%       position vector, 
%       distances are calculated between vectors x_1(i,:) and x_2(j,:)
%
% OUTPUT:
%   D (nonnegative float matrix) - N by N (symmetric) distance matrix,
%       element D(i,j) is SQUARED distance ||x_1(i,:) - x_2(j,:)||^2

arguments
    x_1 (:,:) float
    x_2 (:,:) float = x_1 
end

if size(x_2) ~= size(x_1)
    error("Wrong dimensions of input matrices")
end

d = size(x_1,2);

if size(x_1) ~= size(x_2)
    error("Input matrices have not equal sizes.")
end

switch d
    case 1
        % pairwise differences
        dx = abs(x_1(:,1) - x_2(:,1).');
        
        % calculate Euclidean distance
        D = dx.^2;
    case 2
        % pairwise differences
        dx = abs(x_1(:,1) - x_2(:,1).');
        dy = abs(x_1(:,2) - x_2(:,2).');
        
        % calculate Euclidean distance
        D = dx.^2 + dy.^2;
    case 3
        % pairwise differences
        dx = abs(x_1(:,1) - x_2(:,1).');
        dy = abs(x_1(:,2) - x_2(:,2).');
        dz = abs(x_1(:,3) - x_2(:,3).');
        
        % calculate Euclidean distance
        D = dx.^2 + dy.^2 + dz.^2;
    otherwise
        error("Dimension %i not supported.", d)

% % Initialize the matrices
% % 3rd matrix dimension corresponds to the d-dim space coordinate axes
% % 1st and 2nd correspond to the agents
% x_13(:,1,:) = x_1;
% x_23(1,:,:) = x_2;
% 
% % Pairwise differences
% x_diff = abs(x_13 - x_23);
% 
% % Euclidean distance
% D = (x_diff(:,:,1).^2 + x_diff(:,:,2).^2);

end