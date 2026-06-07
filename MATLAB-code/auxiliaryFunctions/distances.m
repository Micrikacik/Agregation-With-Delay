function D = distances(x_1,x_2)

% Calculates euclidian distances between the positions in 'x_1' and 'x_2',
% returning a distance matrix 'D'.
%
% INPUT:
%   x_1, x_2 (float matrices) - N by d matrices, each row represents
%       position vector, 
%       distances are calculated between vectors x_1(i,:) and x_2(j,:)
%
% OUTPUT:
%   D (nonnegative float matrix) - N by N (symmetric) distance matrix,
%       element D(i,j) is distance ||x_1(i,:) - x_2(j,:)||

arguments
    x_1 (:,:) double
    x_2 (:,:) double = x_1 
end

if size(x_2) ~= size(x_1)
    error("Wrong dimensions of input matrices")
end

% Initialize the matrices
% 3rd matrix dimension corresponds to the d-dim space coordinate axes
% 1st and 2nd correspond to the agents
x_13(:,1,:) = x_1;
x_23(1,:,:) = x_2;

% Pairwise differences
x_diff = x_13 - x_23;

% Euclidean distance
D = sqrt(sum(x_diff.^2,3));

end