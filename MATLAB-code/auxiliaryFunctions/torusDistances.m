function D = torusDistances(x_1, x_2, domainDims)

% Calculates euclidian distances on torus with dimensions 'domainDims' 
% (as a periodic hypercube) between the positions in 'x_1' and 'x_2',
% returning a distance matrix 'D'.
%
% INPUT:
%   x_1, x_2 (float matrices) - N by d matrices, each row represents
%       position vector.
%       Distances are calculated between vectors x_1(i,:) and x_2(j,:).
%   domainDims (float vector) - vector of length d (second dimension of 
%       'x_1, x_2'), values represent dimensions of the torus (as a 
%       periodic hypercube), i.e. 'domainDims(i)' is the length of the 
%       hypercube in the i-th dimension.
%       
% OUTPUT:
%   D (nonnegative float matrix) - N by N distance matrix.
%       Element D(i,j) is distance on torus with dimensions 'domainDims'.


arguments
    x_1 (:,:) double {mustBeFinite, mustBeNonempty, mustBeReal}
    x_2 (:,:) double {mustBeFinite, mustBeNonempty, mustBeReal} = x_1 
    domainDims (1,:) double {mustBeFinite, mustBeReal, mustBePositive} = ones(1, size(x_1, 2))

end

% Initialize the matrices
% 3rd matrix dimension corresponds to the d-dim space coordinate axes
% 1st and 2nd correspond to the agents
x_13(:,1,:) = x_1;
x_23(1,:,:) = x_2;

% Pairwise differences
x_diff = abs(x_13 - x_23);

% Transform the dimensions vector to use it for periodic distances
b_mat(1,1,:) = domainDims;

% Apply periodic boundary conditions
x_diff = min(x_diff, b_mat - x_diff);

% Euclidean distance
D = sqrt(sum(x_diff.^2, 3));

end