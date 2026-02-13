function D = torusDistances(x_1,x_2,dims)
% x_1, x_2 are Nxd, all d coordinates of all N agents in [0,1]^d
% D will be NxN symmetric distance matrix on the torus with dimensions
% defined by dims

arguments
    x_1 
    x_2 = x_1
    dims = ones(size(x_1,2),1)
end

% Initialize the matrices
% 3rd matrix dimension corresponds to the d-dim space coordinate axes
% 1st and 2nd correspond to the agents
x_13(:,1,:) = x_1;
x_23(1,:,:) = x_2;

% Pairwise differences
x_diff = abs(x_13 - x_23);

% Transform the dimensions vector to use it for periodic distances
b_mat(1,1,:) = dims;

% Apply periodic boundary conditions
x_diff = min(x_diff, b_mat - x_diff);

% Euclidean distance
D = sqrt(sum(x_diff.^2,3));

end