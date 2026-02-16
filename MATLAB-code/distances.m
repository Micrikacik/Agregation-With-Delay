function D = distances(x_1,x_2)
% x_1, x_2 are Nxd, all d coordinates of all N agents in [0,1]^d
% D will be NxN (symmetric) distance matrix

arguments
    x_1 
    x_2 = x_1
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