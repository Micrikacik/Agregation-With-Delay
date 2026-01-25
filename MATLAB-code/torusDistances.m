function D = torusDistances(varargin)
% varargin - input either one matrix or two matrices
% Input matrix is Nxd, all d coordinates of all N agents in [0,1]^d
% D is NxN symmetric distance matrix on the unit torus

% Initialize the matrices
% 3rd matrix dimension corresponds to the d-dim space coordinate axes
% 1st and 2nd correspond to the agents
switch nargin
    case 1
        x_13(:,1,:) = varargin{1}(:,:);
        x_23(1,:,:) = varargin{1}(:,:);
    case 2
        x_13(:,1,:) = varargin{1}(:,:);
        x_23(1,:,:) = varargin{2}(:,:);
    otherwise
        error('Invalid number of input arguments.');
end


% Pairwise differences
x_diff = abs(x_13 - x_23);

%dx = abs(x(1,:).' - x(1,:));
%dy = abs(x(2,:).' - x(2,:));

% apply periodic boundary conditions
x_diff = min(x_diff, 1 - x_diff);

%dx = min(dx, 1 - dx);
%dy = min(dy, 1 - dy);

% Euclidean distance
D = sqrt(sum(x_diff.^2,3));

%D = sqrt(dx.^2 + dy.^2)
end