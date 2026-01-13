% Calculation of mutual distances over all pairs of points in 2D
% on a unit torus (unit rectangle with periodic BC)
%
% Input:
% x = array(2,N) - point coordinates in [0,1]x[0,1]
%
% Output:
% y = array(N,N) - symmetric distance matrix on the unit torus

function D = torusDistances(x)

% pairwise differences
dx = abs(x(1,:).' - x(1,:));
dy = abs(x(2,:).' - x(2,:));

% apply periodic boundary conditions
dx = min(dx, 1 - dx);
dy = min(dy, 1 - dy);

% calculate Euclidean distance
D = sqrt(dx.^2 + dy.^2);

end
