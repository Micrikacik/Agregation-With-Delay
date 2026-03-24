function D = distancesSqrd(x_1,x_2)
% x_1, x_2 are Nxd, all d coordinates of all N agents in [0,1]^d
% D will be NxN (symmetric) SQUARED-distance matrix on the torus where rowrs
% correspond to x_1 and columns correspond to x_2

arguments
    x_1
    x_2 = x_1
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