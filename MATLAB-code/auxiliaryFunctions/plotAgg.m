function plotAgg(x, theta, boundConds, domainDims)

% Plots the agents in a new figure and highligts agregation.
%
% INPUT:
%   x (float matrix) - matrix of positions.
%       x(i,:) - position (float vector) in d-dim of the i-th agent.
%   theta (float nonnegative vector) - vector of densities 'theta'
%       corresponding to the agents with positions in rows of 'x'.
%   boundConds (string) - type of the boundary conditions.
%       Must be one of the following strings:
%           "NoBoundary"
%           "Periodic"
%           "Reflective"
%   domainDims (positive float vector) - dimensions of the simulation, 
%       i.e., dimensions of the box, in which the agents move.
%       This row vector must have its length the same as is the second
%       dimension of the matrix 'x'.

arguments
    x (:,:) double {mustBeFinite, mustBeNonempty, mustBeReal}
    theta (:,1) double {mustBeFinite, mustBeNonempty, mustBeReal, mustBeNonnegative} ...
        = zeros(size(x,1),1)
    boundConds (1,1) string = ones(size(x,2),1)
    domainDims (1,:) double {mustBeFinite, mustBeReal, mustBePositive} ...
        = "Periodic"
end

N = size(x,1);
d = size(x,2);

if d > 3
    error("Dimension %i not supported for plots.", d)
end
if ~isequal(size(theta), [N,1])
    error("Input 'theta' has wrong size, which does not correspond to 'x'.")
end
if d ~= length(domainDims)
    error("Input 'domainDims' has wrong size, which does not correspond to 'x'.")
end

volume = prod(domainDims);

%identify clusters:

switch boundConds
    case "Periodic"
        %calculate distances over the torus
        dists = torusDistances(x, x, domainDims);
    case "Reflective"
        %calculate distances
        dists = distances(x);
    otherwise
        dists = distances(x);
end

epsilon = getIntRad(d);
minpts = getMinClusterSize(N, d, volume);
        
%identify clusters
idx = dbscan(dists, epsilon, minpts, 'Distance', 'precomputed');

max_idx_or_0 = max(max(idx), 0);
out_idx = idx <= 0;
clus_idx = idx > 0;

maxDistColors = maxdistcolor(max_idx_or_0 + 1);
colors = idxMap(maxDistColors(1:end-1,:), maxDistColors(end,:), idx);

figure

switch d
    case 1
        hold on
        scatter(x(clus_idx,1), theta(clus_idx), getPointSize()^2, colors(clus_idx,:), ".")
        scatter(x(out_idx,1), theta(out_idx), getPointSize()^2, colors(out_idx,:), "x", 'LineWidth', 3)
        hold off
        axis([0 domainDims(1) 0 (max(theta)*getMaxPlotMultip())]);
    case 2
        hold on
        scatter(x(clus_idx,1), x(clus_idx,2), getPointSize()^2, colors(clus_idx,:), ".")
        scatter(x(out_idx,1), x(out_idx,2), getPointSize()^2, colors(out_idx,:), "x", 'LineWidth', 3)
        hold off
        axis([0 domainDims(1) 0 domainDims(2)]);
    case 3
        hold on
        scatter3(x(clus_idx,1), x(clus_idx,2), x(clus_idx,3), getPointSize()^2, colors(clus_idx,:), ".")
        scatter3(x(out_idx,1), x(out_idx,2), x(out_idx,3), getPointSize()^2, colors(out_idx,:), "x", 'LineWidth', 3)
        hold off
        axis([0 domainDims(1) 0 domainDims(2) 0 domainDims(3)]);
end

fontsize(getFontSize(),'points')
theme(gcf, 'light')

    function result = idxMap(code,default,idx)
        result = zeros(length(idx), size(code, 2));
        for i = 1:length(idx)
            setDefault = true;
            for c = 1:size(code, 1)
                if idx(i) == c
                    result(i,:) = code(c,:);
                    setDefault = false;
                    continue
                end
            end
            if setDefault
                result(i,:) = default;
            end
        end
    end

end