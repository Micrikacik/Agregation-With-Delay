function [] = plotAgg(x, theta, dims, boundConds)

arguments
    x
    theta = zeros(size(x,1),1)
    dims = ones(size(x,2),1)
    boundConds = "Periodic"
end

% Plots the agents and highligts agregation
% Supports dimensions 1, 2, 3

N = size(x,1);
d = size(x,2);
volume = prod(dims);

%identify clusters:

%parameter setting for DBSCAN method
switch d
    case 1
        %epsilon = 0.7/20;
        %minpts = 14;
    case 2
        %epsilon = 1/20;
        %minpts = 12;
    case 3
        %epsilon = 1.8/20;
        %minpts = 9;
    otherwise
        error('Invalid number of dimensions. Only 1, 2, or 3 are supported.');
end

switch boundConds
    case "Periodic"
        %calculate distances over the torus
        dists = torusDistances(x,x,dims);
    case "Reflective"
        %calculate distances
        dists = distances(x);
    otherwise
        dists = distances(x);
end

epsilon = getIntRad(d);
minpts = getMinClusterSize(N, volume);
        
%identify clusters
idx = dbscan(dists, epsilon, minpts, 'Distance', 'precomputed');

max_idx = max(idx);
out_idx = idx <= 0;
clus_idx = idx > 0;

maxDistColors = maxdistcolor(max_idx + 1);
colors = idxMap(maxDistColors(1:end-1,:), maxDistColors(end,:), idx);

figure

switch d
    case 1
        hold on
        scatter(x(clus_idx,1), theta(clus_idx), getPointSize()^2, colors(clus_idx,:), ".")
        scatter(x(out_idx,1), theta(out_idx), getPointSize()^2, colors(out_idx,:), "x", 'LineWidth', 3)
        hold off
        axis([0 dims(1) 0 1]);
    case 2
        hold on
        scatter(x(clus_idx,1), x(clus_idx,2), getPointSize()^2, colors(clus_idx,:), ".")
        scatter(x(out_idx,1), x(out_idx,2), getPointSize()^2, colors(out_idx,:), "x", 'LineWidth', 3)
        hold off
        axis([0 dims(1) 0 dims(2)]);
    case 3
        hold on
        scatter3(x(clus_idx,1), x(clus_idx,2), x(clus_idx,3), getPointSize()^2, colors(clus_idx,:), ".")
        scatter3(x(out_idx,1), x(out_idx,2), x(out_idx,3), getPointSize()^2, colors(out_idx,:), "x", 'LineWidth', 3)
        hold off
        axis([0 dims(1) 0 dims(2) 0 dims(3)]);
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