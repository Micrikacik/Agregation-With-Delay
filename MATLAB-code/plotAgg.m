function [] = plotAgg(x)

% Plots the agents and highligts agregation
% Supports dimensions 1, 2, 3

N = size(x,1);
d = size(x,2);

%identify clusters:

%parameter setting for DBSCAN method
switch d
    case 1
        epsilon = 0.8/sqrt(2*N);
        minpts = 14;
    case 2
        epsilon = 1.3/sqrt(2*N);
        minpts = 12;
    case 3
        epsilon = 2.3/sqrt(2*N);
        minpts = 9;
    otherwise
        error('Invalid number of dimensions. Only 1, 2, or 3 are supported.');
end

%calculate distance over the torus
dist = torusDistances(x);
        
%identify clusters
idx = dbscan(dist,epsilon,minpts,'Distance','precomputed');

switch d
    case 1
        gscatter(x(:,1),ones(N,1).*0.5,idx);
        getframe;
    case 2
        gscatter(x(:,1),x(:,2),idx);
        getframe;
    case 3
        max_idx = max(idx);
        hues = (idx + 1) ./ (max_idx + 1);
        colors = hsv2rgb([hues, ones(N,2)]);
        scatter3(x(:,1),x(:,2),x(:,3),[],colors,'filled');
        getframe;
end

end