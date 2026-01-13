function [] = plotAgg(x)

N = size(x,1);
d = size(x,2);

%identify clusters:

%parameter setting for DBSCAN method
epsilon = 2.2/sqrt(2*N);
%epsilon = 1.3/sqrt(2*N);
minpts = 9;

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
    case 3 % Bad dbscan parameters
        max_idx = max(idx);
        hues = (idx + 1) ./ (max_idx + 1);
        colors = hsv2rgb([hues, ones(N,2)]);
        scatter3(x(:,1),x(:,2),x(:,3),[],colors);
        getframe;
end

end