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

        satDef = 0;
        valDef = 0.3;
        satSwap = 0.5;
        valSwap = 0.5;
    otherwise
        error('Invalid number of dimensions. Only 1, 2, or 3 are supported.');
end

kappa_1 = 2;  % "volume" of a unit 1-ball
kappa = pi^(d/2) / gamma(d / 2 + 1);  % volume of a unit d-ball
intRad_1 = 0.05^2 / 2 * pi;
intRad = (kappa_1 / kappa * intRad_1)^(1/d);         % interaction radius
epsilon = intRad;
minpts = 12 * (N / volume / 400);

switch boundConds
    case "Periodic"
        %calculate distances over the torus
        dist = torusDistances(x,x,dims);
    case "Ricochet"
        %calculate distances
        dist = distances(x);
    otherwise
        dist = distances(x);
end

        
%identify clusters
idx = dbscan(dist,epsilon,minpts,'Distance','precomputed');

switch d
    case 1
        gscatter(x(:,1),theta,idx);
        axis([0 dims(1) 0 1]);
        getframe;
    case 2
        gscatter(x(:,1),x(:,2),idx);
        axis([0 dims(1) 0 dims(2)]);
        getframe;
    case 3
        max_idx = max(idx);
        hues = (idx + 1) ./ (max_idx + 1);
        satsCode = repmat([satSwap;1],[round(max_idx/2),1]);
        valsCode = repmat([1;valSwap],[round(max_idx/2),1]);
        if mod(max_idx,2) == 1
            satsCode = [satsCode;satSwap];
            valsCode = [valsCode;1];
        end
        sats = idxMap(satsCode,satDef,idx);
        vals = idxMap(valsCode,valDef,idx);
        colors = hsv2rgb([hues, sats, vals]);
        scatter3(x(:,1),x(:,2),x(:,3),[],colors,'filled');
        axis([0 dims(1) 0 dims(2) 0 dims(3)]);
        getframe;
end

    function result = idxMap(code,default,idx)
        result = zeros(length(idx),1);
        for i = 1:length(idx)
            setDefault = true;
            for c = 1:length(code)
                if idx(i) == c
                    result(i) = code(c);
                    setDefault = false;
                    continue
                end
            end
            if setDefault
                result(i) = default;
            end
        end
    end

end