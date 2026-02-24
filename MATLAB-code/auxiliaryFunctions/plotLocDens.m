function [] = plotLocDens(x, theta, dims)

arguments
    x
    theta
    dims = ones(size(x,2),1)
end

% Plots the agents and gives them color based on their local density theta
% Supports dimensions 1, 2, 3

N = size(x,1);
d = size(x,2);

%identify clusters:

%parameter setting for DBSCAN method
if d ~= 1 && d ~= 2 && d ~= 3
    error('Invalid number of dimensions. Only 1, 2, or 3 are supported.');
end

switch d
    case 1
        scatter(x(:,1),theta,[],redblue(theta),"filled");
        getframe;
    case 2
        scatter(x(:,1),x(:,2),[],redblue(theta),"filled");
        axis([0 dims(1) 0 dims(2)]);
        getframe;
    case 3
        scatter3(x(:,1),x(:,2),x(:,3),[],redblue(theta),'filled');
        axis([0 dims(1) 0 dims(2) 0 dims(3)]);
        getframe;
end

gradientTicks = 100;
hue = linspace(2/3,1,gradientTicks).';
colormap(hsv2rgb([hue, ones(gradientTicks,2)]))
colorbarTicks = 10;
cTicks = linspace(0,1,colorbarTicks);
thetas = linspace(min(theta),max(theta),colorbarTicks);
colorbar('Ticks',cTicks,'TickLabels',thetas)

    function result = redblue(theta)
        n = size(theta,1);
        minTheta = min(theta);
        maxTheta = max(theta);
        span = maxTheta - minTheta;
        hues = zeros(n,1);
        minHue = 2/3;
        maxHue = 1;
        for i = 1:n
            hues(i) = minHue + (theta(i) - minTheta) / span * (maxHue - minHue);
        end
        result = hsv2rgb([hues, ones(n,2)]);
    end

end