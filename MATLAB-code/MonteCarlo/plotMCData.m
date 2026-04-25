function [] = plotMCData(outlierCountsStat,clusterCountsStat,clusterSizesStat,clusterFreePercentage,plotData,N)

arguments
    outlierCountsStat 
    clusterCountsStat 
    clusterSizesStat 
    clusterFreePercentage 
    plotData = struct()
    N = 400
end

% This function visualizes statistical measures returned by processMCData.m.
% It can be specified, which of them will be plotted in the input struct
% 'plotData' with fields (of logical values): 
%   plotData.minimaMaxima
%   plotData.histograms
%   plotData.means
%   plotData.medians
%   plotData.variances

% Set default plotData values
if ~isfield(plotData,'minimaMaxima')
    plotData.minimaMaxima = true;
else
    validateattributes(plotData.minimaMaxima,{'logical'},{'scalar'})
end

if ~isfield(plotData,'histograms')
    plotData.histograms = false;
else
    validateattributes(plotData.histograms,{'logical'},{'scalar'})
end

if ~isfield(plotData,'means')
    plotData.means = true;
else
    validateattributes(plotData.means,{'logical'},{'scalar'})
end

if ~isfield(plotData,'medians')
    plotData.medians = true;
else
    validateattributes(plotData.medians,{'logical'},{'scalar'})
end

if ~isfield(plotData,'variances')
    plotData.variances = false;
else
    validateattributes(plotData.variances,{'logical'},{'scalar'})
end


figure
MCCount = length(clusterFreePercentage);

% Settings of the graphs
xPoints = (1:MCCount).';
minMaxPlot.color = 'red';
minMaxPlot.lineStyle = 'none';
minMaxPlot.lineWidth = 2;
meanPlot.color = minMaxPlot.color;
meanPlot.lineStyle = '--';
meanPlot.marker = 'diamond';
meanPlot.lineWidth = 2;
medianPlot.lineStyle = ':';
medianPlot.color = 'magenta';
medianPlot.lineWidth = 2;
varPlot.color = 'yellow';
varPlot.lineStyle = 'none';
varPlot.marker = 'x';
varPlot.lineWidth = 2;
histPlot.valOffset = 0.2;
histPlot.color = [0.2,0.3,0.9];


% Plot number of outliers
subplot(2,2,1);
hold on

plotStatData(outlierCountsStat)

hold off
axis([1, MCCount+1, -1, N]);
% Create ylabel
ylabel({'number of outliers'});
% Create xlabel
xlabel({'stepDelay'});
% Create title
title({'min, avg and max number of outliers'});


% Plot number of clusters
subplot(2,2,2);
hold on

plotStatData(clusterCountsStat)

hold off
axis([1, MCCount+1, -1, 16]);
% Create ylabel
ylabel({'number of clusters'});
% Create xlabel
xlabel({'stepDelay'});
% Create title
title({'min, avg and max number of clusters'});


% Plot cluster sizes
subplot(2,2,3);
hold on

plotStatData(clusterSizesStat)

hold off
axis([1, MCCount+1, -1, N]);
% Create ylabel
ylabel({'cluster size'});
% Create xlabel
xlabel({'stepDelay'});
% Create title
title({'min, avg and max cluster size'});


%plot percentage of cluster-free outcomes
subplot(2,2,4);
plot(xPoints, clusterFreePercentage, 'o');
axis([1 MCCount -5 105]);
% Create ylabel
ylabel({'percentage'});
% Create xlabel
xlabel({'stepDelay'});
% Create title
title({'% of cluster-free outcomes'});

function [] = plotStatData(data)

    if plotData.minimaMaxima
        % Plot minima and maxima
        errorbar(xPoints, ...
            data.means, ...
            data.means - data.minima, ...
            - data.means + data.maxima, ...
            'LineWidth',minMaxPlot.lineWidth,'LineStyle',minMaxPlot.lineStyle,'Color',minMaxPlot.color);
    end

    if plotData.histograms
        % Plot horizontal histograms
        for i = 1:length(data.histograms)
            histMax = max(data.histograms{i}) / (1 - histPlot.valOffset);
            for j = 1:length(data.histograms{i})
                if data.histograms{i}(j) == 0
                    continue
                end
                rectangle('Position',[ ...
                    i, ...
                    data.minima(i) + j - 1.5, ...
                    data.histograms{i}(j) / histMax, ...
                    1, ...
                    ], ...
                    'FaceColor',histPlot.color,'EdgeColor',histPlot.color)
            end
        end
    end

    if plotData.means
        % Plot means connected with line
        plot(xPoints,data.means,meanPlot)
    end

    if plotData.medians
        % Plot medians connected with lines
        plot(xPoints,data.medians,medianPlot)
    end

    if plotData.variances
        % Plot variances around means
        plot(xPoints,[ ...
            data.means + sqrt(data.variances); ...
            data.means - sqrt(data.variances) ...
            ].', ...
            varPlot)
    end
end

end