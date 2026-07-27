function [] = plotMCData(outlierCountsStat, clusterCountsStat, clusterSizesStat, clusterFreePercentage, plotSettings)

arguments
    outlierCountsStat (1,1) struct
    clusterCountsStat (1,1) struct
    clusterSizesStat (1,1) struct
    clusterFreePercentage (:,1)
    plotSettings (1,1) struct = struct()
end

% This function visualizes statistical measures returned by processMCData.m.
% TODOOOOOOOOOOOOOOOOOOOOOOOOOO
% It can be specified, which of them will be plotted in the input struct
% 'plotSettings' with fields (of logical values): 
%   plotSettings.minimaMaxima
%   plotSettings.histograms
%   plotSettings.means
%   plotSettings.medians
%   plotSettings.variances
% plotSettings.xLabel
% plotSettings.xTickLabels
% plotSettings.fontSize

% Set default plotSettings values
if ~isfield(plotSettings, 'minimaMaxima')
    plotSettings.minimaMaxima = true;
else
    validateattributes(plotSettings.minimaMaxima, {'logical'}, {'scalar'})
end

if ~isfield(plotSettings, 'histograms')
    plotSettings.histograms = true;
else
    validateattributes(plotSettings.histograms, {'logical'}, {'scalar'})
end

if ~isfield(plotSettings, 'means')
    plotSettings.means = true;
else
    validateattributes(plotSettings.means, {'logical'}, {'scalar'})
end

if ~isfield(plotSettings, 'medians')
    plotSettings.medians = true;
else
    validateattributes(plotSettings.medians, {'logical'}, {'scalar'})
end

if ~isfield(plotSettings, 'variances')
    plotSettings.variances = true;
else
    validateattributes(plotSettings.variances, {'logical'}, {'scalar'})
end

if ~isfield(plotSettings, "xLabel")
    plotSettings.xLabel = "Delay \tau";
else
    validateattributes(plotSettings.xLabel, {'string'}, {'scalar'})
end

if ~isfield(plotSettings, "xTickLabels")
    plotSettings.xTickLabels = 0:0.03:0.42;
else
    validateattributes(plotSettings.xTickLabels, {'string'}, {'vector'})
end

if ~isfield(plotSettings, "fontSize")
    plotSettings.fontSize = getFontSize();
else
    validateattributes(plotSettings.fontSize, {'double'}, {'vector', 'positive'})
end


figure
MCCount = length(clusterFreePercentage);

% Settings of the graphs
xPoints = (1:MCCount).';
minMaxPlot.color = 'red';
minMaxPlot.lineStyle = 'none';
minMaxPlot.lineWidth = 3;
meanPlot.color = 'red';
meanPlot.lineStyle = '-.';
meanPlot.lineWidth = 3;
medianPlot.lineStyle = ':';
medianPlot.color = 'magenta';
medianPlot.lineWidth = 3;
varPlot.color = 'black';
varPlot.lineStyle = 'none';
varPlot.marker = 'x';
varPlot.markerSize = 7;
varPlot.lineWidth = 2;
histPlot.valOffset = 0.2;
histPlot.color = [0.2,0.3,0.9];
percMargin = 0.05;

% Plot number of outliers
subplot(2,2,1);
hold on

plotStatData(outlierCountsStat)
setAxis(outlierCountsStat)

hold off
% Create ylabel
ylabel({'Number of outliers'});
% Create xlabel
xlabel({plotSettings.xLabel});
% Set xticks & xticklabels
xticks(1:length(plotSettings.xTickLabels))
xticklabels(plotSettings.xTickLabels)
% Create title
title({'Number of outliers'});


% Plot number of clusters
subplot(2,2,2);
hold on

plotStatData(clusterCountsStat)
setAxis(clusterCountsStat)

hold off
% Create ylabel
ylabel({'number of clusters'});
% Create xlabel
xlabel({plotSettings.xLabel});
% Set xticks & xticklabels
xticks(1:length(plotSettings.xTickLabels))
xticklabels(plotSettings.xTickLabels)
% Create title
title({'Number of clusters'});


% Plot cluster sizes
subplot(2,2,3);
hold on

plotStatData(clusterSizesStat)
setAxis(clusterSizesStat)

hold off
% Create ylabel
ylabel({'cluster size'});
% Create xlabel
xlabel({plotSettings.xLabel});
% Set xticks & xticklabels
xticks(1:length(plotSettings.xTickLabels))
xticklabels(plotSettings.xTickLabels)
% Create title
title({'Cluster sizes'});


%plot percentage of cluster-free outcomes
subplot(2,2,4);
plot(xPoints, clusterFreePercentage, 'o');
axis([0 MCCount+1 -percMargin*100 100+percMargin*100]);
% Create ylabel
ylabel({'percentage'});
% Create xlabel
xlabel(plotSettings.xLabel);
% Set xticks & xticklabels
xticks(1:length(plotSettings.xTickLabels))
xticklabels(plotSettings.xTickLabels)
% Create title
title({'% of cluster-free outcomes'});

% Font size
fontsize(plotSettings.fontSize, "points");

function [] = setAxis(data)

    upper = 0;
    lower = 0;

    if plotSettings.minimaMaxima
        upper = max(upper, max(data.maxima));
        lower = min(lower, min(data.minima));
    end

    if plotSettings.histograms
        upper = max(upper, max(data.maxima));
        lower = min(lower, min(data.minima));
    end

    if plotSettings.means
        upper = max(upper, max(data.means));
        lower = min(lower, min(data.means));
    end

    if plotSettings.medians
        upper = max(upper, max(data.medians));
        lower = min(lower, min(data.medians));
    end

    if plotSettings.variances
        upper = max(upper, max(data.means + sqrt(data.variances)));
        lower = min(lower, min(data.means - sqrt(data.variances)));
    end

    diff = upper - lower;
    upper = upper + percMargin * diff;
    lower = lower - percMargin * diff;

    axis([0, MCCount+1, lower, upper])

end

function [] = plotStatData(data)

    if plotSettings.minimaMaxima
        % Plot minima and maxima
        errorbar(xPoints, ...
            data.means, ...
            data.means - data.minima, ...
            - data.means + data.maxima, ...
            'LineWidth',minMaxPlot.lineWidth,'LineStyle',minMaxPlot.lineStyle,'Color',minMaxPlot.color);
    end

    if plotSettings.histograms
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

    if plotSettings.means
        % Plot means connected with line
        plot(xPoints,data.means,meanPlot)
    end

    if plotSettings.medians
        % Plot medians connected with lines
        plot(xPoints,data.medians,medianPlot)
    end

    if plotSettings.variances
        % Plot variances around means
        plot(xPoints,[ ...
            data.means + sqrt(data.variances); ...
            data.means - sqrt(data.variances) ...
            ].', ...
            varPlot)
    end
end

end