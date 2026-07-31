function [] = plotMCData(outlierCountsStat, clusterCountsStat, clusterSizesStat, clusterFreePercentages, plotSettings)

arguments
    outlierCountsStat (1,1) struct
    clusterCountsStat (1,1) struct
    clusterSizesStat (1,1) struct
    clusterFreePercentages (:,1)
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
%   plotSettings.deviations
% plotSettings.xLabel
% plotSettings.xTickLabels
% plotSettings.fontSize
% plotSettings.separateFigs

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

if ~isfield(plotSettings, 'deviations')
    plotSettings.deviations = true;
else
    validateattributes(plotSettings.deviations, {'logical'}, {'scalar'})
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

if ~isfield(plotSettings, "separateFigs")
    plotSettings.separateFigs = false;
else
    validateattributes(plotSettings.separateFigs, {'logical'}, {'scalar'})
end

MCCount = length(clusterFreePercentages);

% Settings of the graphs
xPoints = (1:MCCount).';
minMaxPlot.color = [0.85,0.6,0];
minMaxPlot.lineStyle = 'none';
minMaxPlot.lineWidth = 2;
meanPlot.color = 'red';
meanPlot.lineStyle = ':';
meanPlot.lineWidth = 2;
meanPlot.marker = '^';
meanPlot.markerSize = 8;
meanPlot.markerFaceColor = meanPlot.color;
medianPlot.lineStyle = ':';
medianPlot.color = 'blue';
medianPlot.lineWidth = 2;
medianPlot.marker = 'v';
medianPlot.markerSize = 8;
medianPlot.markerFaceColor = medianPlot.color;
devPlot.color = 'black';
devPlot.lineStyle = 'none';
devPlot.marker = 'x';
devPlot.markerSize = 15;
devPlot.lineWidth = 2;
histPlot.valOffset = 0.1;
histPlot.color = [0.3,0.8,0];
clustFreePercPlot.lineStyle = ':';
clustFreePercPlot.color = 'black';
clustFreePercPlot.lineWidth = 2;
clustFreePercPlot.marker = '.';
clustFreePercPlot.markerSize = 30;
percMargin = 0.05;

% Make new figure
figure

% Plot number of outliers
if plotSettings.separateFigs
    % Figure already exists
else
    subplot(2,2,1);
end

plotStatData(outlierCountsStat)
setAxis(outlierCountsStat)

% Create ylabel
ylabel({'Number of outliers'});
% Create title
title({'Number of outliers'});

if plotSettings.separateFigs
    % Set theme and font size
    theme(gcf,"light")
    fontsize(plotSettings.fontSize, "points");
end


% Plot number of clusters
if plotSettings.separateFigs
    figure
else
    subplot(2,2,2);
end

plotStatData(clusterCountsStat)
setAxis(clusterCountsStat)

% Create ylabel
ylabel({'number of clusters'});
% Create title
title({'Number of clusters'});

if plotSettings.separateFigs
    % Set theme and font size
    theme(gcf,"light")
    fontsize(plotSettings.fontSize, "points");
end


% Plot cluster sizes
if plotSettings.separateFigs
    figure
else
    subplot(2,2,3);
end

plotStatData(clusterSizesStat)
setAxis(clusterSizesStat)

% Create ylabel
ylabel({'cluster size'});
% Create title
title({'Cluster sizes'});

if plotSettings.separateFigs
    % Set theme and font size
    theme(gcf,"light")
    fontsize(plotSettings.fontSize, "points");
end


%plot percentage of cluster-free outcomes
if plotSettings.separateFigs
    figure
else
    subplot(2,2,4);
end

plotStatData(clusterFreePercentages)
setAxis(clusterFreePercentages)

% Create ylabel
ylabel({'percentage'});
% Create title
title({'% of cluster-free outcomes'});


% Font size
fontsize(plotSettings.fontSize, "points");

% Set light theme
theme(gcf,"light")

function [] = setAxis(data)

    % First, get the upper and lower bounds for the displayed values
    upper = 0;
    lower = 0;

    isPercentages = isfloat(data);
    
    % Get it for percentages
    if isPercentages  
        upper = 100;
        lower = 0;
    % Get it for statistical measures
    else
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
    
        if plotSettings.deviations
            upper = max(upper, max(data.means + sqrt(data.variances)));
            lower = min(lower, min(data.means - sqrt(data.variances)));
        end
    end

    diff = upper - lower;
    upper = upper + percMargin * diff;
    lower = lower - percMargin * diff;

    % Set the axes
    axis([0, MCCount+1, lower, upper])

    % Create xlabel
    xlabel({plotSettings.xLabel});
    % Set xticks & xticklabels
    xticks(1:length(plotSettings.xTickLabels))
    xticklabels(plotSettings.xTickLabels)

    % Show grid
    grid on

end

function [] = plotStatData(data)

    isPercentages = isfloat(data);
    
    % Plot data as cluster free percentages
    if isPercentages    
        plot(xPoints, clusterFreePercentages, clustFreePercPlot);
        return
    end

    % Otherwise - plot data as statistical measures
    hold on

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
                w = data.histograms{i}(j) / histMax;
                rectangle('Position',[ ...
                    i - w/2, ...
                    data.minima(i) + j - 1.5, ...
                    w, ...
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

    if plotSettings.deviations
        % Plot deviations around means
        plot(xPoints,[ ...
            data.means + sqrt(data.variances); ...
            data.means - sqrt(data.variances) ...
            ].', ...
            devPlot)
    end

    hold off

end

end