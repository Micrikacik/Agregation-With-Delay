function [outlierCounts,outlierCountsStat, ...
            clusterCounts,clusterCountsStat, ...
            clusterSizes,clusterSizesStat, ...
            clusterFreePercentage] = processMCData(xData, plotData, epsilon, minpts)

% Processes the data from Monte Carlo simulations, returning statistical 
% measures and plotting their graphs (if user wants). Uses a metric 
% on a unit torus to measure distances, by calling auxiliary function 
% torusDistances.m.
%
%--------------------------------------------------------------------------
%
% INPUT
%   xData (float matrix) - N by d by MCSize by MCCount - matrix of positions 
%       of agents, where 'N' is the number of agents, 'd' is the dimension,
%       'MCSize' is the number of all simulations in one Monte Carlo simulation 
%       and 'MCCount' is the count of all different Monte Carlo simulations.
%       In other words, xData(j_Sim,j,k,l) is the j-th coordinate of the 
%       postition of the j_Sim-th agent from the k-th simulation from the l-th 
%       set of simulations.
%   plotData (logical) - if true, then all the statistical measures are
%       plotted into the current figure.
%   epsilon (positive float) - parameter for dbscan function, it is a
%       cluster search radius. For more details, see the documentation of 
%       dbscan.
%   minpts (positive integer) - parameter for dbscan function, it
%       represents a minimal number of agents that can form a cluster. For
%       more details, see the documentation of dbscan.

arguments
    xData (:,:,:,:) double
    plotData (1,1) logical = true
    epsilon (1,1) double {mustBePositive} = 0.05
    minpts (1,1) double {mustBePositive} = 17
end

minpts = round(minpts);

wBar = waitbar(0,"Processing data...");

%number of agents
N = size(xData,1);

% Dimension
d = size(xData,2);

% Number of runs in each Monte Carlo simulation
MCSize = size(xData,3);

% Number of different Monte Carlo simulations
MCCount = size(xData,4);

% Preallocation
outlierCounts = zeros(MCSize,MCCount);
clusterCounts = zeros(MCSize,MCCount);
clusterSizes = cell(1,MCCount);
clusterFreePercentage = zeros(MCCount,1);

% Go through all Monte Carlo simulations
for i_MC=1:MCCount    
    % Go through all individual simulations in the Monte Carlo simulation
    for j_Sim=1:MCSize
        % Extract the data
        x = xData(:,:,j_Sim,i_MC);

        % Calculate distances over the torus
        dists = torusDistances(x);
        
        % Identify clusters using the DBSCAN method with distances 'dists'
        idx = dbscan(dists,epsilon,minpts,'Distance','precomputed');

        % Number of outliers
        outlierCounts(j_Sim,i_MC) = nnz(idx == -1);

        % Number of clusters
        clusterCounts(j_Sim,i_MC) = numel(unique(idx(idx > 0)));

        % Cluster sizes (how many clusters of give size appeared in i_MC-th Monte Carlo)
        for cluster = 1:clusterCounts(j_Sim,i_MC)
            clusterSizes{i_MC}(end + 1) = nnz(idx == cluster);
        end
    end

    % Percentage of cluster-free outcomes
    clusterFreePercentage(i_MC) = 100 * nnz(clusterCounts(:,i_MC) == 0) / MCSize;

    waitbar(i_MC / MCCount,wBar,"Processing data...")
end

% Fix clusterSizes format
emptyIndices = cellfun(@isempty,clusterSizes);
clusterSizes(emptyIndices) = {0};

% Calculate statistical measures for outlierCountsStat
outlierCountsStat.means = mean(outlierCounts);
outlierCountsStat.maxima = max(outlierCounts);
outlierCountsStat.minima = min(outlierCounts);
outlierCountsStat.variances = var(outlierCounts);
outlierCountsStat.medians = median(outlierCounts);
for col = 1:MCCount
    outlierCountsStat.histograms{col} = histcounts(outlierCounts(:,col), ...
        outlierCountsStat.minima(col)-0.5:1:outlierCountsStat.maxima(col)+0.5);
end

% Calculate statistical measures for clusterCountsStat
clusterCountsStat.means = mean(clusterCounts);
clusterCountsStat.maxima = max(clusterCounts);
clusterCountsStat.minima = min(clusterCounts);
clusterCountsStat.variances = var(clusterCounts);
clusterCountsStat.medians = median(clusterCounts);
for col = 1:MCCount
    clusterCountsStat.histograms{col} = histcounts(clusterCounts(:,col), ...
        clusterCountsStat.minima(col)-0.5:1:clusterCountsStat.maxima(col)+0.5);
end

% Calculate statistical measures for clusterSizesStat
clusterSizesStat.means = cellfun(@mean,clusterSizes);
clusterSizesStat.maxima = cellfun(@max,clusterSizes);
clusterSizesStat.minima = cellfun(@min,clusterSizes);
clusterSizesStat.variances = cellfun(@var,clusterSizes);
clusterSizesStat.medians = cellfun(@median,clusterSizes);
for col = 1:MCCount
    clusterSizesStat.histograms{col} = histcounts(clusterSizes{col}, ...
        clusterSizesStat.minima(col)-0.5:1:clusterSizesStat.maxima(col)+0.5);
end


if plotData == false
    close(wBar)
    return
end

waitbar(1,wBar,"Plotting results...")

plotMCData(outlierCountsStat,clusterCountsStat,clusterSizesStat,clusterFreePercentage,[],N)

close(wBar)