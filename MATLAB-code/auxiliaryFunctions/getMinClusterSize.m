function [minpts] = getMinClusterSize(N, volume)

baseMinpts = 17;
minpts = round(baseMinpts * (N / volume / 400));