function [minpts] = getMinClusterSize(N, d, volume)

baseMinpts = 17;
minpts = round(baseMinpts * (N / volume / 400));