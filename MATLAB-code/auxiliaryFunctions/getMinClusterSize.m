function minpts = getMinClusterSize(N, d, volume)

% Returns the default 'minpts' parameter used in the 'dbscan' method
%
% INPUT:
%   N (positive integer) - number of agents.
%   d (positive integer) - dimension. (Currently not used in the implementation.)
%   volume (positive float) - volume of the simulation domain.
%
% OUTPUT:
%   minpts (positive integer) - default value of the parameter 'minpts' for
%       the method 'dbscan'.

arguments
    N (1,1) double {mustBeInteger, mustBePositive, mustBeReal, mustBeFinite}
    d (1,1) double {mustBeInteger, mustBePositive, mustBeReal, mustBeFinite}
    volume (1,1) double {mustBePositive, mustBeReal, mustBeFinite}
end

baseMinpts = 17;
minpts = round(baseMinpts * (N / volume / 400));