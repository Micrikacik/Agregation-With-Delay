function [A, distM] = prepareMatrices_optimized(gridDims, domainDims, stencil)
%PREPAREMATRICES_OPTIMIZED Faster equivalent of prepareMatrices.
% MATLAB R2025b. Uses sparse triplet assembly for A and pairwise
% periodic coordinate differences for distM.

arguments
    gridDims (1,:) double {mustBeInteger, mustBePositive}
    domainDims (1,:) double {mustBePositive}
    stencil double
end

if isscalar(gridDims)
    gridDims = [gridDims, 1];
end
if isscalar(domainDims)
    domainDims = [domainDims, 1];
end

if length(size(stencil)) ~= length(gridDims)
    error("Wrong stencil format.");
end
if length(domainDims) ~= length(gridDims)
    error("'domainDims' and 'gridDims' do not match.");
end

nd = numel(gridDims);
N = prod(gridDims);
K = numel(stencil);
stencilDims = size(stencil);

% Coordinates of all grid points in MATLAB's column-major linear order.
% We only use 'coords' from this, for quick coordinate conversion.
gridSubs = cell(1, nd);
[gridSubs{:}] = ind2sub(gridDims, (1:N).');
coords = zeros(N, nd);
for d = 1:nd
    coords(:, d) = gridSubs{d};
end

% Coordinates of stencil points and their offsets from the stencil center.
% We only use 'offsets' from this, for quick coordinate conversion.
stencilSubs = cell(1, nd);
[stencilSubs{:}] = ind2sub(stencilDims, (1:K).');
offsets = zeros(K, nd);
for d = 1:nd
    offsets(:, d) = stencilSubs{d} - ceil(stencilDims(d) / 2);
end

% Assemble A in one sparse construction. Repeated (row,column) pairs are
% summed by SPARSE, as in the original repeated updates to A(g_index, j).
rows = repelem((1:N).', K, 1);
cols = zeros(N*K, 1);
stencVals = repmat(stencil(:), N, 1);

for g = 1:N
    idx = (g-1)*K + (1:K);
    neighborSubs = mod(offsets + coords(g, :) - 1, gridDims) + 1;
    neighborCells = num2cell(neighborSubs, 1);
    cols(idx) = sub2ind(gridDims, neighborCells{:});
end
A = sparse(rows, cols, stencVals, N, N);

% Pairwise periodic distances. This directly computes the same toroidal
% distance represented by the original shifted distance array.
dDim = domainDims ./ gridDims;
distM = zeros(N, N);
for d = 1:nd
    delta = abs(coords(:, d) - coords(:, d).');
    delta = min(delta, gridDims(d) - delta) .* dDim(d);
    distM = distM + delta.^2;
end
distM = sqrt(distM);
end
