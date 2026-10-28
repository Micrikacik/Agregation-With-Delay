function [A, distM] = makeFPMatrices(gridDims, domainDims, stencil)

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

stencilPointCount = numel(stencil);
stencilDims = size(stencil);
c_gridDims = num2cell(gridDims);
c_stencilDimsHalves = num2cell(ceil(stencilDims/2));
gridPointCount = prod(gridDims);

gridMod = @(indices) cellfun(@(i, d) mod(i-1,d)+1, indices, c_gridDims, UniformOutput=false);

function index = gridToVec(indices)
    index = sub2ind(gridDims, indices{:});
end

function [indices] = vecToGrid(index)
    indices = cell(size(gridDims));
    [indices{:}] = ind2sub(gridDims, index);
end

function [indices] = vecToStenc(index)
    indices = cell(size(gridDims));
    [indices{:}] = ind2sub(stencilDims, index);
end

function [indices] = stencToStencCen(indices)
    indices = cellfun(@(i, s_d) i - s_d, indices, c_stencilDimsHalves, UniformOutput=false);
end

%Prepare the matrix A
A = sparse(gridPointCount, gridPointCount);

for g_index = 1:gridPointCount
    g_indices = vecToGrid(g_index);
    for s_index = 1:stencilPointCount
        s_indices = vecToStenc(s_index);
        s_indicesCen = stencToStencCen(s_indices);
        n_indices = gridMod(cellfun(@plus, s_indicesCen, g_indices, UniformOutput=false));
        A(g_index, gridToVec(n_indices)) = A(g_index, gridToVec(n_indices)) + stencil(s_indices{:});
    end
end

%Prepare the distance matrix
distM = zeros(gridPointCount, gridPointCount);

dDim = domainDims ./ gridDims;

dist = zeros(gridDims);
for g_index = 1:gridPointCount
    g_indices = vecToGrid(g_index);
    distsFrom1 = cellfun(@(i) i-1, g_indices);
    dist(g_indices{:}) = sqrt( sum( ( min( distsFrom1, gridDims - distsFrom1 ) .* dDim ) .^ 2 ) );
end

lastShiftIDs = ones(size(gridDims));

for l = 1:gridPointCount
    distM(l,:) = dist(:);
    g_indices = vecToGrid(l+1);
    newShiftIDs = [g_indices{:}];
    shift = newShiftIDs - lastShiftIDs;
    lastShiftIDs = newShiftIDs;
    dist = circshift(dist, shift);
end

end