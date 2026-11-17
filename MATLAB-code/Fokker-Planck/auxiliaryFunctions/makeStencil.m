function stencil = makeStencil(gridDims, domainDims)

% Creates stencil for finite volumes/differences scheme on domain with
% 'domainDims' discretized by grid with 'gridDims'. Scheme is simple second
% difference (x(i+1) - 2x(i) + x(i-1)) / dx^2

% gridDims - number of steps in each dimension = [N1, N2, ..., Nd]

d = length(gridDims);

if length(domainDims) ~= d
    error("Input vectors do not have the same size.")
end

if d == 1
    stencil = zeros(3,1);
else
    stencil = zeros(3*ones(1, d));
end

center = cell(1,d);

for i = 1:length(center)
    center{i} = 2;
end

multip = (gridDims ./ domainDims) .^ 2; % dDims = domainDims ./ gridDims and we divide by it

stencil(center{:}) = sum(-2 * multip);

for dim = 1:d
    indices = center;
    indices{dim} = indices{dim} + 1;
    stencil(indices{:}) = multip(dim);
    indices = center;
    indices{dim} = indices{dim} - 1;
    stencil(indices{:}) = multip(dim);
end