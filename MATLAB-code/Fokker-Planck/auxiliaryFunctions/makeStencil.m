function stencil = makeStencil(d, gridDims)

% d - dimension
% gridDims - number of steps in each dimension = [N1, N2, ..., Nd]

if d == 1
    stencil = zeros(3,1);
else
    stencil = zeros(repmat(3,1,d));
end
middle = cell(1,d);

for index = 1:length(middle)
    middle{index} = 2;
end

gridDimsSqrd = gridDims .^ 2;

stencil(middle{:}) = sum(2 * gridDimsSqrd);

for dim = 1:d
    indices = middle;
    indices{dim} = indices{dim} + 1;
    stencil(indices{:}) = -1 * gridDimsSqrd(dim);
    indices = middle;
    indices{dim} = indices{dim} - 1;
    stencil(indices{:}) = -1 * gridDimsSqrd(dim);
end