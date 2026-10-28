function stencil = makeStencil(d)

if d == 1
    stencil = zeros(3,1);
else
    stencil = zeros(repmat(3,1,d));
end
middle = cell(1,d);

for index = 1:length(middle)
    middle{index} = 2;
end

stencil(middle{:}) = 2*d;

for dim = 1:d
    indices = middle;
    indices{dim} = indices{dim} + 1;
    stencil(indices{:}) = -1;
    indices = middle;
    indices{dim} = indices{dim} - 1;
    stencil(indices{:}) = -1;
end