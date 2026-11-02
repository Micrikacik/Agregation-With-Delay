function s = smoothRand(gridDims)
    
w = min(gridDims);

if isscalar(gridDims)
    gridDims = [gridDims, 1];
end

nd = numel(gridDims);
s = zeros(gridDims); 
i = 0;

indices = cell(1, nd);
for d = 1:nd
    indices{d} = (1:gridDims(d));
end

while w > 3
    i = i + 1;
    spl = interpn(rand(gridDims), i-1, 'spline');
    if isrow(spl)
        spl = spl.';
    end
    s = s + i * spl(indices{:});
    w = w - ceil(w/2 - 1);
end

s = (s - min(s(:))) ./ (max(s(:)) - min(s(:)));