function noise = smoothRand(gridDims)

nd = numel(gridDims);
gridPeriod = min(gridDims);

if isscalar(gridDims)
    gridDims = [gridDims, 1];
end

noise = zeros(gridDims); 
noisePart = ones(gridDims); 

while gridPeriod > 4
    for d = 1:nd
        ls = (1:gridDims(d)).' / gridPeriod * 2*pi;
        rOffset = rand(1) * gridDims(d) / gridPeriod * 2*pi;
        noiseElem = sin(ls + rOffset) * log(gridPeriod);
        noisePart = noisePart .* reshape(noiseElem, [ones(1, d-1), numel(noiseElem), 1]);
    end
    noise = noise + noisePart;
    gridPeriod = gridPeriod / 2;
end

return