function discrVol = discretizedBallVolume(rad, dDim)

arguments
    rad (1,1) double {mustBePositive}
    dDim (1,:) double {mustBePositive}
end

d = numel(dDim);

if isscalar(dDim)
    grid = zeros([ceil(rad ./ dDim) + 1, 1]);
else
    grid = zeros(2 * ceil(rad ./ dDim) + 1); % might be done more memory effectively
end

gridDims = size(grid);
N = prod(gridDims);

indices = cell(1, d);
[indices{:}] = ind2sub(gridDims, (1:N).');
coords_dist = zeros(N, d);
for d = 1:d
    coords_dist(:, d) = indices{d} - 1;
end

for i = 1:N
    delta = min(coords_dist(i,:), gridDims - coords_dist(i,:));
    grid(i) = sum((delta .* dDim) .^ 2);
end

grid = sqrt(grid);

grid = grid < rad;

discrVol = sum(grid, "all");