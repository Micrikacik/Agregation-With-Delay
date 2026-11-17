function noise = randSmoothFFT(gridDims, centralFreqs)

% Creates a periodic "continuous" noise with dimensions 'gridDims'.
%
% INPUT:
%   gridDims (positive integer) - a vector of dimensions (or just size) of
%       the resulting noise.
%   centralFreqs (non-negative float) - a vector of most influencial
%       frequencies in the noise. The larger the 'centralFreqs(d)', the more
%       oscillating will the noise be in dimension 'd'.
%       If 'centralFreqs(d) == 1', then there is usually one larger hill 
%       and one larger canyon in that direction.
%       If 'centralFreqs(d) > 1', then the number of hills and canyons
%       grows and their size decreases.
%       If 'centralFreqs(d) < 1', then only a part of a hill or a canyon is
%       created in direction 'd', and the smaller the centralFreqs(d) the 
%       more "stretched"/"closer to a constant" that direction becomes.
%       If 'centralFreqs(d) == 0', then the noise is constant in the
%       direction 'd'.
%       If the length  of 'centralFreqs' is not the same as the length of
%       'gridDims', then it is either truncated or filled up with ones.
%
% OUTPUT:
%   noise (non-negative float) - matrix with size 'gridDims' containing values
%       between 0 and 1, representing "continuous" periodic noise
%       (periodic only in the user defined directions, see 'centralFreqs').

arguments
    gridDims double {mustBeVector, mustBeInteger, mustBePositive, mustBeFinite, mustBeReal, mustBeNonempty}
    centralFreqs double {mustBeVector, mustBeNonnegative, mustBeFinite, mustBeReal} = ones(size(gridDims ))
end

nd = numel(gridDims);
ndf = numel(centralFreqs);

% Make them same length
if ndf > nd
    centralFreqs = centralFreqs(1:nd);
end
if ndf < nd
    centralFreqs = [centralFreqs(:), ones(1, nd - ndf)];
end

% Compute dimensions of the noise we use FFT on
noiseDims = gridDims;
cfe0_idxs = centralFreqs <= 0;
noiseDims(cfe0_idxs) = 1;
cflt1_idxs = (centralFreqs < 1) & ~cfe0_idxs;
noiseDims(cflt1_idxs) = round(gridDims(cflt1_idxs) ./ centralFreqs(cflt1_idxs));

% Construct frequency coordinates
freq_sqrd = cell(1, nd);
for d = 1:nd
    n = noiseDims(d);
    f = [0:floor((n-1)/2), -floor(n/2):-1];
    freq_sqrd{d} = f.^2;
end

% Construct squared radial frequency 
% (i.e. each array element contains its distance squared from (0,0))
F = cell(1, nd);
[F{:}] = ndgrid(freq_sqrd{:});

if nd == 1
    freqAddup = zeros(noiseDims, 1);
else
    freqAddup = zeros(noiseDims);
end

for d = 1:nd
    c = centralFreqs(d);
    if c == 0
        continue
    end
    freqAddup = freqAddup + F{d} / c^2;
end

% Construct the filtering matrix
H = 1 ./ (abs(freqAddup - 1) + 1);

% Make noise
if nd == 1
    whiteNoise = rand(noiseDims, 1);
else
    whiteNoise = rand(noiseDims);
end


% Filter frequencies
noise = real(ifftn(fftn(whiteNoise).*H));

% Normalize
noise = noise - min(noise, [], "all");
noiseMax = max(abs(noise), [], "all");
if noiseMax > 0
    noise = noise ./ noiseMax;
end

% Expand noise in the constant dimensions
if ~isequal(noiseDims, gridDims)
    repeat = gridDims;
    repeat(~cfe0_idxs) = 1;
    noise = repmat(noise, repeat);
end

% Cut noise in the nonperiodic (fractional) dimensions
idx = arrayfun(@(n) 1:n, gridDims, 'UniformOutput', false);
noise = noise(idx{:});

end
