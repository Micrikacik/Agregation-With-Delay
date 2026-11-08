function noise = randSmoothFFT(gridDims, centralFreqs)

nd = numel(gridDims);

if isscalar(gridDims)
    gridDims = [gridDims, 1];
end

whiteNoise = randn(gridDims);

% Construct normalized frequency coordinates
freq = cell(1, nd);

for d = 1:nd
    n = gridDims(d);
    f = mod((0:(n-1)) + n/2, n) - n/2; 
    freq{d} = f.^2;
end

% Construct radial frequency
F = cell(1, nd);
[F{:}] = ndgrid(freq{:});

freqAddup = zeros(gridDims);

for d = 1:nd
    c = centralFreqs(d);
    freqAddup = freqAddup + F{d} / max(c,1)^2 - c^2;
end

% id = num2cell(centralFreqs + 1);
% freqAddup(id{:}) = 1;

figure
imagesc(freqAddup)

% Gaussian low-pass filter
H = 1 ./ (abs(freqAddup) + 1);
figure
imagesc(H)

% Filter in frequency domain
noise = real(ifftn(fftn(whiteNoise).*H));

% Normalize
noise = noise - min(noise, [], "all");
noise = noise ./ max(noise, [], "all");

end

% function noise = randSmoothFFT(gridDims, sigma)
% 
% if isscalar(gridDims)
%     gridDims = [gridDims, 1];
% end
% 
% nd = numel(gridDims);
% 
% % Generate white noise
% whiteNoise = randn(gridDims);
% 
% % Construct normalized frequency coordinates
% freq = cell(1, nd);
% 
% for d = 1:nd
%     n = gridDims(d);
%     freq{d} = [0:ceil(n/2)-1, -floor(n/2):-1] / n;
% end
% 
% % Construct squared radial frequency
% F = cell(1, nd);
% [F{:}] = ndgrid(freq{:});
% 
% freqSquared = zeros(gridDims);
% 
% for d = 1:nd
%     freqSquared = freqSquared + F{d}.^2;
% end
% 
% % Gaussian low-pass filter
% H = exp(-2*pi^2*sigma^2*freqSquared);
% imagesc(H)
% 
% % Filter in frequency domain
% noise = real(ifftn(fftn(whiteNoise).*H));
% 
% % Normalize
% noise = noise - min(noise, [], "all");
% noise = noise ./ max(noise, [], "all");
% 
% end
