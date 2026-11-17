function flucts = makeFluctsSmooth(maxFlucts, gridDims, centralFreqs)

flucts = randSmoothFFT(gridDims, centralFreqs);
flucts = flucts - mean(flucts, "all");
maxf = max(abs(flucts), [], "all");
if maxf > 0
    flucts = flucts / maxf * maxFlucts;
end
