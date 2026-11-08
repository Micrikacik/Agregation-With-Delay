function flucts = makeFluctsSmooth(maxFlucts, gridDims)

flucts = randSmooth(gridDims);
flucts = flucts - mean(flucts);
flucts = flucts / max(abs(flucts), [], "all") * maxFlucts;
