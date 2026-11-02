function [flucts] = makeSmoothFlucts(maxFluc, gridDims)

flucts = rands(gridDims);
flucts = flucts - mean(flucts);
flucts = maxFluc * flucts / max(flucts, [], "all");
