function [flucts] = makeFlucts(maxFlucts, N)

flucts = maxFlucts * 2 * (rand(N, 1) - 0.5);
flucts = flucts - mean(flucts);
