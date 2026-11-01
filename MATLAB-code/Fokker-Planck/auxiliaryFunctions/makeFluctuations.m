function [fluctuations] = makeFluctuations(maxFluc, V)

fluctuations = maxFluc * 2 * (rand(V, 1) - 0.5);
fluctuations = fluctuations - mean(fluctuations);
