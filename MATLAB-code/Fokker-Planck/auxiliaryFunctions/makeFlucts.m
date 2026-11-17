function [flucts] = makeFlucts(maxFlucts, gridDims)

if isscalar(gridDims)
    flucts = rand(gridDims, 1);
else 
    flucts = rand(gridDims);
end

flucts = flucts - mean(flucts, "all");
maxf = max(abs(flucts), [], "all");
if maxf > 0
    flucts = flucts / maxf * maxFlucts;
end
