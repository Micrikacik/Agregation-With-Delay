function xInitHist = generInitHist(x,dt,stepDelay)

[N, d] = size(x);

xInitHist = zeros([N, d, stepDelay]);

xInitHist(:,:,1) = x - sqrt(dt) * rand(N, d);

for i = 2:stepDelay
    xInitHist(:,:,i) = xInitHist(:,:,i-1) - sqrt(dt) * randn(N, d);
end
