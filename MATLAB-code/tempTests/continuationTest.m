p.rngSeed = 1;
p.stepPlotMod = -2;
p.waitForConf = false;
p.T = 14;
q = p;
q.T = 8;
Xp = aggWithDelay(p);
[Xq,histq,rngq] = aggWithDelay(q);
r = q;
r.T = 6;
r.x0 = Xq;
r.xInitHist = histq;
r.rngSetts = rngq;
temp = rand(1); % to offset randomness
Xr = aggWithDelay(r);
norm(Xp-Xr)

