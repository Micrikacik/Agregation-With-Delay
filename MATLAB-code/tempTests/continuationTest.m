rng('shuffle')
p.stepPlotMod = -2;
p.waitForConf = false;
p.stepCount = 14;
q = p;
q.stepCount = 8;
Xp = aggWithDelay(p);
[Xq,~,histq,rngq] = aggWithDelay(q);
r = q;
r.stepCount = 6;
r.x0 = Xq;
r.xInitHist = histq;
r.rngSetts = rngq;
temp = rand(1); % to offset randomness
Xr = aggWithDelay(r);
norm(Xp-Xr)

