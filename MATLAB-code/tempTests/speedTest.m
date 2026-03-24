delayType = "Reaction";
T = 5e4;
dt = 1e-3;
d = 2;
p = struct("T",T,"delayType",delayType,"dt",dt,"stepDelay",30,"stepPlotMod",-1,"stepRecMod",-1,"waitForConf",false,"thetaOccurMod",1);

tic
[~,~,thOc] = aggWithDelaySpeedup(p);
time = toc; 