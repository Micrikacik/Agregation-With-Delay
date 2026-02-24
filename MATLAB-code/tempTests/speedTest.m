delayType = "Reaction";
T = 1e5;
dt = 1e-3;
d = 2;
p = struct("T",T,"delayType",delayType,"dt",dt,"stepPlotMod",-1,"stepRecMod",-1,"waitForConf",false);

tic
aggWithDelay(p);
time = toc; 