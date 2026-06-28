clear
p.stepPlotMod = 100;
p.stepCount = 10000;
p.stepDelay = 240;
figure(1)
x_1 = aggWithDelaySpeedup(p);

p.N = 3600;
p.dims = [3,3];
p.boundConds = "NoBoundary";
figure(2)
x_2 = aggWithDelaySpeedup(p);