params.stepDelay = 0;
params.x0 = 0.5 * ones(2,2);
%params.x0 = rand(400,2);
params.xInitHist = repmat(params.x0,[1,1,params.stepDelay]);
params.stepPlotMod = 1;
params.stepCount = 50000;
params.dt = 0.001;
params.intRad = 0.05;
params.rngSeed = 1;
figure(1)
aggWithDelaySpeedupAdapt(params); 
% figure(2)
% aggWithDelaySpeedup(params);