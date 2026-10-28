clearvars
clc

[par,fn1,fn2] = makeDelayExpsParamsSmallstep(1, 0.09);
xData1 = gatherMCData(par,fn1,fn2,78);
[par,fn1,fn2] = makeDelayExpsParams(1, 90);
xData2 = gatherMCData(par,fn1,fn2,78);
