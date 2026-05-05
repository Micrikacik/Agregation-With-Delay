% Navrhuji tedy pro vsechny nadchazejici simulace zafixovat N=400, a pro
% 1D zvolit int_r=0.025, pro 2D zvolit int_r=0.05. To jsou stejne volby
% jako v clanku s memory.

%% Reaction test - delay length
jt = 10;
N = 400;
params = cell(jt,1);
seed = 1;
dt = 1e-3;
delayType = "Reaction";
stepCount = 200000;
d = 2;
for j = 1:jt
    params{j} = struct("N", N, "d", d, "stepDelay", 20 * (j+jt), ...
        "dt", dt, "stepRecMod", -1, "rngSeed", seed, "waitForConf", false, ...
        "stepPlotMod", -2, "delayType", delayType, "stepCount", stepCount);
end

% X = cell(jt,1);
% parfor j = 1:jt
%     X{j} = aggWithDelay(params{j});
% end
% 
% fprintf("Experiments finished. \n")
% 
% save("tempTests/testResults.mat","X","params")

%% Reaction test limit

params.x0 = X;
params.xInitHist = hist;
params.rngSetts = sett;
params.stepCount = 50000; % sum = 300 000

F = parfeval(backgroundPool,@aggWithDelay,3,params);

%% 1D aggregation test - radius of density

N = 400;
d = 1;
seed = 1;
dt = 1e-3;
delayType = "Reaction";
stepDelay = 10;
stepCount = 5000;
params{1} = struct("N", N, "d", d, "stepDelay", stepDelay, ...
        "dt", dt, "stepRecMod", 5, "rngSeed", seed, "waitForConf", false, ...
        "stepPlotMod", 5, "delayType", delayType, "stepCount", stepCount);

int_r = 0.025^(1/d);         % interaction radius
kappa = pi^(d/2) / gamma(d / 2 + 1);  % volume of a unit d-ball
W_norm = kappa * int_r^d;   % W_norm to normalize W
W = @(D) (D <= int_r) ./ W_norm;    % default handle of the weight function
params{2} = params{1};
params{2}.W = W;

X = cell(2,1);
for j = 1:2
    figure(j)
    X{j} = aggWithDelay(params{j});
end

fprintf("Experiments finished. \n")

save("tempTests/radiusTestResults.mat","X","params")

%% Compare larger box

N = 400;
params = cell(2,1);
seed = 1;
dt = 1e-2;
delayType = "None";
stepCount = 500;
d = 2;
params{1} = struct("N", N, "d", d, ...
        "dt", dt, "stepRecMod", -1, "rngSeed", seed, "waitForConf", false, ...
        "stepPlotMod", 10, "delayType", delayType, "stepCount", stepCount);
params{2} = params{1};
params{2}.dims = [2,2];
params{2}.N = 1600;

X = cell(2,1);

figure(1)
X{1} = aggWithDelay(params{1});

figure(2)
X{2} = aggWithDelay(params{2});

save("tempTests/sizeTestResults.mat","X","params")

%% Compare BCs

N = 400;
params = cell(2,1);
seed = 1;
dt = 1e-2;
delayType = "None";
stepCount = 5000;
d = 2;
interRad = 0.15;
kappa = pi^(d/2) / gamma(d / 2 + 1);  % volume of a unit d-ball
W_norm = kappa * interRad^d;   % W_norm to normalize W
W = @(D) (D <= interRad) ./ W_norm;    % default handle of the weight function
params{1} = struct("N", N, "d", d, ...
        "dt", dt, "stepRecMod", -1, "rngSeed", seed, "waitForConf", false, ...
        "stepPlotMod", 100, "delayType", delayType, "stepCount", stepCount, "W", W);
params{2} = params{1};
params{2}.boundConds = "Reflective";

X = cell(2,1);

figure(1)
X{1} = aggWithDelay(params{1});

figure(2)
X{2} = aggWithDelay(params{2});

save("tempTests/BCsTestResults.mat","X","params")

%%
N = 400;
dt = 1e-3;
delayType = "Reaction";
stepCount = 300000;
d = 2;
params = struct("N", N, "d", d, "stepDelay", 100, ...
        "dt", dt, "stepRecMod", -1, "waitForConf", false, ...
        "stepPlotMod", -2, "delayType", delayType, "stepCount", stepCount);

aggWithDelay(params)