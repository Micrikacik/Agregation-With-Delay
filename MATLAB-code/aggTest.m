% Initialize experiments settings
iN = 3;
jt = 10;
params = cell(iN,jt);
seed = int16(1);
dt = 1e-2;
delayType = "Transmission";
T = int16(300);
d = int16(2);
for i = 1:iN
    for j = 1:jt+1
        params{i,j} = struct("N", int16(100 * (5 * i - 4)), "d", d, "stepDelay", int16(20 * (j-1)), ...
            "dt", dt, "stepRecMod", int16(-1), "rngSeed", seed, "waitForConf", false, ...
            "stepPlotMod", int16(-2), "delayType", delayType, "T", T);
    end
end

X = cell(iN,jt);
for i = 1:iN
    parfor j = 1:jt+1
        X{i,j} = aggWithDelay(params{i,j});
    end
    fprintf("FINISHED %i. out of %i.\n\n", i, iN)
end

fprintf("Experiments finished. \n")

save("testResults.mat","X","params")