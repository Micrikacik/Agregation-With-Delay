
if isempty(gcp('nocreate'))
    parpool(14);
end

p = gcp("nocreate");

[pars,fn1,fn2] = makeDelayExpsParamsFP(2,30:30:420);
results = cell(1,length(pars));

% rho0AndInitHist = genRhoInitHistSmooth(1e-3,ones(10000,1),[100,100],100001,[2,2,2*100001*1e-3]);
% rho0AndInitHist = reshape(rho0AndInitHist, [100,100,100001]);
% 
% for i = 1:length(pars)
%     expParams = pars(i);
%     pars(i).rhoInitHist = rho0AndInitHist(:,:,2:expParams.stepDelay + 1);
%     pars(i).rho0 = rho0AndInitHist(:,:,1);
% end

group = 0;

while group * p.NumWorkers < length(pars)
    parfor i = (group * p.NumWorkers + 1):((group+1) * p.NumWorkers)
        expParams = pars(i);
        results{i} = FokkerPlanck(expParams)
    end
    
    for i = (group * p.NumWorkers + 1):((group+1) * p.NumWorkers)
        expParams = pars(i);
        folder = fn1(expParams, i);
        postfix = fn2(expParams, i);
        if ~isfolder(folder)
            mkdir(folder)
        end
        result = results{i};
        save(FPFilePath(folder,"FPData",postfix), "result", "expParams")
    end
    group = group + 1;
end