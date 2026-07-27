%%%%%%%%%%%%%%%%%%%
%%% TEMP SCRIPT %%%
%%%%%%%%%%%%%%%%%%%

function path = folderPathFunc(d, intRad)
    path = "MonteCarlo/stabilityByRadius";
    path = sprintf("%s/%iD/intRad0p%.0f", path, d, intRad * 100);
end

function postfix = filePostfixFunc(d, intRad)
    postfix = sprintf("%iD_0p%.0f", d, intRad * 100);
end

fn1 = @(params, i_exp) folderPathFunc(params.d, params.intRad);
fn2 = @(params, i_exp) filePostfixFunc(params.d, params.intRad);

for i = 1:length(par)
    param = par(i);
    path = MCFilePath(fn1(param, 1), "MCData", fn2(param, 1));
    dirp = fnb1(param, 1);
    mkdir(dirp)
    pathb = MCFilePath(fnb1(param, 1), "MCData", fnb2(param, 1));
    isf = isfile(path);
    fprintf("%s is %i\n", path, isf)

    movefile(path,pathb)
end