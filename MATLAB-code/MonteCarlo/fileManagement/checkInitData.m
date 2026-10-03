%%%%%%%%%%%%%%%%%%%
%%% TEMP SCRIPT %%%
%%%%%%%%%%%%%%%%%%%

fileName = "MCData";

[params, folderPathFunc, filePostfixFunc] = makeDelayExpsParams(2,0:30:420,"Reaction");

totalSame = 0;

for i_exp = 1:length(params)
    par = params(i_exp);

    folderPath = folderPathFunc(params, i_exp);
    postfix = filePostfixFunc(params, i_exp);
    
    path = MCFilePath(folderPath, fileName, postfix);
    results = load(path).results;
    fprintf("Started \n")
    disp(par)
    same = 0;
    for i = 1:length(results)
        for j = i+1:length(results)
            if norm(results(i).xRec(:,:,1) - results(j).xRec(:,:,1), "fro") < 0.001
                fprintf("(%i,%i) ", i, j)
                same = same + 1;
                break
            end
        end
    end
    fprintf("\nFinished \n")
    disp(par)
    fprintf("same: %i\n", same)
    totalSame = totalSame + same;
end

fprintf("\nTotal same: %i\n", totalSame)