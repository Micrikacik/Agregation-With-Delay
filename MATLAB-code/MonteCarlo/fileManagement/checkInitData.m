d = 2;
delayType = "Transmission";
fileName = "MCData";



for stepDelay = [10000:10000:100000]
    path = MCFilePath(MCFolderPath(delayType,d,stepDelay),fileName,MCFilePostfix(delayType,d,stepDelay));
    results = load(path).results;
    fprintf("Started %i\n", stepDelay)
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
    fprintf("\nfinished %i, same: %i\n", stepDelay, same)
end

% for stepDelay1 = 0:30:420
%     fprintf("Start %i\n", stepDelay1)
% for stepDelay2 = stepDelay1+30:30:420
%     path1 = MCFilePath(MCFolderPath(delayType,d,stepDelay1),fileName,MCFilePostfix(delayType,d,stepDelay1));
%     results1 = load(path1).results; % Load the data from the specified path
%     path2 = MCFilePath(MCFolderPath(delayType,d,stepDelay2),fileName,MCFilePostfix(delayType,d,stepDelay2));
%     results2 = load(path2).results; % Load the data from the specified path
%     same(stepDelay/30+1) = 0;
%     for i = 1:100
%         for j = 1:100
%             if norm(results1{i}.xRec(:,:,1) - results2{j}.xRec(:,:,1), "fro") < 0.001
%             fprintf("(%i,%i,%i,%i)", stepDelay1, i, stepDelay2, j)
%             same(stepDelay/30+1) = same(stepDelay/30+1) + 1;
%             end
%         end
%     end
%     fprintf("\nfinished %i, same: %i\n", stepDelay2, same(stepDelay/30+1))
% end
% end