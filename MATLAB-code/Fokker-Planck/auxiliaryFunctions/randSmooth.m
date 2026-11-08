function noise = randSmooth(gridDims)

nd = numel(gridDims);
minDim = min(gridDims);

if isscalar(gridDims)
    gridDims = [gridDims, 1];
end

noise = zeros(gridDims);
randCount = 2;
dimOnes = ones([1,nd]);

while randCount < minDim / 3
    randDims = dimOnes * randCount;
    if nd == 1
        randDims = [randDims, 1];
    end
    baseRand = rand(circshift(randDims, nd-1));

    for d = 1:nd

        % baseRand is already shifted so its 2nd dimension is d
        baseRand = paddata(baseRand, size(baseRand, nd) + 5, "Pattern", "circular", "Dimension", nd);

        interpX = 1:size(baseRand, nd); % = randDims(d)
        interpQuery = linspace(3, size(baseRand, nd) - 2, gridDims(d) + 1);
        interpQuery = interpQuery(1:end-1);

        newRand = makima(interpX, baseRand, interpQuery);
        
        if d < nd
            baseRand = shiftdim(newRand, nd-1);
        end
    end
    if nd == 1
        noise = noise + shiftdim(newRand, 1) / randCount;
    else
        noise = noise + newRand / randCount;
    end
    randCount = randCount * 2;
end

noise = noise - min(noise, [], "all");
noise = noise ./ max(noise, [], "all");


% function noise = randSmooth(gridDims)
% 
% nd = numel(gridDims);
% minDim = min(gridDims);
% 
% if isscalar(gridDims)
%     gridDims = [gridDims, 1];
% end
% 
% noise = zeros(gridDims);
% randCount = 2;
% dimOnes = ones([1,nd]);
% 
% indices = cell(1,nd);
% indices_d = cell(1, nd-1);
% 
% while randCount < minDim / 3
%     randDims = dimOnes * randCount;
%     if nd == 1
%         randDims = [randDims, 1];
%     end
%     rV = prod(randDims);
%     baseRand = rand(randDims);
% 
%     for d = 1:nd
%         rV_d = rV/randDims(d);
%         nd_d = nd-1;
%         [indices_d{:}] = ind2sub(randDims, ( 1:(rV_d) ).');
%         coords = zeros(rV_d, nd_d);
%         for dd = 1:nd_d
%             coords(:,dd) = indices_d{dd};
%         end
% 
%         baseRand = paddata(baseRand, size(baseRand, d) + 5, "Pattern", "circular", "Dimension", d);
%         newRandDims = dimOnes * randCount;
%         if nd == 1
%             newRandDims = [newRandDims, 1];
%         end
%         newRandDims(d) = gridDims(d);
%         newRand = zeros(newRandDims);
% 
%         interpX = 1:size(baseRand, d);
%         interpQuery = linspace(3, size(baseRand, d) - 2, gridDims(d) + 1);
% 
%         for index = 1:rV_d
%             for dd = 1:nd_d
%                 if dd < d
%                     indices{dd} = coords(index, dd);
%                 else
%                     indices{dd+1} = coords(index, dd);
%                 end
%             end
%             indices{d} = 1:randDims(d) + 5;
% 
%             randVec = baseRand(indices{:});
%             noiseVec = makima(interpX, randVec, interpQuery);
%             indices{d} = 1:gridDims(d);
% 
%             if d == 1
%                 randShift = randi([0, gridDims(d)-1]);
%                 newRand(indices{:}) = circshift(noiseVec(1:end-1), randShift);
%             else
%                 newRand(indices{:}) = noiseVec(1:end-1);
%             end
%         end
%         baseRand = newRand;
%         randDims = size(baseRand);
%         rV = prod(randDims);
%     end
%     noise = noise + baseRand / randCount;
%     randCount = randCount * 2;
% end
% 
% noise = noise - min(noise, [], "all");
% noise = noise ./ max(noise, [], "all");