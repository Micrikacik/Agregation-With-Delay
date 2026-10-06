function [laplace, distance] = prepareMatrices(gridDims, stencil)

arguments
    gridDims (1,:) double {mustBeInteger, mustBePositive}
    stencil double
end

if length(size(stencil)) ~= length(gridDims)
    error("Wrong stencil format.");
end

stencilPointCount = numel(stencil);
stencilDims = size(stencil);
c_gridDims = num2cell(gridDims);
c_stencilDimsHalves = num2cell(ceil(stencilDims/2));
gridPointCount = prod(gridDims);
stepLengths = 1 ./ gridDims;

gridMod = @(indices) mod(indices-1, gridDims)+1;
vecMod = @(index) mod(index-1, gridPointCount)+1;

gridMod = @(indices) cellfun(@(i, d) mod(i-1,d)+1, indices, c_gridDims, UniformOutput=false);

function index = gridToVec(indices)
    index = sub2ind(gridDims, indices{:});
end

function [indices] = vecToGrid(index)
    indices = cell(size(gridDims));
    [indices{:}] = ind2sub(gridDims, index);
end

function [indices] = vecToStenc(index)
    indices = cell(size(gridDims));
    [indices{:}] = ind2sub(stencilDims, index);
end

function [indices] = stencToStencCen(indices)
    indices = cellfun(@(i, s_d) i - s_d, indices, c_stencilDimsHalves, UniformOutput=false);
end

% N = 100;
% M = 100;
% L = N*M;
% 
% dx = 1/N;

%The periodic coordinate functions
fi=@(i) mod(i-1,N)+1;
fj=@(j) mod(j-1,M)+1;
fl=@(l) mod(l-1,M*N)+1;

lij=@(i,j) fi(i)+(fj(j)-1)*N; % converts (i,j) into (l), which is than used to acces matrix in the array sense

jl=@(l) ceil(fl(l)/N); % converts (l) into (j) in (i,j)
il=@(l) fl(l) - (jl(l)-1)*N; % converts (l) into (i) in (i,j)

disti=@(i1,i2) min(abs(fi(i1)-fi(i2)), N-abs(fi(i1)-fi(i2))); % torus distance in MATRIX coordinates
distj=@(j1,j2) min(abs(fj(j1)-fj(j2)), M-abs(fj(j1)-fj(j2))); % torus distance in MATRIX coordinates


%Prepare the matrix A
SeyeL = sparse(1:gridPointCount, 1:gridPointCount, 1);
A = 4 * SeyeL;

for g_index = 1:gridPointCount
    for s_index = 1:stencilPointCount
        s_indices = vecToStenc(s_index);
        s_indicesCen = stencToStencCen(s_indices);
        g_indices = vecToGrid(g_index);
        n_indices = gridMod(cellfun(@plus, s_indicesCen, g_indices, UniformOutput=false));
        A(g_index, gridToVec(n_indices)) = A(g_index, gridToVec(n_indices)) + stencil(s_indices{:});
    end
end

laplace = A;
return
%Save the matrix A
% fname = sprintf('A_%dx%d.mat',N,M);
% save(fname,'A');
% 
% clear A

%Prepare the distance matrix
distM = zeros(L,L);

for l=1:L
    for k=l:L
        distM(l,k) = sqrt((disti(il(l),il(k))^2 + distj(jl(l),jl(k))^2))*dx;
        distM(k,l) = distM(l,k);
    end
end

%Save the matrix W
fname = sprintf('distM_%dx%d.mat',N,M);
save(fname,'distM');

end