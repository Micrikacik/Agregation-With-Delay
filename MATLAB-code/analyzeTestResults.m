X = load("testResults.mat").X;

iN = size(X,1);
jt = size(X,2);
r = cell(size(X));
for i = 1:iN
    for j = 1:jt
        r{i,j} = unique(dbscan(torusDistances(X{i,j}),1.3/sqrt(2*(size(X{i,j},1))),12,'Distance','precomputed'));
    end
end

r