%%
X = load("tempTests/testResults.mat").X;

jt = size(X,1);
r = cell(size(X));
for j = 1:jt
    r{j} = unique(dbscan(torusDistances(X{j}),1.3/sqrt(2*(size(X{j},1))),12,'Distance','precomputed'));
end 

r

%%
X = load("tempTests/radiusTestResults.mat").X;

jt = size(X,1);
r = cell(size(X));
for j = 1:jt
    r{j} = unique(dbscan(torusDistances(X{j}),1.3/sqrt(2*(size(X{j},1))),12,'Distance','precomputed'));
end

r
