clearvars
clc

%% default
[par,fn1,fn2] = makeDelayExpsParamsSmallstep(1, 9000);
xData1 = gatherMCData(par,fn1,fn2,78);
[par,fn1,fn2] = makeDelayExpsParams(1, 90, "Reaction");
xData2 = gatherMCData(par,fn1,fn2,100);

ttest2_data(xData1, xData2);

%% noise to positions
noise = 0.001;

xData = xData2;
xData2 = xData + randn(size(xData2)) * noise;

ttest2_data(xData1, xData2);

%% noise to stochastic measures
noise = 0.001;

xData2 = xData;

[oC1,~,cC1,~,cS1] = processMCData(xData1, false);
[oC2,~,cC2,~,cS2] = processMCData(xData2, false);
cS1 = cS1{1};
cS2 = cS2{1};

oC2 = oC2 + randn(size(oC2)) * noise;
cC2 = cC2 + randn(size(cC2)) * noise;
cS2 = cS2 + randn(size(cS2)) * noise;

[h_oC, p_oC] = ttest2(oC1, oC2, "Vartype", "unequal");
[h_cC, p_cC] = ttest2(cC1, cC2, "Vartype", "unequal");
[h_cS, p_cS] = ttest2(cS1, cS2, "Vartype", "unequal");

results = table([h_oC; h_cC; h_cS], [p_oC; p_cC; p_cS], ...
    'VariableNames', {'TestResult', 'P-Value'}, ...
    'RowNames', {'oC', 'cC', 'cS'});

fig = figure;
uitable(fig, "Data", results{:,:}, ...
    "ColumnName", results.Properties.VariableNames, ...
    "RowName", results.Properties.RowNames,    ...
    "Position", [200 200 350 170]);
fontsize(getFontSize(), "points")

%% bounded gauss sensitivity
N = 1000;
mean1 = 0;
mean2 = 0;
varSqrt1 = 1;
varSqrt2 = 1;
G1 = randn([1,N]) * varSqrt1 + mean1;
G2 = randn([1,N]) * varSqrt2 + mean2;

%% control group
figure
hold on
histogram(G1, "Normalization", "pdf")
histogram(G2, "Normalization", "pdf")
hold off

[h_G, p_G] = ttest2(G1, G2);

results = table(h_G, p_G, ...
    'VariableNames', {'TestResult', 'P-Value'}, ...
    'RowNames', {'G'});

fig = figure;
uitable(fig, "Data", results{:,:}, ...
    "ColumnName", results.Properties.VariableNames, ...
    "RowName", results.Properties.RowNames,    ...
    "Position", [200 200 350 170]);
fontsize(getFontSize(), "points")

%% test group
upperBound = 2;
lowerBound = -2;

G2(G2 > upperBound) = upperBound;
G2(G2 < lowerBound) = lowerBound;

figure
hold on
histogram(G1, "Normalization", "pdf")
histogram(G2, "Normalization", "pdf")
hold off

[h_G, p_G] = ttest2(G1, G2, "Vartype", "unequal");

results = table(h_G, p_G, ...
    'VariableNames', {'TestResult', 'P-Value'}, ...
    'RowNames', {'G'});

fig = figure;
uitable(fig, "Data", results{:,:}, ...
    "ColumnName", results.Properties.VariableNames, ...
    "RowName", results.Properties.RowNames,    ...
    "Position", [200 200 350 170]);
fontsize(getFontSize(), "points")

%% 

rng(1);

N = 100;

x = randn(N,1);
y = randn(N,1);

[h1,p1,ci1,stats1] = ttest2(x,y,'Vartype','equal');
[h2,p2,ci2,stats2] = ttest2(x,y,'Vartype','unequal');

fprintf('Equal variance:\n');
fprintf('t = %.6f, df = %.6f, p = %.6f\n', ...
    stats1.tstat, stats1.df, p1);

fprintf('Unequal variance:\n');
fprintf('t = %.6f, df = %.6f, p = %.6f\n', ...
    stats2.tstat, stats2.df, p2);