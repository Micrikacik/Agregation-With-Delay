function [h_oC,p_oC,h_cC,p_cC,h_cS,p_cS] = ttest2_data(xData1, xData2)

[oC1,~,cC1,~,cS1] = processMCData(xData1, false);
[oC2,~,cC2,~,cS2] = processMCData(xData2, false);
cS1 = cS1{1};
cS2 = cS2{1};

[h_oC,p_oC,ci_oC] = ttest2(oC1, oC2, "Vartype", "unequal");
[h_cC,p_cC,ci_cC] = ttest2(cC1, cC2, "Vartype", "unequal");
[h_cS,p_cS,ci_cS] = ttest2(cS1, cS2, "Vartype", "unequal");

results = table([h_oC; h_cC; h_cS], [p_oC; p_cC; p_cS], ...
    [ci_oC(1); ci_cC(1); ci_cS(1)], [mean(oC1) - mean(oC2); mean(cC1) - mean(cC2); mean(cS1) - mean(cS2)] ,[ci_oC(2); ci_cC(2); ci_cS(2)], ...
    'VariableNames', {'TestResult', 'PValue', 'Start of conf. int.', 'Mean difference', 'End of conf. int.'}, ...
    'RowNames', {'oC', 'cC', 'cS'});

fig = figure;
uitable(fig, "Data", results{:,:}, ...
    "ColumnName", results.Properties.VariableNames, ...
    "RowName", results.Properties.RowNames,    ...
    "Position", [200 200 600 170]);
fontsize(getFontSize(), "points")