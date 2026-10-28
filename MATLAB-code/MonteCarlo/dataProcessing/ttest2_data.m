function [h_oC,p_oC,h_cC,p_cC,h_cS,p_cS] = ttest2_data(xData1, xData2)

[oC1,~,cC1,~,cS1] = processMCData(xData1, false);
[oC2,~,cC2,~,cS2] = processMCData(xData2, false);

[h_oC,p_oC] = ttest2(oC1, oC2, "Vartype", "unequal");
[h_cC,p_cC] = ttest2(cC1, cC2, "Vartype", "unequal");
[h_cS,p_cS] = ttest2(cS1{1}, cS2{1}, "Vartype", "unequal");

results = table([h_oC; h_cC; h_cS], [p_oC; p_cC; p_cS], ...
    'VariableNames', {'TestResult', 'PValue'}, ...
    'RowNames', {'oC', 'cC', 'cS'});

fig = figure;
uitable(fig, "Data", results{:,:}, ...
    "ColumnName", results.Properties.VariableNames, ...
    "RowName", results.Properties.RowNames,    ...
    "Position", [200 200 350 170]);
fontsize(getFontSize(), "points")