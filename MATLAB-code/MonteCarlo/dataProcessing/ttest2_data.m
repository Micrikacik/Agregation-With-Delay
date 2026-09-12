[oC1,~,cC1,~,cS1] = processMCData(xData1);
[oC2,~,cC2,~,cS2] = processMCData(xData2);

[h_oC,p_oC] = ttest2(oC1, oC2, "Vartype", "unequal");
[h_cC,p_cC] = ttest2(cC1, cC2, "Vartype", "unequal");
[h_cS,p_cS] = ttest2(cS1{1}, cS2{1}, "Vartype", "unequal");