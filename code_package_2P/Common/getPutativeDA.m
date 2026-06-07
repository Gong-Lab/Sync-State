function [cells_pDA,n_cells_pDA] = getPutativeDA(spikeTimes, pumpOn, lickT)

fs = 200;
kernelWidth = 0.5;

baselineT=[-600, pumpOn];
[r_b, ~] = spikesToRates(spikeTimes, fs, kernelWidth, baselineT, 'gaussian');

consumT=[lickT(1), lickT(1)+180];
[r_c, ~] = spikesToRates(spikeTimes, fs, kernelWidth, consumT, 'gaussian');

r_b_mean=mean(r_b, 2);
r_b_max=prctile(r_b, 99.7, 2);
r_c_mean=mean(r_c, 2);
r_c_max=prctile(r_c, 99.7, 2);

cells_pDA=find(r_b_mean<=10 & r_b_max<=60 & r_c_mean-r_b_mean>0 & r_c_mean<=15 & r_c_max<=90);%ids of putative DAN
n_cells_pDA=length(cells_pDA);
end