function [lag, ccg, ccg_cr, results]=getCCGs(spkt1, spkt2, binSize, maxLag, t_range, n_shfl, smth, plt)

% spkt2 ref

% jitter window
jitter_windowSize=0.025;


if nargin<7
    smth=false;
end
if nargin<8
    plt=false;
end

% t_range=[780 880];
% t_range=[100 200];
% binSize=0.001;
% maxLag=0.2;

smoothKernel=[0.05, 0.25, 0.40, 0.25, 0.05];

% [lag, ccg] = cross_correlogram(Target_Index.Spike_AbsTime{51}, Target_Index.Spike_AbsTime{68}, 0.001, 0.5, [650 700]);
[lag, ccg] = ccg_diffHist(spkt1, spkt2, binSize, maxLag, t_range);
if smth
    ccg=conv(ccg, smoothKernel, 'same');
end



% n_shfl=5000;
ccg_shfl=zeros(length(lag), n_shfl);
for i=1:n_shfl
    spkt_i=jitter_spike_times(spkt1, t_range, jitter_windowSize);
    [~, ccg_shfl(:,i)]=ccg_diffHist(spkt_i, spkt2, binSize, maxLag, t_range);

%     spkt_i=jitter_spike_times(spkt1, t_range, 0.025);
%     spkt_j=jitter_spike_times(spkt2, t_range, 0.025);
%     [~, ccg_shfl(:,i)]=ccg_diffHist(spkt_i, spkt_j, binSize, maxLag, t_range);
    if smth
        ccg_shfl(:,i)=conv(ccg_shfl(:,i), smoothKernel, 'same');
    end
end

% plot(lag, ccg_shfl, 'm');
p99=prctile(ccg_shfl, 99, 2);
p1=prctile(ccg_shfl, 1, 2);

ccg_shfl_mean=mean(ccg_shfl, 2)';
ccg_cr=ccg-ccg_shfl_mean;


i_b1=dsearchn(lag', -0.1);
i_b2=dsearchn(lag', -0.05);
i_b3=dsearchn(lag', 0.05);
i_b4=dsearchn(lag', 0.1);
ccg_shfl_b=[ccg_shfl(i_b1:i_b2, :); ccg_shfl(i_b3:i_b4, :)];
ccg_b_mean=mean(mean(ccg_shfl_b, 1));
ccg_b_std=mean(std(ccg_shfl_b, 0, 1));
% ccg_b_std=std(ccg_shfl_b, 0, "all");

results.shflUL=[p99, p1];
results.ccg_b_mean=ccg_b_mean;
results.ccg_b_std=ccg_b_std;

i_p1=dsearchn(lag', -0.01);
i_p2=dsearchn(lag', 0.01);
[ccgPeak, i_ccgPeak]=max(ccg);
% [ccgPeak, i_ccgPeak]=max(ccg_cr);
% if i_ccgPeak>=i_p1 && i_ccgPeak<=i_p2 && ccgPeak>ccg_b_mean+7*ccg_b_std
if i_ccgPeak>=i_p1 && i_ccgPeak<=i_p2 && ccgPeak>p99(i_ccgPeak)
    results.ccg_sigPeak=true;
else
    results.ccg_sigPeak=false;
end

if plt
    figure;hold on
    plot(lag, ccg, 'k');
    plot(lag, p99, 'm-');
    plot(lag, p1, 'm--');
    plot(lag, ccg_cr, 'g');
    plot([lag(1), lag(end)], [ccg_b_mean+7*ccg_b_std, ccg_b_mean+7*ccg_b_std], 'r-');

    xlabel('Lag (s)');
    ylabel('Coincidences per spike');
%     title('Cross-Correlogram');
end