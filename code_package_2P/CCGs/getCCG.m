function [ccg_norm, lags, thr_shfl, thr_b, pkloc]=getCCG(event1, event2, tRange2Check, bin_size, max_lag, b_range, p_range, plt, n_b, shflMethod, parRun)

if nargin<8
    plt=true;
end
if nargin<9
    n_b=3;
end
if nargin<10
    shflMethod='cirShft';
end
if nargin<11
    parRun=false;
end


event1=event1(event1>=tRange2Check(1) & event1<tRange2Check(2));
event2=event2(event2>=tRange2Check(1) & event2<tRange2Check(2));

if isempty(event1) || isempty(event2)
    ccg_norm=[];
    lags=[];
    thr_shfl=[];
    thr_b=[]; 
    pkloc=[];
    return;
end

%%
[lags, ccg_norm] = ccg_diffHist(event1, event2, bin_size, max_lag, tRange2Check); %event2-event1

n_shfl=10000;
ccg_norm_shfl=nan(n_shfl, length(lags));
if parRun
    parfor s=1:n_shfl
        switch shflMethod
            case 'perm'
                [~, ccg_temp]=ccg_diffHist(event1, shflEventTrain_perm(event2, tRange2Check), bin_size, max_lag, tRange2Check);
            case 'cirShft'
                [~, ccg_temp]=ccg_diffHist(event1, shflEventTrain_cirShft(event2, tRange2Check), bin_size, max_lag, tRange2Check);
            case 'jitterSpks'
                jitter_windowSize=0.025;
                [~, ccg_temp]=ccg_diffHist(event1, jitter_spike_times(event2, tRange2Check, jitter_windowSize), bin_size, max_lag, tRange2Check);

        end
        ccg_norm_shfl(s, :)=ccg_temp;
    end
else
    for s=1:n_shfl
        switch shflMethod
            case 'perm'
                [~, ccg_temp]=ccg_diffHist(event1, shflEventTrain_perm(event2, tRange2Check), bin_size, max_lag, tRange2Check);
            case 'cirShft'
                [~, ccg_temp]=ccg_diffHist(event1, shflEventTrain_cirShft(event2, tRange2Check), bin_size, max_lag, tRange2Check);
            case 'jitterSpks'
                jitter_windowSize=0.025;
                [~, ccg_temp]=ccg_diffHist(event1, jitter_spike_times(event2, tRange2Check, jitter_windowSize), bin_size, max_lag, tRange2Check);

        end
        ccg_norm_shfl(s, :)=ccg_temp;
    end
end

if strcmp(shflMethod, 'jitterSpks')
    ccg_norm=ccg_norm-mean(ccg_norm_shfl, 1);
end

%%
alpha=0.05;
thr_shfl.ptwise=prctile(ccg_norm_shfl, (1-alpha)*100, 1);

target_globalAlpha=0.05;
alpha_adj=findAlphaBelowThreshold(ccg_norm_shfl, target_globalAlpha, 1/n_shfl, 0.05, [], 'greater', n_b);
% alpha_adj=0.05;
thr_shfl.global=prctile(ccg_norm_shfl, 100*(1 - alpha_adj), 1);

%%
i_b1=dsearchn(lags', -b_range(2));
i_b2=dsearchn(lags', -b_range(1));
i_b3=dsearchn(lags', b_range(1));
i_b4=dsearchn(lags', b_range(2));
ccg_b=[ccg_norm(i_b1:i_b2), ccg_norm(i_b3:i_b4)];
thr_b=mean(ccg_b)+5*std(ccg_b,0);

i_p1=dsearchn(lags', -p_range);
i_p2=dsearchn(lags', p_range);
[ccgPeak, i_ccgPeak]=max(ccg_norm);

if i_ccgPeak>i_p1 && i_ccgPeak<i_p2 && ccgPeak>thr_b && ccgPeak>thr_shfl.global(i_ccgPeak)
    pkloc=lags(i_ccgPeak);
else
    pkloc=[];
end

if plt
    figure;hold on
    plot(lags, ccg_norm);
    plot(lags, thr_shfl.ptwise, 'r--');
    plot(lags, thr_shfl.global, 'm--');
    plot([lags(1), lags(end)], [thr_b, thr_b], 'g--');
    if ~isempty(pkloc)
        plot(pkloc, ccgPeak, 'rv');
    end

    plot([-b_range(2), -b_range(2)], [0, max(ccg_norm)], 'k--');
    plot([-b_range(1), -b_range(1)], [0, max(ccg_norm)], 'k--');
    plot([b_range(1), b_range(1)], [0, max(ccg_norm)], 'k--');
    plot([b_range(2), b_range(2)], [0, max(ccg_norm)], 'k--');

    plot([-p_range, -p_range], [0, max(ccg_norm)], 'k--');
    plot([p_range, p_range], [0, max(ccg_norm)], 'k--');
    if isempty(pkloc)
        title('no significant peak');
    else
        title(['peak at ', num2str(pkloc, 3)]);
    end
end