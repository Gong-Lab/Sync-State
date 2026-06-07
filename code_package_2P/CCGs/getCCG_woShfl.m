function [ccg_norm, lags, ccgPeak, pkloc]=getCCG_woShfl(event1, event2, tRange2Check, bin_size, max_lag, p_range, plt)

if nargin<7
    plt=false;
end


event1=event1(event1>=tRange2Check(1) & event1<tRange2Check(2));
event2=event2(event2>=tRange2Check(1) & event2<tRange2Check(2));

if isempty(event1) || isempty(event2)
    ccg_norm=[];
    lags=[];
    ccgPeak=[]; 
    pkloc=[];
    return;
end

%%
[lags, ccg_norm] = ccg_diffHist(event1, event2, bin_size, max_lag, tRange2Check); %event2-event1

%%

i_p1=dsearchn(lags', -p_range);
i_p2=dsearchn(lags', p_range);
i_pkRange=i_p1:i_p2;
[ccgPeak, i_ccgPeak_temp]=max(ccg_norm(i_pkRange));

pkloc=lags(i_pkRange(i_ccgPeak_temp));

if plt
    figure;hold on
    plot(lags, ccg_norm);
    plot(pkloc, ccgPeak, 'rv');
    plot([-p_range, -p_range], [0, max(ccg_norm)], 'k--');
    plot([p_range, p_range], [0, max(ccg_norm)], 'k--');
end