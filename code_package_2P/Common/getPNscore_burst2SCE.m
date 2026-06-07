function [FN_rate, FP_rate]=getPNscore_burst2SCE(pBurs, sces, tRange, t_lag_thr)
if nargin<4
    t_lag_thr=0.8;
end

pBurs=pBurs(pBurs<=tRange(2) & pBurs>=tRange(1));
sces=sces(sces<=tRange(2) & sces>=tRange(1));

n_pBurs=length(pBurs);
FN=0;
for i=1:n_pBurs
    if ~any(sces-pBurs(i)>0 & sces-pBurs(i)<=t_lag_thr)
        FN=FN+1;
    end
end
FN_rate=FN/n_pBurs;

n_sces=length(sces);
FP=0;
for i=1:n_sces
    if ~any(sces(i)-pBurs>0 & sces(i)-pBurs<=t_lag_thr)
        FP=FP+1;
    end
end
FP_rate=FP/n_sces;