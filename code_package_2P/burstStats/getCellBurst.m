function [bursts, p, burStat]=getCellBurst(spkt, baselineT, statMethod)
if nargin<3
    statMethod='PS';
end

switch statMethod
    case 'PS'
        [bursts_cand, p_cand]=getPS_allSeeds(spkt, baselineT);
    case 'RGS'
        [bursts_cand, p_cand]=getRGS_allSeeds(spkt, baselineT);
end

if ~isempty(p_cand) && ~all(isnan(p_cand))
    [bursts, p]=filterBur(bursts_cand, p_cand, statMethod);
    burStat=getBurStat(bursts, spkt);
else
    bursts=[];
    p=[];
    burStat=[];
end