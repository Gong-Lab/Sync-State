function burStat=getBurStat(bursts, spkt)

n_bur=length(bursts);
burStartT=nan(1, n_bur);
burEndT=nan(1, n_bur);
burMedT=nan(1, n_bur);
burSpkNum=nan(1, n_bur);

for i=1:n_bur
    burStartT(i)=spkt(bursts{i}(1));
    burEndT(i)=spkt(bursts{i}(end));
    burMedT(i)=median(spkt(bursts{i}));
    burSpkNum(i)=length(bursts{i});
end

burStat.burStartT=burStartT;
burStat.burEndT=burEndT;
burStat.burMedT=burMedT;
burStat.spkNum=burSpkNum;
