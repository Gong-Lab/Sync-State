function spkTimes=getSpkNewTimes2t0(spkTimes_old, t0)
n_cells=length(spkTimes_old);
spkTimes=cell(1, n_cells);
for i=1:n_cells
    spkt=spkTimes_old{i}-t0;
    spk2remove=find(diff(spkt)<0.003)+1;
    spkt(spk2remove)=[];
    spkTimes{i}=spkt;
end