function [bursts, p]=getRGS_allSeeds(spkt, baselineT)

[NLISI, ~]=getNLISI(spkt);

NLISI_b=NLISI(spkt>=baselineT(1) & spkt<=baselineT(2));
mad_b=median(abs(NLISI_b-median(NLISI_b)));

binEdges=-3:0.1:3;
binCenters=binEdges(1:end-1)+0.1/2;
counts_NLISI=histcounts(NLISI_b, binEdges);
% bar(binCenters, counts_NLISI);
[~, maxCount_i]=max(counts_NLISI);
if abs(binCenters(maxCount_i))>0.2
    thr=0.3;
else
    thr=2.58;
end

centralSet_NLISI_b=NLISI_b(abs(NLISI_b-median(NLISI_b))<=thr*mad_b);
mu=mean(centralSet_NLISI_b);
sigma=std(centralSet_NLISI_b);

% figure;histogram(NLISI_b, 100, 'Normalization','pdf');
% hold on;plot([median(NLISI_b)-thr*mad_b, median(NLISI_b)-thr*mad_b], [0, 1], 'k--');
% plot([median(NLISI_b)+thr*mad_b, median(NLISI_b)+thr*mad_b], [0, 1], 'k--');
% xlim([-3, 3]);
% title('b');
% 
% figure;histogram(NLISI(spkt>=0 & spkt<300), 100, 'Normalization','pdf');
% hold on;plot([median(NLISI_b)-thr*mad_b, median(NLISI_b)-thr*mad_b], [0, 1], 'k--');
% plot([median(NLISI_b)+thr*mad_b, median(NLISI_b)+thr*mad_b], [0, 1], 'k--');
% xlim([-3, 3]);
% title('dur');

i_burSeeds=find((NLISI-median(NLISI_b))<-thr*mad_b);
n_burSeeds=length(i_burSeeds);

n_isi=length(NLISI);

bursts=cell(1,n_burSeeds);
p=nan(1, n_burSeeds);

% parfor i_seed=1:n_burSeeds
for i_seed=1:n_burSeeds
    Tq_ini=NLISI(i_burSeeds(i_seed));
    q_ini=1;
    p_bur_ini=normcdf(Tq_ini, q_ini*mu, q_ini^0.5*sigma);
    %look forward
    i_isi_r_now=i_burSeeds(i_seed);
    q_r=q_ini;
    Tq_r=Tq_ini;
    p_bur_r=p_bur_ini;
    cand_bur_r=[];
    % inBurFlag=false;
    while i_isi_r_now<n_isi
        i_isi_r_now=i_isi_r_now+1;
        Tq_r=Tq_r+NLISI(i_isi_r_now);
        q_r=q_r+1;
        p_bur_temp=normcdf(Tq_r, q_r*mu, q_r^0.5*sigma);
        if p_bur_temp<p_bur_r
            % inBurFlag=true;
            cand_bur_r=[cand_bur_r, i_isi_r_now];
            p_bur_r=p_bur_temp;
        else
            % inBurFlag=false;
            break;
        end
    end

    %look backward
    i_isi_l_now=i_burSeeds(i_seed);
    q_l=q_ini;
    Tq_l=Tq_ini;
    p_bur_l=p_bur_ini;
    cand_bur_l=[];
    % inBurFlag=false;
    while i_isi_l_now>1
        i_isi_l_now=i_isi_l_now-1;
        Tq_l=Tq_l+NLISI(i_isi_l_now);
        q_l=q_l+1;
        p_bur_temp=normcdf(Tq_l, q_l*mu, q_l^0.5*sigma);
        if p_bur_temp<p_bur_l
            % inBurFlag=true;
            cand_bur_l=[i_isi_l_now, cand_bur_l];
            p_bur_l=p_bur_temp;
        else
            % inBurFlag=false;
            break;
        end
    end


    % cand_bur=[cand_bur_l, i_burSeeds(i_seed), cand_bur_r];
    % q=length(cand_bur);
    % Tq=sum(NLISI(cand_bur));
    % p(i_seed)=normcdf(Tq, q*mu, q^0.5*sigma);
    % bursts{i_seed}=[cand_bur, cand_bur(end)+1];

    if ~isempty(cand_bur_l) || ~isempty(cand_bur_r)
        cand_bur=[cand_bur_l, i_burSeeds(i_seed), cand_bur_r];
        q=length(cand_bur);
        Tq=sum(NLISI(cand_bur));
        p(i_seed)=normcdf(Tq, q*mu, q^0.5*sigma);
        bursts{i_seed}=[cand_bur, cand_bur(end)+1];
    else
        bursts{i_seed}=[];
    end

end

% lens = cellfun(@length, bursts);
% idx_del  = find(lens > 20);
% bursts(idx_del)=[];
% p(idx_del)=[];