function [bursts, p]=getPS_allSeeds(spkt, baselineT)
spkt_b=spkt(spkt>=baselineT(1) & spkt<baselineT(2));
isi_b_mean=mean(diff(spkt_b));
r_b=1/isi_b_mean;

isi=diff(spkt);

isi_thr_seed=0.5;

i_burSeeds=find(isi<isi_b_mean*isi_thr_seed);
n_burSeeds=length(i_burSeeds);

n_isi=length(isi);

bursts=cell(1,n_burSeeds);
p=nan(1, n_burSeeds);

for i_seed=1:n_burSeeds
    k_ini=2;
    Tk_ini=isi(i_burSeeds(i_seed));
    p_bur_ini=exp(-r_b*Tk_ini)*(r_b*Tk_ini)^k_ini/factorial(k_ini);
    
    %look forward
    i_isi_r_now=i_burSeeds(i_seed);
    k_r=k_ini;
    Tk_r=Tk_ini;
    p_bur_r=p_bur_ini;
    cand_bur_r=[];
    % inBurFlag=false;
    while i_isi_r_now<n_isi
        i_isi_r_now=i_isi_r_now+1;
        k_r=k_r+1;
        dt=isi(i_isi_r_now);
        Tk_r=Tk_r+dt;
        p_bur_temp=exp(-r_b*Tk_r)*(r_b*Tk_r)^k_r/factorial(k_r);
        if dt<=isi_b_mean*isi_thr_seed && p_bur_temp<p_bur_r
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
    k_l=k_ini;
    Tk_l=Tk_ini;
    p_bur_l=p_bur_ini;
    cand_bur_l=[];
    % inBurFlag=false;
    while i_isi_l_now>1
        i_isi_l_now=i_isi_l_now-1;
        k_l=k_l+1;
        dt=isi(i_isi_l_now);
        Tk_l=Tk_l+dt;
        p_bur_temp=exp(-r_b*Tk_l)*(r_b*Tk_l)^k_l/factorial(k_l);
        if dt<isi_b_mean*isi_thr_seed && p_bur_temp<p_bur_l
            % inBurFlag=true;
            cand_bur_l=[i_isi_l_now, cand_bur_l];
            p_bur_l=p_bur_temp;
        else
            % inBurFlag=false;
            break;
        end
    end


    % cand_bur=[cand_bur_l, i_burSeeds(i_seed), cand_bur_r];
    % k=length(cand_bur)+1;
    % Tk=sum(isi(cand_bur));
    % p(i_seed)=exp(-r_b*Tk)*(r_b*Tk)^k/factorial(k);
    % bursts{i_seed}=[cand_bur, cand_bur(end)+1];

    if ~isempty(cand_bur_l) || ~isempty(cand_bur_r)
        cand_bur=[cand_bur_l, i_burSeeds(i_seed), cand_bur_r];
        k=length(cand_bur)+1;
        Tk=sum(isi(cand_bur));
        p(i_seed)=exp(-r_b*Tk)*(r_b*Tk)^k/factorial(k);
        bursts{i_seed}=[cand_bur, cand_bur(end)+1];
    else
        bursts{i_seed}=[];
    end

end

lens = cellfun(@length, bursts);
idx_del  = find(lens >= 20);
bursts(idx_del)=[];
p(idx_del)=[];