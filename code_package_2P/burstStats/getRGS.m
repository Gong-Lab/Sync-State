function [bursts, pauses]=getRGS(spkt, baselineT)

[NLISI, ~]=getNLISI(spkt);

NLISI_b=NLISI(spkt>=baselineT(1) & spkt<=baselineT(2));
mad_b=median(abs(NLISI_b-median(NLISI_b)));

thr=2.58;
centralSet_NLISI_b=NLISI_b(abs(NLISI_b)<=thr*mad_b);
mu=mean(centralSet_NLISI_b);
sigma=std(centralSet_NLISI_b);

i_burSeeds=find(NLISI<-thr*mad_b);
n_burSeeds=length(i_burSeeds);

bursts={};
i_seed_now=1;
inBurFlag=false;
i_b=1;

n_isi=length(NLISI);
i_isi_now=0;

while i_seed_now<=n_burSeeds && i_isi_now<=n_isi
%     if mod(i_isi_now, 100)==0
%         disp([num2str(i_isi_now), ' of ', num2str(n_isi)]);
%     end
    if ~inBurFlag
        q=1;
        Tq_ini=NLISI(i_burSeeds(i_seed_now));
        p_bur_ini=normcdf(Tq_ini, q*mu, q^0.5*sigma);
        i_isi_now=i_burSeeds(i_seed_now)+1;

        if i_isi_now>n_isi
            break;
        else
            q=q+1;
            Tq=Tq_ini+NLISI(i_isi_now);
            p_bur=normcdf(Tq, q*mu, q^0.5*sigma);

            if p_bur<p_bur_ini
                inBurFlag=true;
                cand_bur=[i_burSeeds(i_seed_now), i_isi_now];
                i_isi_now=i_isi_now+1;
            else
                i_seed_now=i_seed_now+1;
            end
        end
    else
        q=q+1;
        Tq=Tq+NLISI(i_isi_now);
        p_bur_temp=normcdf(Tq, q*mu, q^0.5*sigma);
        if p_bur_temp<p_bur
            p_bur=p_bur_temp;
            cand_bur=[cand_bur, i_isi_now];
            i_isi_now=i_isi_now+1;
        else
            inBurFlag=false;
            if q-1>=1
                bursts{i_b}=[cand_bur, cand_bur(end)+1];
                i_b=i_b+1;
            else
                cand_bur=[];
            end
            i_seed_now=find(i_burSeeds>=i_isi_now, 1);
            if isempty(i_seed_now)
                break;
            end
        end
    end    
end

pauses={};