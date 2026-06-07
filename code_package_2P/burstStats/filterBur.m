function [bursts, p]=filterBur(bursts_cand, p_cand, statMethod)

i_cand=find(~isnan(p_cand));
p_cand=p_cand(i_cand);
n_cand=length(i_cand);

bursts={};
p=nan(1, n_cand);
b=1;
lastSpkSaved=1;

pt=1;
bur_ref=bursts_cand{i_cand(pt)};
p_ref=p_cand(pt);
bur_ref_saveFlag=false;
while pt<n_cand
    pt=pt+1;
    bur_now=bursts_cand{i_cand(pt)};
    p_now=p_cand(pt);
    if bur_now(1)>lastSpkSaved
        spks_overlap=intersect(bur_now, bur_ref);
        if isempty(spks_overlap)
            if ~bur_ref_saveFlag
                bursts{b}=bur_ref;
                p(b)=p_ref;
                lastSpkSaved=bur_ref(end);
                b=b+1;
            end
            bur_ref=bur_now;
            p_ref=p_now;
            bur_ref_saveFlag=false;
        else
            if p_now>p_ref
                bursts{b}=bur_ref;
                p(b)=p_ref;
                lastSpkSaved=bur_ref(end);
                bur_ref_saveFlag=true;
                b=b+1;
            else
                bur_ref=bur_now;
                p_ref=p_now;
                bur_ref_saveFlag=false;
            end
        end
    end
end

p(isnan(p))=[];

%%
n_bur=length(bursts);
iBurStart=zeros(1, n_bur);
for i=1:n_bur
    iBurStart(i)=bursts{i}(1);
end

b2remove=find(diff(iBurStart)<=0)+1;
bursts(b2remove)=[];
p(b2remove)=[];

%% Adjustment for multiple testing
switch(statMethod)
    case 'RGS'
        p=p*length(p)*1.0;

        P_=0.05;
        b2exclude=find(p>=P_);
        if ~isempty(b2exclude)
            p(b2exclude)=[];
            bursts(b2exclude)=[];
        end
    case 'PS'
        b2exclude=find(p>=0.01);
        if ~isempty(b2exclude)
            p(b2exclude)=[];
            bursts(b2exclude)=[];
        end
end
