function [NLISI, mu_i]=getNLISI(spkt)
log_isi=log10(diff(spkt));
n_isi=length(log_isi);
Q=max(20, floor(n_isi*0.1));

mu_i=zeros(1, n_isi);

mu_i(1:Q)=getISICentralLocation(log_isi(1:2*Q+1));
for i=Q+1:n_isi-Q
    mu_i(i)=getISICentralLocation(log_isi(i-Q:i+Q));
end
mu_i(n_isi-Q+1:n_isi)=getISICentralLocation(log_isi(n_isi-2*Q:n_isi));

NLISI=log_isi-mu_i;