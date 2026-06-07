function mu=getISICentralLocation(log_isi)

p=0.05;
Ecenter=(prctile(log_isi, p*100)+prctile(log_isi, (1-p)*100))/2;
eMAD=median(abs(log_isi - median(log_isi)));

centralSet=log_isi(abs(log_isi-Ecenter)<=1.68*eMAD);
cl1=median(centralSet);

centralSet_2=log_isi(abs(log_isi-cl1)<=1.68*eMAD);
mu=median(centralSet_2);