function fr_break=getShflFrBrk(y_shfl, alpha)

pointwise_thr_h=prctile(y_shfl, 100*(1 - alpha), 1);
pointwise_thr_l=prctile(y_shfl, 100*alpha, 1);

fr_break=sum(any(y_shfl-pointwise_thr_h>0 | y_shfl-pointwise_thr_l<0, 2))/size(y_shfl,1);