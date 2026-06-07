function fr_break=getShflFrBrk(y_shfl, alpha, side, n_b)

if nargin<3
    side='both';
end
if nargin<4
    n_b=3;
end


switch side
    case 'both'
        pointwise_thr_h=prctile(y_shfl, 100*(1 - alpha), 1);
        pointwise_thr_l=prctile(y_shfl, 100*alpha, 1);
        [~, rowIdx_u]=countRowsWithConsecutiveTrue(y_shfl-pointwise_thr_h>0, n_b);
        [~, rowIdx_d]=countRowsWithConsecutiveTrue(y_shfl-pointwise_thr_l<0, n_b);
        fr_break=sum(rowIdx_u|rowIdx_d, 2)/size(y_shfl,1);
    case 'greater'
        pointwise_thr_h=prctile(y_shfl, 100*(1 - alpha), 1);
        [nRows, ~] = countRowsWithConsecutiveTrue(y_shfl-pointwise_thr_h>0, n_b);
        fr_break=nRows/size(y_shfl,1);
    case 'less'
        pointwise_thr_l=prctile(y_shfl, 100*alpha, 1);
        [nRows, ~] = countRowsWithConsecutiveTrue(y_shfl-pointwise_thr_l<0, n_b);
        fr_break=nRows/size(y_shfl,1);
end