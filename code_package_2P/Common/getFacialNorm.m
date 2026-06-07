function y_norm=getFacialNorm(y, nrange)
if nargin<2
    nrange='01';
end

switch nrange
    case '01'
        ymax=prctile(y, 99.7);
        ymin=prctile(y, 0.3);
        y_norm=(y-ymin)/(ymax-ymin);
        y_norm(y_norm>1)=1;
        y_norm(y_norm<0)=0;
end
