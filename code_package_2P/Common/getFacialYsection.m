function [y, sce, fs_bhv]=getFacialYsection(ID, tRange2check, t_aroundSce, part)

if nargin<4
    part='Tongue';
end

fDir=['D:\YL\OneDrive\Code\Facial\Facial\', ID, '\SCE_MatchFacial\'];
fName=[fDir, ID, '_facialwithSCEs_raw.mat'];
f=load(fName);
T=f.bpTime_trim-f.actTime_trim(f.lickT(1));
dt_bhv=mean(diff(T));
fs_bhv=1/dt_bhv;

t1i=dsearchn(T', tRange2check(1))+floor(t_aroundSce(1)*fs_bhv);
t2i=dsearchn(T', tRange2check(2))+ceil(t_aroundSce(2)*fs_bhv);
tRange2checki=t1i:t2i;

sce_Time0=f.sce_Time-f.actTime_trim(f.lickT(1));
sce2check=sce_Time0(sce_Time0>=tRange2check(1) & sce_Time0<=tRange2check(2));

switch part
    case 'Tongue'
        y=f.Facial.Tongue(tRange2checki);
        % y=highpass(y, 0.05, fs_bhv);
        y=getFacialNorm(y);
    case 'Jaw'
        y=f.Facial.Jaw(tRange2checki);
        % y=highpass(y, 0.05, fs_bhv);
        y=getFacialNorm(y);
    case 'Nose'
        y=f.Facial.Nose(tRange2checki);
        % y=highpass(y, 0.05, fs_bhv);
        y=getFacialNorm(y);

end


sce=sce2check-T(t1i);

end