function [lickT, pumpOn, t0]=getLickT(Target_Index, LickFile)

pumpOn=Target_Index.Trial_Start(11);
LickEvent=LickFile.LickEvent;
lickT=pumpOn+LickEvent;

t0=lickT(1);

pumpOn=pumpOn-t0;
lickT=lickT-t0;
lick2remove=find(diff(lickT)<0.1)+1;
lickT(lick2remove)=[];