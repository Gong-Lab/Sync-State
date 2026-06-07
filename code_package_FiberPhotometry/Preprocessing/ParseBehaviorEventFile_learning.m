function BehaviorEventTime_S= ParseBehaviorEventFile_learning(name)
         Temp_BehaviorE = tdmsread(name);
         BehaviorEvents_raw=table2array(Temp_BehaviorE{1});
%% lick event time (column 7). 
         BehaviorEventTime_S.Lick=findRisingEdgeOn(BehaviorEvents_raw(:,1));

%% pump event time (including all the pumps,column2). 
         BehaviorEventTime_S.pumpAll=findRisingEdgeOn(BehaviorEvents_raw(:,6));
   
%% fourth pump in channel 4 standby event. The standby duration is 5s(trial duration)
         BehaviorEventTime_S.Pump4=findRisingEdgeOn(BehaviorEvents_raw(:,12));
        BehaviorEventTime_S.Pump4off=findRisingEdgeOff(BehaviorEvents_raw(:,12));
%        BehaviorEventTime_S.Pump1=findRisingEdgeOn(BehaviorEvents_raw(:,9));
%        BehaviorEventTime_S.Pump1off=findRisingEdgeOff(BehaviorEvents_raw(:,9));
 %% Response Time 470nm
        BehaviorEventTime_S.ResponseTime470=findRisingEdgeOn(BehaviorEvents_raw(:,14));
 %% Camera
        BehaviorEventTime_S.Camera=findRisingEdgeOn(BehaviorEvents_raw(:,8));
 %% label
        BehaviorEventTime_S.Label=find(BehaviorEvents_raw(:,16)==1)./1000;
 %% opto BehaviorEvents_raw(:,2)
        