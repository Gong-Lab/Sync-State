function eventTrain_shfl=shflEventTrain_cirShft(eventTrain, tRange)
eventTrain=eventTrain(eventTrain>=tRange(1) & eventTrain<=tRange(2));

eventTrain_shfl=eventTrain+rand*(tRange(2)-tRange(1));
i_overRange=eventTrain_shfl>tRange(2);
eventTrain_shfl(i_overRange)=tRange(1)+(eventTrain_shfl(i_overRange)-tRange(2));

eventTrain_shfl=sort(eventTrain_shfl, 'ascend');
