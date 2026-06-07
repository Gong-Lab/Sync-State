function eventTrain_shfl=shflEventTrain_perm(eventTrain, tRange)
eventTrain=eventTrain(eventTrain>=tRange(1) & eventTrain<=tRange(2));

if size(eventTrain,2)==1
    eventTrain=eventTrain';
end
IEI=[eventTrain(1)-tRange(1), diff(eventTrain), tRange(2)-eventTrain(end)];
permIEI = IEI(randperm(length(IEI))); % random permutation
eventTrain_shfl = tRange(1) + cumsum(permIEI);
eventTrain_shfl(end)=[]; 