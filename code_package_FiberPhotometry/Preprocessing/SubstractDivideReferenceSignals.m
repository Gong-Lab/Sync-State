%% this function is used for the signal correction based on the fluctuation of the reference channel signals.
function [ResponseTraceProcessed]=SubstractDivideReferenceSignals(ResponseTraceAnalysis,RefChannel)
ROIs=fieldnames(ResponseTraceAnalysis); % get the ROI numbers

for i=1:length(ROIs)
    z=polyfit(ResponseTraceAnalysis.(ROIs{i})(:,RefChannel),ResponseTraceAnalysis.(ROIs{i})(:,2),1); % linear regression to fit reference signals to the 470nm channel signals.
    ResponseTraceAnalysis.(ROIs{i})(:,RefChannel)=polyval(z,ResponseTraceAnalysis.(ROIs{i})(:,RefChannel));
    ResponseTraceProcessed_s.(ROIs{i})=ResponseTraceAnalysis.(ROIs{i})(:,2)-ResponseTraceAnalysis.(ROIs{i})(:,RefChannel);% use recorded 470nm signals to substract fitted reference signals.
    ResponseTraceProcessed.(ROIs{i})=ResponseTraceProcessed_s.(ROIs{i})./ResponseTraceAnalysis.(ROIs{i})(:,RefChannel);
end