%% this function is used to extract recorded ROIs and remove artifacts from the signals.
function [ResponseTraceAnalysis]=ExtractRemoveArtifact(ResponseTrace,ROIs)

% Define the artifact as the value is more than mean+20std

% get the used ROIs
idx=find(ROIs==1);

% get the names of all the ROIs
ROIn=fieldnames(ResponseTrace);

for i=1:length(idx)
    for channels=1:size(ResponseTrace.(ROIn{idx(i)}),2)
        ChannelLine=find(isoutlier(ResponseTrace.(ROIn{idx(i)})(:,channels),'gesd')==1);
%         if isempty(ChannelLine)
            ResponseTraceAnalysis.(ROIn{idx(i)})(:,channels)=ResponseTrace.(ROIn{idx(i)})(:,channels);
%         else
           for indexChannel=1:length(ChannelLine)
               if ChannelLine(indexChannel)==1
                  ResponseTraceAnalysis.(ROIn{idx(i)})(ChannelLine(indexChannel),channels)=mean(ResponseTrace.(ROIn{idx(i)})(:,channels));
               else
                  ResponseTraceAnalysis.(ROIn{idx(i)})(ChannelLine(indexChannel),channels)=ResponseTraceAnalysis.(ROIn{idx(i)})(ChannelLine(indexChannel)-1,channels);
               end
           end
        end
    end
end
