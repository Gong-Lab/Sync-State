function ResponseTrace= ParseResponseFile_v3(name)
       Temp_CalciumR = tdmsread(name); % read tdms files.
       ResponseRaw= table2array(Temp_CalciumR{1}); % convert the raw data into a matrix

%% check if any channel is missing during recording. Original channels are 407nm, 470nm and 546nm. 
       numberOfZeros = sum(ResponseRaw(:,1)==0); % column 1-4 are 405nm channels. 
       channelJudge=(numberOfZeros==length(ResponseRaw));
       if channelJudge==1
           channels=0;
           disp('405nm is missing during recording');
       else 
       numberOfZeros = sum(ResponseRaw(:,5)==0); % column 5-8 are 470nm channels. 
       channelJudge=(numberOfZeros==length(ResponseRaw));
       if channelJudge==1
          channels=1;
          disp('470nm is missing during recording'); 
       else
       numberOfZeros = sum(ResponseRaw(:,9)==0); % column 9-12 are 546nm channels. 
       channelJudge=(numberOfZeros==length(ResponseRaw));
       if channelJudge==1
          channels=2;
          disp('546nm is missing during recording');   
       else
          channels=3;
       end
       end
       end
       %% get the raw data trace for each channels
%        switch channels
%            case 0
%                for i=5:12 % get columns 5-12
%                    array(:,i)=nonzeros(ResponseRaw(:,i));
%                end
%                    array=[zeros(size(array,1),4) array];
%            case 1
%                for i=1:4 % get columns 1-4 and columns 9-12
%                    array(:,i)=nonzeros(ResponseRaw(:,i));
%                    array(:,i+8)=nonzeros(ResponseRaw(:,i+8));
%                end
%            case 2
%                for i=[1 2] % get columns 1,2,5,6
%                    array(:,i)=nonzeros(ResponseRaw(:,i));
%                end
%                for i=[5 6] % get columns 1,2,5,6
%                    temp=nonzeros(ResponseRaw(:,i));
%                    if length(temp)<size(array,1)
%                        array(:,i)=[temp;zeros(size(array,1)-length(temp),1)];
%                    else
%                        array(:,i)=temp(1:size(array,1));
%                    end                   
%                end
%                for i=[7 8 9 10 11 12] % get columns 1,2,5,6
%                    array(:,i)=zeros(size(array,1),1);
%                end
%            case 3
%                array(:,1)=nonzeros(ResponseRaw(:,1));
%                for i=2:12 % get columns 1-12
%                    Temp=nonzeros(ResponseRaw(:,i));
%                    if length(Temp)>=size(array,1)
%                       array(:,i)=Temp(1:size(array,1),:);
%                    else
%                       array(:,i)=[Temp;0];
%                    end
%                end
%            otherwise
%        end

       array=[];
for i=1:12
    for j=1:length(ResponseRaw(:,i))/2
        z=[ResponseRaw((j-1)*2+1,i);ResponseRaw((j-1)*2+2,i)];
        array(j,i)=max(z);
    end
end
 %% divide the data according ROIs . the fiber photometry setup in the RongLab currently contains four ROIs (four fibers), which could simutaneously record four brain regions.
    ResponseTrace.ROI_0=array(:,[1 5 9]);
    ResponseTrace.ROI_1=array(:,[2 6 10]);
    ResponseTrace.ROI_2=array(:,[3 7 11]);
    ResponseTrace.ROI_3=array(:,[4 8 12]);
    ResponseTrace.Time=ResponseRaw(2:2:end,13);

end
