%% analysis behavior data and fiber photometry --- DA sensor
%% 保存中间文件，

clc;
clear all;

%% the following info must be mannually input.
%MouseID='BLA9'
MouseID={'BLA9','BLA10','BLA11'};%'BLA9','BLA10','BLA11'
Date={'20250614','20250615','20250616','20250617'}; %'20250301','20250302','20250303','20250304'
Flavor={'m','sl'}; %'m','sl'
Setup='setup1';
for mouse=1:length(MouseID)
    for date=1:length(Date)
        for flavor=1:length(Flavor)
            stirDir=['F:\WRJ\DA sensor\Data\202506\Learning\' Setup '\' Date{date} '\' MouseID{mouse} '\' Flavor{flavor} ''];
            %stirDir=['\\10.50.7.153\Ruijie_backup\DA sensor\202505\Cake\vshNAc2'];
            RefChannel=1; % there are three channels in this setup: 405nm, 470nm and 546nm. While 470nm is used for recording, 405nm or 546nm channels
            % can be used as reference channels for movement correction. So, 1 or 3 is
            % required to be input according to the experiment. 1 means 405nm channel
            % is used as a reference, and 3 means 546nm channel is used as the
            % reference.
            ROIs=[0 1 0 1]; % there are four ROIs in this fiber photometry setup (ROI0, ROI1,ROI2,and ROI3). 1 means the current ROI is used, and 0 means
            % the current ROI is not used. Eg,[1 1 0 0] means ROI0 and ROI1 are used in this recording sessions.
            % trial_Duration=5; % Trial duration. The unit is second. Currently, it is usually 5 second.
            %% read the raw data, including both fiber photometry data and behavior recording data

            mkdir([stirDir '\processing_DAsensor']);
            saveDir=[stirDir '\processing_DAsensor'];
            trials_R=1;
            trials_E=1;
            Files=dir( [stirDir '\*.tdms']);

            for FileN=1:length(Files)
                name=[stirDir filesep Files(FileN).name];
                [strDir_sub, strImgFn_prefix, ext] = fileparts(name);
                % read fiber photometry raw data
                if any(strfind(strImgFn_prefix,'Event'))==0
                    ResponseTrace{trials_R}= ParseResponseFile_v3(name);
                    trials_R=trials_R+1;
                else
                    % read behavior event raw data
                    BehaviorEventTime_S{trials_E}= ParseBehaviorEventFile_learning(name);
                    trials_E=trials_E+1;
                end

                clear strDir_sub strImgFn_prefix ext
            end

            save([saveDir '\ResponseTrace'],"ResponseTrace");
            save([saveDir '\BehaviorEventTime_S'],"BehaviorEventTime_S");

            %% Processing fiber photometry data
            ROI_idx=find(ROIs==1);
            for i=1:length(ResponseTrace)
                ROI=fieldnames(ResponseTrace{i});
%                 for ROIn=1:length(ROI_idx)
%                     figure;
%                     plot(ResponseTrace{i}.(ROI{ROI_idx(ROIn)})(:,2)); % plot out raw data in 470nm channels
%                     title([ROI{ROI_idx(ROIn)} 'rawData 470']);
%                     figure;
%                     plot(ResponseTrace{i}.(ROI{ROI_idx(ROIn)})(:,1)); % plot out raw data in 405nm channels
%                     title([ROI{ROI_idx(ROIn)} 'rawData 405']);
%                 end
                [ResponseTraceAnalysis]=ExtractRemoveArtifact(ResponseTrace{i},ROIs); % extract the used ROIs and remove and fill the outliners for used ROIs.
                [ResponseTraceProcessed]=SubstractReferenceSignals(ResponseTraceAnalysis,RefChannel); % substract the reference signals.
            end
            save([saveDir '\ResponseTraceProcessed'],"ResponseTraceProcessed");
            X = [Date{date},' - ',MouseID{mouse},' - ',Flavor{flavor}];
            disp(X)

        end
    end
end
