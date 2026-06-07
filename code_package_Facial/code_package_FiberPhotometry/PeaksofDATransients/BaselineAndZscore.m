%% 计算整体zscore
% 新方法是先找到traced 30%下分位线，作为baseline, 再计算zscore;
clc;
clear all;

%% the following info must be mannually input.
%MouseID='vshNAc6';
MouseID={'vshNAc9','vshNAc10','vshNAc12'}; % 'vshNAc1','vshNAc2','vshNAc3','vshNAc4','vshNAc5','vshNAc6',%'vshNAc9','vshNAc10','vshNAc12'
Date={'20250516','20250517','20250518','20250519'}; %'20250614','20250615','20250616','20250617' ,'20250516','20250517','20250518','20250519'
Flavor={'m','sl'}; %'m','sl'
Chamber='setup2';
for date=1:length(Date)
    for mouse=1:length(MouseID)
        for flavor=1:length(Flavor)
            stirDir=['F:\WRJ\DA sensor\Data\202505\Learning\' Chamber '\' Date{date} '\' MouseID{mouse} '\' Flavor{flavor} '\processing_DAsensorNew'];
            load([stirDir '\BehaviorEventTime_S.mat']);
            load([stirDir '\ResponseTrace.mat']);
            load([stirDir '\ResponseTraceProcessed.mat']);
            %% 计算baseline 
            ROI=fieldnames(ResponseTraceProcessed);
            fs=20;   % - 采样频率(Hz)
            window_size = 10 * fs;       % 10秒滑动窗口(转换为采样点数)
            for roi=1:length(ROI)
                SensorTrace.(ROI{roi})=ResponseTraceProcessed.(ROI{roi});
                F0 = zeros(size(SensorTrace.(ROI{roi})));
                for i = 1:length(SensorTrace.(ROI{roi}))
                    % 确定滑动窗口范围
                    win_start = max(1, i - window_size + 1);
                    win_end = i;
                    window_indices = win_start:win_end;
                    if win_start>1
                        % 计算窗口内有效数据的第10百分位数作为F0
                        F0(i) = prctile(SensorTrace.(ROI{roi})(window_indices), 10);
                    else
                        F0(i) = 0; % 无有效数据时设为0
                    end
                end
                BaselineTrace.(ROI{roi})=F0;
            end

            %% 计算baseline的z-score
            begintime=BehaviorEventTime_S{1, 1}.Pump4(1);
            endtime=BehaviorEventTime_S{1, 1}.Pump4off(end);
            Time1=min(abs(BehaviorEventTime_S{1, 1}.ResponseTime470-begintime));
            SeparationTime1=find(abs(BehaviorEventTime_S{1, 1}.ResponseTime470-begintime)==Time1);
            
            % 用于排除 SeparationTime1 不唯一的bug
            if length(SeparationTime1)>1   
               SeparationTime1 = SeparationTime1(1);
            else
            end

            for roi=1:length(ROI)
                baesline=BaselineTrace.(ROI{roi})(2001:SeparationTime1);     %  removed the first 2000 frames at each session, and defined 2001 frames to pump4-on time as baseline
                Mean_baseline=mean(baesline);
                Std_baseline=std(baesline);
                BaselineDrift.(ROI{roi}).Trace=(BaselineTrace.(ROI{roi})(2001:end)-Mean_baseline)/Std_baseline;
                BaselineDrift.(ROI{roi}).Time=BehaviorEventTime_S{1, 1}.ResponseTime470(2001:end)-BehaviorEventTime_S{1, 1}.ResponseTime470(2000);
                % 用于排除记录问题导致的数目不相等
                minNum = min(numel(BaselineDrift.(ROI{roi}).Trace), numel(BaselineDrift.(ROI{roi}).Time));
                BaselineDrift.(ROI{roi}).Trace = BaselineDrift.(ROI{roi}).Trace(1:minNum);
                BaselineDrift.(ROI{roi}).Time = BaselineDrift.(ROI{roi}).Time(1:minNum);
            end
            
            save([stirDir '\BaselineTrace'],"BaselineTrace");
            save([stirDir '\BaselineDrift'],"BaselineDrift");
            
            %% 计算整个trace zscore
            for N=1:length(ROI)
                Baesline=BaselineTrace.(ROI{N})(2001:SeparationTime1);     %  removed the first 2000 frames at each session, and defined 2001 frames to pump4-on time as baseline
                Mean_baseline=mean(Baesline);
                Std_baseline=std(Baesline);
                ResponseTraceProcessed_zscore.(ROI{N})=(ResponseTraceProcessed.(ROI{N})(2001:end)-Mean_baseline)/Std_baseline;
            end
            ResponseTraceProcessed_zscore.Time=BehaviorEventTime_S{1, 1}.ResponseTime470(2001:end)-BehaviorEventTime_S{1, 1}.ResponseTime470(2000);
            save([stirDir '\ResponseTraceProcessed_zscore'],"ResponseTraceProcessed_zscore");

            clearvars -except MouseID mouse Date date Flavor flavor StirDir name  ResponseTraceProcessed BehaviorEventTime_S Chamber
            disp([Date{date} '/' MouseID{mouse} '/' Flavor{flavor}]);
        end
    end
end

            





