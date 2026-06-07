%% 分析self-paced consumption 数据，找到baseline，并在去除baseline之后通过findpeaks分析amplitude, frequency

clc;
clear;
%% setting
MouseID={'vshNAc12'}; % setup1:'vshNAc1','vshNAc2','vshNAc3','vshNAc4','vshNAc5','vshNAc6'
                                  %'BLA9','BLA10','BLA11'  % setup2:'vshNAc9','vshNAc10','vshNAc12'
Date={'20250516'}; %'20250614','20250615','20250616','20250617'% '20250516','20250517','20250518','20250519'
Chamber='setup2'; 
Flavor={'sl'}; %'m','sl'
%saveDir='F:\WRJ\DA sensor\Data\202506\Learning\Combine\Figure\FindPeakTrace\New';
for date=1:length(Date)
    for mouse=1:length(MouseID)
        for flavor=1:length(Flavor)
            stirDir=['F:\WRJ\DA sensor\Data\202505\Learning\' Chamber '\' Date{date} '\' MouseID{mouse} '\' Flavor{flavor} '\processing_DAsensorNew'];
            load([stirDir '\BehaviorEventTime_S.mat']);
            load([stirDir '\ResponseTrace.mat']);
            load([stirDir '\ResponseTraceProcessed.mat']);
            load([stirDir '\BaselineTrace.mat']);
            load([stirDir '\ResponseTraceProcessed_zscore.mat']);
            mkdir([stirDir '\Figure']);
            saveDir=[stirDir '\Figure'];
            %% F0
            ROI=fieldnames(ResponseTraceProcessed);
            fs=20;   % - 采样频率(Hz)
            window_size = 10 * fs;       % 10秒滑动窗口(转换为采样点数)
            for roi=1:length(ROI)
                SensorTrace.(ROI{roi})=ResponseTraceProcessed.(ROI{roi}); 
            end
            
            % 计算z-score 
            begintime=BehaviorEventTime_S{1, 1}.Pump4(1);
            %endtime=BehaviorEventTime_S{1, 1}.Pump4off(end);

            Time1=min(abs(BehaviorEventTime_S{1, 1}.ResponseTime470-begintime));
            SeparationTime1=find(abs(BehaviorEventTime_S{1, 1}.ResponseTime470-begintime)==Time1);
            
            % 用于排除 SeparationTime1 不唯一的bug
            if length(SeparationTime1)>1   
               SeparationTime1 = SeparationTime1(1);
            else
            end
           
            %% find peaks
            for roi=1:length(ROI)
                diff.(ROI{roi})=SensorTrace.(ROI{roi})-BaselineTrace.(ROI{roi});
                dF0 = zeros(size(diff.(ROI{roi})));
                for i = 1:length(diff.(ROI{roi}))
                    % 确定滑动窗口范围
                    win_start = max(1, i - window_size + 1);
                    win_end = i;
                    window_indices = win_start:win_end;
                    if win_start>1
                        % 计算窗口内有效数据的第30百分位数作为F0
                        dF0(i) = prctile(diff.(ROI{roi})(window_indices), 30);
                    else
                        dF0(i) = 0; % 无有效数据时设为0
                    end
                end
                DiffBaseline.(ROI{roi})=dF0;
                %baesline=BaselineTrace.(ROI{roi})(2001:SeparationTime1);     %  removed the first 2000 frames at each session, and defined 2001 frames to pump4-on time as baseline
                baesline=DiffBaseline.(ROI{roi})(2001:SeparationTime1); 
                Mean_baseline=mean(baesline);
                Std_baseline=std(baesline);
                Transient.(ROI{roi}).Trace=(diff.(ROI{roi})(2001:end)-Mean_baseline)/Std_baseline;
                Transient.(ROI{roi}).Time=BehaviorEventTime_S{1, 1}.ResponseTime470(2001:end)-BehaviorEventTime_S{1, 1}.ResponseTime470(2000);
                % 用于排除记录问题导致的数目不相等
                minNum = min(numel(Transient.(ROI{roi}).Trace), numel(Transient.(ROI{roi}).Time));
                Transient.(ROI{roi}).Trace = Transient.(ROI{roi}).Trace(1:minNum);
                Transient.(ROI{roi}).Time = Transient.(ROI{roi}).Time(1:minNum);         
                SmoothTransient.(ROI{roi})=smoothdata(Transient.(ROI{roi}).Trace,'gaussian',7);
                base=SmoothTransient.(ROI{roi})(1:SeparationTime1);
                window_size = 10 * fs;       % 10秒滑动窗口(转换为采样点数)
                F0_low = zeros(size(base));
                for i = 1:length(base)
                    % 确定滑动窗口范围
                    win_start = max(1, i - window_size + 1);
                    win_end = i;
                    window_indices = win_start:win_end;
                    if win_start>1
                        % 计算窗口内有效数据的第25百分位数作为F0
                        F0_low(i) = prctile(base, 25);
                    else
                        F0_low(i) = 0; % 无有效数据时设为0
                    end
                end
                F0_high = zeros(size(base));
                for i = 1:length(base)
                    % 确定滑动窗口范围
                    win_start = max(1, i - window_size + 1);
                    win_end = i;
                    window_indices = win_start:win_end;
                    if win_start>1
                        % 计算窗口内有效数据的第25百分位数作为F0
                        F0_high(i) = prctile(base, 75);
                    else
                        F0_high(i) = 0; % 无有效数据时设为0
                    end
                end
                minH=2.5*(mean(F0_high-F0_low));
                [pks,locs]=findpeaks(SmoothTransient.(ROI{roi}), 'MinPeakDistance',12, 'MinPeakProminence',2.5, 'MinPeakHeight',minH, 'MinPeakWidth',6);
                ResponseTraceProcessed_zscore.(ROI{roi}) = ResponseTraceProcessed_zscore.(ROI{roi})(1:minNum);
                ResponseTraceProcessed_zscore.Time = ResponseTraceProcessed_zscore.Time(1:minNum);
                SmoothResopnse.(ROI{roi})=smoothdata(ResponseTraceProcessed_zscore.(ROI{roi}),'gaussian',7);
                pks=SmoothResopnse.(ROI{roi})(locs);
                % find peaks figure
                f1=figure;
                hold on;
                Time=Transient.(ROI{roi}).Time;
                plot(Time, SmoothResopnse.(ROI{roi})(1:length(Time)));
                %plot(Time, SmoothTransient.(ROI{roi})(1:length(Time)));
                begintime=BehaviorEventTime_S{1, 1}.Pump4(1);
                endtime=BehaviorEventTime_S{1, 1}.Pump4off(end);
                %endtime=BehaviorEventTime_S{1, 1}.pumpAll(346);
                Time1=min(abs(BehaviorEventTime_S{1, 1}.ResponseTime470-begintime));
                Time2=min(abs(BehaviorEventTime_S{1, 1}.ResponseTime470-endtime));
                SeparationTime1=find(abs(BehaviorEventTime_S{1, 1}.ResponseTime470-begintime)==Time1)-2030; % consumption往前多算30 frames
                SeparationTime1=SeparationTime1(1);
                SeparationTime2=find(abs(BehaviorEventTime_S{1, 1}.ResponseTime470-endtime)==Time2)-2000; 
                SeparationTime2=SeparationTime2(1);
                plot([Transient.(ROI{roi}).Time( SeparationTime1) Transient.(ROI{roi}).Time( SeparationTime1)],[min(SmoothTransient.(ROI{roi}))-2 max(SmoothTransient.(ROI{roi}))+2],'--g','Linewidth',1.5); 
                plot([Transient.(ROI{roi}).Time( SeparationTime2) Transient.(ROI{roi}).Time( SeparationTime2)],[min(SmoothTransient.(ROI{roi}))-2 max(SmoothTransient.(ROI{roi}))+2],'--g','Linewidth',1.5); 
                plot(Transient.(ROI{roi}).Time(locs),pks,"om"); % plot out the peaks for manual inspection.
                %yline(minH,'--k','LineWidth',1.5);
                %plot(Transient.(ROI{roi}).Time(locs),pks,'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r','MarkerSize', 10);
                xlabel('Time (s)','fontsize',12,'fontweight','bold');
                ylabel('zscore','fontsize',12,'fontweight','bold');
                %ylim([-5 20]);
                title(sprintf('%s-%s-Trace',string(MouseID{mouse}),string(ROI{roi})),'Fontsize',12,'Fontweight','bold');
                A=[MouseID{mouse} '-' ROI{roi} '-Trace'];
                figPath = fullfile(saveDir, [A '.tif']);
                saveas(gcf, figPath);
                figPath = fullfile(saveDir, [A '.fig']);
                saveas(gcf, figPath);
                %close;

                % 划分phase
                Peaks_base.(ROI{roi}).be.locs=locs(Transient.(ROI{roi}).Time(locs)<Transient.(ROI{roi}).Time(SeparationTime1));
                Peaks_base.(ROI{roi}).be.pks=pks(Transient.(ROI{roi}).Time(locs)<Transient.(ROI{roi}).Time(SeparationTime1));
                Peaks_base.(ROI{roi}).con.locs=locs(Transient.(ROI{roi}).Time(locs)>=Transient.(ROI{roi}).Time(SeparationTime1)-200 & Transient.(ROI{roi}).Time(locs)<=Transient.(ROI{roi}).Time(SeparationTime2));
                Peaks_base.(ROI{roi}).con.pks=pks(Transient.(ROI{roi}).Time(locs)>=Transient.(ROI{roi}).Time(SeparationTime1)-200 & Transient.(ROI{roi}).Time(locs)<=Transient.(ROI{roi}).Time(SeparationTime2));
                Peaks_base.(ROI{roi}).af.loc=locs(Transient.(ROI{roi}).Time(locs)>Transient.(ROI{roi}).Time(SeparationTime2));
                Peaks_base.(ROI{roi}).af.pks=pks(Transient.(ROI{roi}).Time(locs)>Transient.(ROI{roi}).Time(SeparationTime2));
            end    
            
            save([stirDir '\Peaks_base'],"Peaks_base");
            save([stirDir '\Transient'],"Transient");
            

        disp([Date{date} '/' MouseID{mouse} '/' Flavor{flavor}]);
        end
    end
end
