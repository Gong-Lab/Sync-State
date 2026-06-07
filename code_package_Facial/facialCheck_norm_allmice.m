clear; 
clc; 
close all; 
MouseID = {'Urey', 'Young', 'Zwicky', 'Biscuit','Cinnamon','Fubao','HaHa','LaLa'};
FacialFeature = {'Tongue', 'Nose', 'Pupil', 'WhiskerMovement', 'Jaw'};
Path = '\\10.50.7.153\Functional_Imaging\Conditioning_consumption_task\Facial\';
FacialID = 4;

ff = [];
sce2check_all = [];
lickT2check_all = [];
lickBout2check_all = [];
y_eventTrig_all = [];
fs_bhv_all = [];
T_all = [];

for id = 1: length(MouseID)
    id
    load([Path MouseID{id} '\SCE_MatchFacial\' MouseID{id} '_facialwithSCEs_raw.mat']);
    
    tRange2check=[0, 200]; % s
    T=bpTime_trim-actTime_trim(lickT(1));
    t1i=dsearchn(T', tRange2check(1));
    t2i=dsearchn(T', tRange2check(2));
    tRange2checki=t1i:t2i;
    T_single = T(t1i:t2i) + (id-1)*200;
    T_all = [T_all, T_single];

    sce_Time0=sce_Time-actTime_trim(lickT(1));
    sce2check=sce_Time0(sce_Time0>=tRange2check(1) & sce_Time0<=tRange2check(2));
    sce2check_single=sce2check + (id-1)*200;
    sce2check_all = [sce2check_all, sce2check_single];
    %% Choose Facial Feature !!!
    y_raw=Facial.WM; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Normalization
    y = y_raw/prctile(abs(y_raw(:)),99.7);
    %% Remove outliner
    idx_olh = find(y>1);
    y(idx_olh) = 1;
    idx_oll = find(y<-1);
    y(idx_oll) = -1;

    ff = [ff, y(tRange2checki)];

    lickT0=actTime_trim(lickT)-actTime_trim(lickT(1));
    lickT2check=lickT0(lickT0>=tRange2check(1) & lickT0<=tRange2check(2));
    lickT2check_single = lickT2check + (id-1)*200;
    lickT2check_all = [lickT2check_all, lickT2check_single];
    
    dt_bhv=mean(diff(T));
    fs_bhv=1/dt_bhv;
    fs_bhv_all = [fs_bhv_all, fs_bhv];
    %% Plot Raw Feature Trace & SCE & Lick event
    figure;hold on
    plot(T(tRange2checki), y(tRange2checki));
    plot(sce2check, 40, 'm|', 'MarkerSize',50, 'LineWidth',1);
    plot(lickT2check, 45, 'r|', 'MarkerSize',20, 'LineWidth',1);
    
    lickBout_indices=getBurst(lickT0, [0.2, 0.4]);
    n_lickBout=length(lickBout_indices);
    lickBout_start=zeros(1, n_lickBout);
    for i=1:n_lickBout
        lickBout_start(i)=lickT0(lickBout_indices{i}(1));
    end
    lickBout2check=lickBout_start(lickBout_start>=tRange2check(1) & lickBout_start<=tRange2check(2));
    lickBout2check_single = lickBout2check + (id-1)*200;
    lickBout2check_all = [lickBout2check_all, lickBout2check_single];
    plot(lickBout2check-0.01, 45, 'c|', 'MarkerSize',30, 'LineWidth',1);
    
    %% Sanity Check (set putative sce as lick bout start)
    % sce2check=lickBout2check-0.01;
    %%
    % figure;hold on
    % plot(Facial.Jaw(tRange2checki), Facial.Tongue(tRange2checki), 'k.');
    % [r, p]=corr(Facial.Jaw(tRange2checki)', Facial.Tongue(tRange2checki)');
    % xlabel('Jaw');
    % ylabel('Tongue');
    % title(['r=', num2str(r, 3), ' p=', num2str(p, 3)]);
    % 
    % figure;hold on
    % plot(Facial.Jaw(tRange2checki), Facial.Nose(tRange2checki), 'k.');
    % [r, p]=corr(Facial.Jaw(tRange2checki)', Facial.Nose(tRange2checki)');
    % xlabel('Jaw');
    % ylabel('Nose');
    % title(['r=', num2str(r, 3), ' p=', num2str(p, 3)]);
    % 
    % figure;hold on
    % plot(Facial.Jaw(tRange2checki), Facial.Pupil(tRange2checki), 'k.');
    % [r, p]=corr(Facial.Jaw(tRange2checki)', Facial.Pupil(tRange2checki)');
    % xlabel('Jaw');
    % ylabel('Pupil');
    % title(['r=', num2str(r, 3), ' p=', num2str(p, 3)]);
    % 
    % figure;hold on
    % plot(Facial.Jaw(tRange2checki), Facial.WM(tRange2checki), 'k.');
    % [r, p]=corr(Facial.Jaw(tRange2checki)', Facial.WM(tRange2checki)');
    % xlabel('Jaw');
    % ylabel('WM');
    % title(['r=', num2str(r, 3), ' p=', num2str(p, 3)]);
    
    %%
    % figure;
    % plot(Facial.Nose_x(tRange2checki), Facial.Nose_y(tRange2checki), 'k.');
    % xlabel('Nose_x');
    % ylabel('Nose_y');
    
    %%
    % y=highpass(y, 0.05, fs_bhv);
    % y=lowpass(y, 0.7, fs);
    % figure;
    % plot(T(tRange2checki), y(tRange2checki));
     %% Plot Spectrum
%     y_=y(tRange2checki);
%     [s,f,t]=spectrogram(y_,round(fs_bhv*2), round(fs_bhv*1.5), 4096*16,fs_bhv,'yaxis');
%     t=t+T(tRange2checki(1));
%     m=10*log10(abs(s));
%     
%     fmin=dsearchn(f, 0.1);
%     fmax=dsearchn(f, 20);
%     
%     figure;
%     % imagesc(t, f, m);
%     imagesc(t, f(fmin:fmax), m(fmin:fmax, :));
%     axis xy
    
    %% Mean Trace
    n_sce2check=length(sce2check);
    t_aroundSce=[-1, 1];
    ti_aroundSce=round(t_aroundSce(1)*100):round(t_aroundSce(2)*100);
    y_eventTrig=zeros(n_sce2check, length(ti_aroundSce));
    for i=1:n_sce2check
        i1=dsearchn(T', sce2check(i))+ti_aroundSce(1);
        i2=dsearchn(T', sce2check(i))+ti_aroundSce(end);
        % y_eventTrig(i,:)=y_hp(i1:i2);
        y_eventTrig(i,:)=y(i1:i2);
    end
    y_eventTrig_mean=mean(y_eventTrig, 1);
    y_eventTrig_all = [y_eventTrig_all; y_eventTrig];
    figure(id*3);hold on
    plot(ti_aroundSce/fs_bhv, y_eventTrig', 'LineStyle','--', 'Color',[0.3,0.3,0.3]);
    plot(ti_aroundSce/fs_bhv, y_eventTrig_mean, 'k-', 'LineWidth',2);
    
    %% 
    n_shfl=10000;
    y_eventTrig_mean_shfl=zeros(n_shfl, length(ti_aroundSce));
    parfor s=1:n_shfl
        sce2check_shfl=shflEventTrain_cirShft(sce2check, tRange2check); % Shuffle Method 1
    %     sce2check_shfl=shflEventTrain_perm(sce2check, tRange2check); % Shuffle Method 2
        y_eventTrig_shfl=zeros(n_sce2check, length(ti_aroundSce));
        for i=1:n_sce2check
            i1=dsearchn(T', sce2check_shfl(i))+ti_aroundSce(1);
            i2=dsearchn(T', sce2check_shfl(i))+ti_aroundSce(end);
            % y_eventTrig(i,:)=y_hp(i1:i2);
            y_eventTrig_shfl(i,:)=y(i1:i2);
        end
        y_eventTrig_mean_shfl(s,:)=mean(y_eventTrig_shfl, 1);
        
    end
%     y_eventTrig_mean_shfl_all = [y_eventTrig_mean_shfl_all, y_eventTrig_mean_shfl];
    %% Shuffled Trace Thr
    target_globalAlpha=0.05;
    alpha_adj=findAlphaBelowThreshold(y_eventTrig_mean_shfl, target_globalAlpha, 1/n_shfl, 0.05);
    thr_h=prctile(y_eventTrig_mean_shfl, 100*(1 - alpha_adj), 1);
    thr_l=prctile(y_eventTrig_mean_shfl, 100*alpha_adj, 1);
    
    
    % plot(ti_aroundSce/fs_bhv, y_eventTrig_mean_shfl(1:20, :)', 'b--', 'LineWidth',1);
    plot(ti_aroundSce/fs_bhv, thr_h, 'r--', 'LineWidth',1.5);
    plot(ti_aroundSce/fs_bhv, thr_l, 'r--', 'LineWidth',1.5);
    
    title([MouseID(id) FacialFeature{FacialID} '(' num2str(tRange2check) ' s)']);
%     savefig(figure(id*3),['\\10.50.7.153\Functional_Imaging\Conditioning_consumption_task\Facial\FacialFeature_alignSCE\' MouseID{id} '_' FacialFeature{FacialID} '.fig'])
end

%% Total Shuffle

    y_eventTrig_mean_all=mean(y_eventTrig_all, 1);
    t_aroundSce=[-1, 1];
    ti_aroundSce=round(t_aroundSce(1)*100):round(t_aroundSce(2)*100);
    figure;hold on
    plot(ti_aroundSce/100, y_eventTrig_all', 'LineStyle','--', 'Color',[0.3,0.3,0.3]);
    plot(ti_aroundSce/100, y_eventTrig_mean_all, 'k-', 'LineWidth',2);

    n_shfl=10000;
    tRange2check_all = [0 length(MouseID)*200];

    n_sce2check_all=length(sce2check_all);
    y_eventTrig_mean_shfl_all=zeros(n_shfl, 201);
    parfor s=1:n_shfl
        sce2check_shfl=shflEventTrain_cirShft(sce2check_all, tRange2check_all); % Shuffle Method 1
    %     sce2check_shfl=shflEventTrain_perm(sce2check, tRange2check); % Shuffle Method 2
        y_eventTrig_shfl_all=zeros(n_sce2check_all, 201);
        for i=1:n_sce2check_all
            i1=dsearchn(T_all', sce2check_shfl(i))+ti_aroundSce(1);
            i2=dsearchn(T_all', sce2check_shfl(i))+ti_aroundSce(end);
            if i1 < 1 || i2 > length(ff)
                y_eventTrig_shfl_all(i,:)=nan;
            else
                % y_eventTrig(i,:)=y_hp(i1:i2);
                y_eventTrig_shfl_all(i,:)=ff(i1:i2);
            end
        end
        y_eventTrig_mean_shfl_all(s,:) = nanmean(y_eventTrig_shfl_all, 1);
    end
%     y_eventTrig_mean_shfl_all = [y_eventTrig_mean_shfl_all, y_eventTrig_mean_shfl];
    %% Shuffled Trace Thr
    target_globalAlpha=0.05;
    alpha_adj=findAlphaBelowThreshold(y_eventTrig_mean_shfl_all, target_globalAlpha, 1/n_shfl, 0.05);
    thr_h=prctile(y_eventTrig_mean_shfl_all, 100*(1 - alpha_adj), 1);
    thr_l=prctile(y_eventTrig_mean_shfl_all, 100*alpha_adj, 1);
    
    
    % plot(ti_aroundSce/fs_bhv, y_eventTrig_mean_shfl(1:20, :)', 'b--', 'LineWidth',1);
    plot(ti_aroundSce/100, thr_h, 'r--', 'LineWidth',1.5); hold on;
    plot(ti_aroundSce/100, thr_l, 'r--', 'LineWidth',1.5);
    
    title([FacialFeature{FacialID}]);