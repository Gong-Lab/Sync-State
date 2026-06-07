function [sce_locs, sce_amps, sce_pCellsID, sce_pCellsID_pk]=findSCEs_New(cellM_zs, MinPeakDist, plotSCEs, compThr) 
% sce_locs: locations of the defined synchronized activities
% sce_amps: mean amplitude of the defined syncrhonized activities across the whole ensemble.
% sce_pCellsID：the index of the cells that participate in the current syncrhonized activities
% sce_pCellsID_pk: amplitude of the peak responses for cells involved in the current syncrhonized activities


if nargin<2
    MinPeakDist=7;
end

if nargin<3
    plotSCEs=false;
end

if nargin<4
    compThr=false;
end


cellM_zs=smoothdata(cellM_zs, 2, 'gaussian',7);

t=[(1:size(cellM_zs,2))/11.4];
tRange=[t(1),t(end)];

[n_cells, n_T]=size(cellM_zs);
pksM=false(size(cellM_zs));
        for c=1:n_cells
            [pks{c}, locs{c}]=findpeaks(cellM_zs(c,:), t, 'MinPeakHeight',2.0, 'MinPeakProminence',1, 'MinPeakDistance',1.0);
            if ~isempty(locs{c})
        pksM(c,knnsearch(t',locs{c}'))=true;
    end
        end
 [sce_locs, ~, ~]=getPE(locs, tRange, 0.4, 0.1, false, 1, false, 'cirShft');
 sce_locs=knnsearch(t',sce_locs);
% actCellCount=getActCellCount_pksM(pksM, wd_size); % bin the calcium transient peak data
 n_sces=length(sce_locs);
% for i=1:n_sces
%     if actCellCount(sce_locs(i))==actCellCount(sce_locs(i)+1)
%         sce_locs(i)=sce_locs(i)+1; % The second plateau point is defined as the synchronized calcium transient event.
%     end
% end

%%
% cellM_mean_shfl=getCellMeanCirShfl_simple(cellM, 500);
% thr_mean=prctile(cellM_mean_shfl, 99, 1);
% cellM_mean=mean(cellM, 1);
% sce_locs=sce_locs(cellM_mean(sce_locs)>thr_mean(sce_locs));
% n_sces=length(sce_locs);

% disp(length(i_rm));
% if ~isempty(i_rm)
%     n_rm=length(i_rm);
%     for i=1:n_rm
%         plot([sce_locs(i_rm(i)), sce_locs(i_rm(i))], [-5, 5+10*(n_cells)], 'm--', 'LineWidth',0.5);
%         % rectangle('Position',[sce_locs(i_rm(i))-wd_hw, -5, wd_size-1, 10*(n_cells+1)], 'EdgeColor','m', 'FaceColor','none', 'LineStyle','--', 'LineWidth',2);
%     end
% end
% figure;
% plot(cellM_mean_shfl', 'color',[0.5, 0.5, 0.5, 0.5]);hold on
% plot(cellM_mean);
% plot(thr_mean);


%% 
sce_amps=nan(1, n_sces);
sce_pCellsID=cell(1, n_sces);
sce_pCellsID_pk=cell(1, n_sces);
cellM_pks=cellM_zs.*pksM;
wd_hw=2;
for i=1:n_sces
    i1=max(sce_locs(i)-wd_hw, 1);
    i2=min(sce_locs(i)+wd_hw, n_T);
    [~, sce_pCellsID{i}, sce_pCellsID_pk{i}]=find(cellM_pks(:, i1:i2)');
    sce_amps(i)=mean(sce_pCellsID_pk{i})*length(sce_pCellsID{i})/n_cells;
    % sce_amps(i)=sum(sce_pCellsID_pk{i});
end

%%
if plotSCEs

 figure; hold on
    offset=0;
    for i=1:n_cells
        c=squeeze(cellM_zs(i,:));
        plot(c+offset);
        plot(knnsearch(t',locs{i}'), pks{i}+offset, 'r.');
        offset=offset+10;
    end

    for i=1:n_sces
%         rectangle('Position',[sce_locs(i)-wd_hw, -5, wd_size-1, 10*(n_cells+1)], 'EdgeColor',[.7 .7 .7], 'FaceColor','none', 'LineStyle','--', 'LineWidth',0.05);
         plot([sce_locs(i), sce_locs(i)], [-5, 5+10*(n_cells)], '-','color',[.7 .7 .7]);
    end
    ylim([-5, 5+10*(n_cells)]);
end
