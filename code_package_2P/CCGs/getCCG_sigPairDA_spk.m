function [sigPair, pkLocs]=getCCG_sigPairDA_spk(spkt, cells_pDA, tRange2Check, bin_size, max_lag, b_range, p_range)

n_cells_pDA=length(cells_pDA);
ij_combos = nchoosek(1:n_cells_pDA, 2);
n_combos=size(ij_combos,1);

sigPair=nan(n_combos, 2);
pkLocs=nan(1, n_combos);
% for i=1:n_combos
parfor i=1:n_combos
    cellPair=cells_pDA(ij_combos(i, :));
    cell1=cellPair(1);
    cell2=cellPair(2);
    spkt1=spkt{cell1};
    spkt2=spkt{cell2};
    % [ccg_norm, lags, thr_shfl, thr_b, pkloc]=getCCG(event1, event2, tRange2Check, bin_size, max_lag, b_range, p_range, plt, n_b, shflMethod)
    [ccg_norm, lags, thr_shfl, thr_b, pkloc]=getCCG(spkt1, spkt2, tRange2Check, bin_size, max_lag, b_range, p_range, false, 3);
    if ~isempty(pkloc)
        sigPair(i, :)=[cell1, cell2];
        pkLocs(i)=pkloc;
        disp([num2str([i, cell1, cell2]), ' sync']);
        % figure;hold on
        % plot(lags, ccg_norm);
        % plot(lags, thr_shfl.ptwise, 'r--');
        % plot(lags, thr_shfl.global, 'm--');
        % plot([lags(1), lags(end)], [thr_b, thr_b], 'g--');
        % title([num2str([i, cell1, cell2]), ' sync']);
    end
end
idx2exclude=isnan(pkLocs);
sigPair(idx2exclude, :)=[];
pkLocs(idx2exclude)=[];
