function MI = mutual_info_binned(binned1, binned2)
    % Ensure logical
    x = logical(binned1(:));
    y = logical(binned2(:));
    N = length(x);

    % Joint probabilities
    p11 = sum(x & y)/N;
    p10 = sum(x & ~y)/N;
    p01 = sum(~x & y)/N;
    p00 = sum(~x & ~y)/N;

    % Marginals
    px1 = p11 + p10;
    px0 = p01 + p00;
    py1 = p11 + p01;
    py0 = p10 + p00;

    % Flatten joint and corresponding marginals
    joint = [p11 p10 p01 p00];
    px = [px1 px1 px0 px0];
    py = [py1 py0 py1 py0];

    % Compute MI, ignore zero probabilities
    idx = joint > 0;
    MI = sum(joint(idx) .* log2(joint(idx) ./ (px(idx).*py(idx))));
end




