function M = getDragMetric(linklengths, drag_ratio, J_full)
%Calculate drag metric from individual link drag matrices

    %Pull back forces on link into chosen body frame using jacobians
    M = zeros(size(J_full{1},2));
    for idx = 1:numel(J_full)
        L = linklengths(idx);
        d = diag([L,drag_ratio*L,drag_ratio*L^3/12]);
        M = M + J_full{idx}'*d*J_full{idx};
    end

end