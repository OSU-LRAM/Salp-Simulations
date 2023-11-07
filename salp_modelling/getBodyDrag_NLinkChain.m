function M_full = getBodyDrag_NLinkChain(J_full, h, drag_ratio)
% Calculate the matrix that maps from body forces to forces acting on the
% base frame of the system

    M_full = zeros(size(J_full{1}));  % Storage for metric tensor
    
    for idx = 1:numel(J_full)
    
        % Create drag matrix for current link
        L = h.lengths(idx);  % Link length for current link
        drag_mat = diag([L, drag_ratio*L, drag_ratio/12*L^3]);    
        M_full = M_full +  transpose(J_full{idx}) * drag_mat * J_full{idx};
    
    end

end

