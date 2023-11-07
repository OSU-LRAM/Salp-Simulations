function J_dual_full = getDualJacobians_discrete(J_full)
% Get distribution matrix for a set of forces
%
% Inputs:
%
%   J_full: Cell array with full Jacobians from body velocity of the base
%           frame and joint angle velocities to body velocity of each link
%           relative to a non-moving frame.
%
% Outputs:
%
%   B: Dual Jacobians for 

    J_dual_full = cellfun(@transpose, J_full, 'UniformOutput', false);
    J_dual_full = [J_dual_full{:}];

end