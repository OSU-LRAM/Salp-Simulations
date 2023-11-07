function local_force_func = getLocalFixedForceFunction(jet_force_norms, theta, alternate)
%% Get fixed force function for a set of given 

    if nargin == 2
        alternate = true;
    end
    
    % Get force vector 
    sign = 1;
    jet_forces = cell(numel(jet_force_norms), 1);
    for idx = 1:numel(jet_force_norms)
        rot_mat = getRotMat(theta);
        jet_norm = jet_force_norms(idx);
        
        % Alternate sign for each of the forces
        if alternate
            sign = -1 * sign;
        end

        jet_forces{idx + 1} = sign * jet_norm * rot_mat * [1 0 0]';
    end

    % Get total force vector
    jet_forces = vertcat(jet_forces{:});
    local_force_func = @(t, x) local_fixed_force(jet_forces);
end


function rot_mat = getRotMat(theta)
    rot_mat = [
        cos(theta) -sin(theta) 0;
        sin(theta) cos(theta)  0;
        0          0           1];
end


function [local_forces, jet_forces] = local_fixed_force(forces)
    local_forces = forces;  % Total local force at each link
    jet_forces = forces;  % Jet forces for each link
end
