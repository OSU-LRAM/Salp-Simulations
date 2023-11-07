function local_force_func = getMinForceFunction(n, B_aug, domain, q_dot_des, thrust_angle, omega, bias_constant_jets)
    %% Output force function for the jets 

    % TODO: Maybe add more stuff to precalculate

    base_field_num = size(B_aug{1}, 1) - 2;  % Get the number of independent vector fields 

    %% Get big rotation matrix that maps thrust input thrust vector to forces
    R = zeros([3*n, n]);
    for idx = 1:n
        if mod(idx, 2) == 0
            angle = thrust_angle;
        else
            angle = -thrust_angle;
        end
        R(3*(idx-1)+1: 3*(idx-1)+3, idx) = [cos(angle); sin(angle); 0];
    end

    local_force_func = @(t, q) min_jet_force(t, q, B_aug, domain, q_dot_des, omega, R, base_field_num, bias_constant_jets);
end


function [local_forces, jet_forces] = min_jet_force(t, q, B_aug, domain, q_dot_des, omega, R, base_field_num, bias_constant_jets)
    %% Generate ideal jet forces that satisfy constraints

    config = q(3: end);  % Get current world angle + shape

    % Get optimal jet inputs
    B_aug_q = interpAugmentedDistribution(B_aug, config, domain);
    [c_i, A_ij] = getOptimalJetInputs(B_aug_q, q_dot_des, base_field_num, bias_constant_jets);
    
    % Solve for parameters
    lb = zeros([2 * base_field_num - 1, 1]);
    ub = cat(1, c_i, repmat(2*pi, base_field_num - 1, 1));
    params = lsqnonlin( ...
        @(x) parameterErrorFunction(x, omega, A_ij, base_field_num), ...
        zeros([5, 1])); ...
        %lb, ub);

    jet_forces = getJetForces(t, R, omega, c_i, params, base_field_num);
    local_forces = jet_forces;

end


function B_aug_q = interpAugmentedDistribution(B_aug, config, domain)
    %% Interpolate augmented distribution for a given configuration

    config = num2cell(config);
    B_aug_q = zeros([size(B_aug{1}, 1), numel(B_aug)]);
    for row_idx = 1:size(B_aug{1}, 1)

        % Specify slicing for the nd-array
        S.type = '()';
        S.subs = cat(2, {row_idx}, repmat({':'}, 1, ndims(B_aug{1}) - 1));

        % Iterate over vector fields and interpolate
        for field_idx = 1:numel(B_aug)
            B_aug_q_row_field = subsref(B_aug{field_idx}, S);
            size(B_aug_q_row_field)
            B_aug_q(row_idx, field_idx) = interpn(domain{:}, squeeze(B_aug_q_row_field), config{:}, 'cubic');
        end

    end
end


function error = parameterErrorFunction(x, omega, A_ij, base_field_num)
    %% Temporary measure to evaluate error function 

    % Assign variables
    k_1 = x(1);
    k_2 = x(2);
    k_3 = x(3);
    phi_2 = x(4);
    phi_3 = x(5);

    % Area terms it must satisfy
    area_cons = [
        pi * k_1 * k_2 * sin(phi_2);
        pi * k_1 * k_3 * sin(phi_3);
        pi * k_2 * k_3 * sin(phi_3 - phi_2)
        ];

    error = A_ij - omega * area_cons;  % Error term to minimize
end


function [c_i, A_ij] = getOptimalJetInputs(B_aug_q, q_dot_des, base_field_num, bias_constant_jets)
    %% Optimal jet inputs constants for the given desired generalized velocity

    if bias_constant_jets
        c_i = lsqminnorm(B_aug_q(:, 1:base_field_num), q_dot_des);
        q_dot_res = q_dot_des - B_aug_q(:, 1:base_field_num) * c_i;
        A_ij = lsqminnorm(B_aug_q(:, base_field_num + 1:end), q_dot_res);
    else
        size(B_aug_q)
        size(q_dot_des)
        u = lsqminnorm(B_aug_q, q_dot_des);
        c_i = u(1:base_field_num);
        A_ij = u(base_field_num + 1: end);
    end

end


% function [c_i, A_ij] = getOptimalJetInputs(B_aug_q, q_dot_des, base_field_num, bias_constant_jets)
%     %% Optimal jet inputs constants for the given desired generalized velocity
% 
%     if bias_constant_jets
%         lie_field_num = size(B_aug_q, 2) - base_field_num;
%         c_i = lsqnonlin(@(x) B_aug_q(:, 1:base_field_num)*x - q_dot_des, ...
%             zeros([base_field_num 1]), ...
%             zeros([base_field_num 1]), []);
%         q_dot_res = q_dot_des - B_aug_q(:, 1:base_field_num) * c_i;
%         A_ij = lsqminnorm(B_aug_q(:, base_field_num + 1:end), q_dot_res);
% %         A_ij = lsqnonlin(@(x) B_aug_q(:, base_field_num + 1:end)*x - q_dot_res, ...
% %             zeros([lie_field_num 1]), ...
% %             zeros([lie_field_num 1]), []);
%     else
%         size(B_aug_q)
%         size(q_dot_des)
%         u = lsqminnorm(B_aug_q, q_dot_des);
%         c_i = u(1:base_field_num);
%         A_ij = u(base_field_num + 1: end);
%     end
% 
% end


function jet_forces = getJetForces(t, R, omega, c_i, params, base_link_num)
    %% Get jet/local forces from control parameters

    k_i = params(1: base_link_num);
    phi_i = cat(1, 0, params(base_link_num + 1: end));

    u_norm = cell2mat(arrayfun(@(c, k, phi) k*cos(omega*t + phi) + c, c_i, k_i, phi_i, 'UniformOutput', false));
    jet_forces = R * u_norm;
end