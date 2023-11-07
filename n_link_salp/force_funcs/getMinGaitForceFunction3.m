function local_force_func = getMinGaitForceFunction3(num_of_links, q_init, B_aug, domain, q_dot_des, thrust_angle, omega, max_thrust)
    %% Output force function for the jets

    addpath('SalpUtils')

    % Run input force procedure
    u = getOptimalJetContributions(omega, B_aug, domain, q_init, q_dot_des, num_of_links, max_thrust);
    force_funcs = buildJetForceFunctions(omega, u, num_of_links, thrust_angle);
    local_force_func = createConsecutiveForces(T, force_funcs);
end


function cyclic_force_func = createConsecutiveForces(T, force_funcs)

    periods_passed = 0;  % Numbers of periods passed
    next_period_time = T;  % Variable to keep track if a period has passed
    force_func = force_funcs{1};

    function [local_forces, jet_forces] = cyclicForceFunction(t, q)
        %% Actual force function executor

        if next_period_time < t
            next_period_time = (periods_passed + 2) * T;
            periods_passed = periods_passed + 1;
            force_func = force_funcs{mod(periods_passed, numel(force_funcs)) + 1};
        end

        [local_forces, jet_forces] = force_func(t, q);
    end

    cyclic_force_func = @cyclicForceFunction;
end


function u = getOptimalJetContributions(omega, B_aug, domain, q_init, q_dot_des, num_of_links, max_thrust)
    %% Get a linear combination of c's and A's that results in the desired
    % average velocity 

    % Initialize parameters
    counter = 0
    u = {};
    q = q_init;
    T = (2*pi)/omega;
    q_dot_res = q_dot_des;

    while norm(q_dot_res) > 0.001

        % Parameters for current iteration
        u_i = 0;
        q_i = 0;
        r = q(4:end);
        q_dot_res_i = q_dot_res;
        [B, B_Lie] = InputOptimUtils.interpAugB(B_aug,r,domain,num_of_links);

        % Get most significant input displacer
        for idx = 1:numel(B_aug)

            % Get relevant values for current subiteration
            if idx <= num_of_links
                u_ik = lsqminnorm(B(:, idx), q_dot_res) ;
%                 u_ik = lsqlin(B(:, idx), q_dot_res, [], [], [], [], 0, max_thrust);
                q_dot_res_ik = q_dot_res - B(:, idx)*u_ik;
                q_ik = T*q_dot_res_ik + q;
            else
                col_idx = idx - num_of_links;
                B_Lie_k = B_Lie(:, col_idx)/norm(B_Lie(:, col_idx));
                u_ik = lsqminnorm(B_Lie_k, q_dot_res);
%                 u_ik = lsqlin(B_Lie(:, col_idx), q_dot_res, [], [], [], [], 0, pi*max_thrust^2);
                q_ik = q + u_ik*B_Lie_k;
                q_dot_res_ik = q_dot_res - u_ik*B_Lie_k/T;
            end

            % Update parameters if these have a better displacment
            if norm(q_dot_res_ik) < norm(q_dot_res_i)
                q_i = q_ik;
                u_i = {idx, u_ik};
                q_dot_res_i = q_dot_res_ik;
            end
        end
        
        % Update general parameters
        q = q_i;
        u{end + 1} = u_i;
        q_dot_res = q_dot_res_i; 
    end
end


function force_funcs = buildJetForceFunctions(omega, u, num_of_links, thrust_angle)

    force_funcs = cell(numel(u), 1);

    for idx = 1:numel(u)
        u_idx = u{1};
        u_i_const = u{2};
        if u_idx <= num_of_links
            k_vec = zeros(num_of_links, 1);
            phi_vec = zeros(num_of_links, 1);
            c_vec = zeros(num_of_links, 1);
            c_vec(u_idx) = u_i_const;
            u_i = InputOptimUtils.createInputFunction(0, k_vec, phi_vec, c_vec);
        else
            k_vec = zeros(num_of_links, 1);
            phi_vec = zeros(num_of_links, 1);
            c_vec = zeros(num_of_links, 1);

            % This is temporary, but this get correct indexes to overwrite
            if u_idx == 4
                mid_idxs = [1 3];
            else
                mid_idxs = [2 3];
            end
            [k_vec_i, phi_vec_i] = getOptimalSinusoidParameters(omega, u_i_const, max_thrust);

            for idx2 = 1:2
                mid_idx = mid_idxs(idx2);
                k_vec(mid_idx) = k_vec_i(idx2);
                phi_vec(mid_idx) = phi_vec_i(idx2);
            end
            u_i = InputOptimUtils.createInputFunction(omega, k_vec, phi_vec, c_vec);
        end
        force_funcs{idx} = InputOptimUtils.getJetOnlyForceFunction(u_i, thrust_angle, num_of_links);
    end
end


function [k_vec, phi_vec] = getOptimalSinusoidParameters(omega, A_vec, max_thrust)
    %% Optimal jet input functions for the given the jet contributions

    % Set bounds for amplitude oscillations (these conditions lock the salp
    % jets)
    ub = zeros(2, 1);
    lb = zeros(2, 1);
%     for idx = 1:num_of_links
%         if c_vec(idx) > 0
%             ub(idx) = min(max_thrust, c_vec(idx));
%         elseif c_vec(idx) < 0
%             lb(idx) = max(-max_thrust, c_vec(idx));
%         end
%     end

    % Set parameters for optimizer
    x0 = zeros([3, 1]);
    ub = vertcat(ub, zeros(1, 1) + pi);
    lb = vertcat(lb, zeros(1, 1) - pi);

    % Create and run optimization problem
    ms = MultiStart('UseParallel', true);
    error_min = @(x) parameterErrorFunction(x, omega, A_vec, num_of_links);
    problem = createOptimProblem('lsqnonlin', 'x0', x0, 'objective', error_min, ...
        'ub', ub, 'lb', lb); %, ...
%         'nonlcon', @(x) evalNonlinearConditions(x, num_of_links));
    parpool
    params = run(ms, problem, 100);
    delete(gcp('nocreate'))

    % Get input function from paramters
    [k_vec, phi_vec] = InputOptimUtils.parseParameterInputs(params, num_of_links);
end


function [c, ceq] = evalNonlinearConditions(x, num_of_links)
    
    [k_vec, ~, c_vec] = InputOptimUtils.parseParameterInputs(x, num_of_links);

    c = k_vec - c_vec;
    ceq = [];
end


function error = parameterErrorFunction(x, A_vec)
    %% Error function to minimize for a given time

    [k_vec, phi_vec] = InputOptimUtils.parseParameterInputs(x, 2);
    area_vect = InputOptimUtils.generateAreaVector(k_vec, phi_vec);
    
    error = zeros(size(A_vec));
    for idx = 1:numel(A_vec)
        if A_vec > 0 
            error(idx) = A_vec(idx) + area_vect(idx);
        else
            error(idx) = A_vec(idx) - area_vect(idx);
        end
    end
end
