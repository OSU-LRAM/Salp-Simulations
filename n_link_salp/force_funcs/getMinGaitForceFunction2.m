function local_force_func = getMinGaitForceFunction2(num_of_links, q_init, B_aug, domain, q_dot_des, thrust_angle, omega, max_thrust, bias_jets)
    %% Output force function for the jets

    addpath('SalpUtils')

    % Default to not bias jet contributions
    if nargin == 8
        bias_jets = false;
    end

    % Run input force procedure
    r_init = q_init(4:end);
    [B, B_Lie] = InputOptimUtils.interpAugB(B_aug,r_init,domain,num_of_links);
    [c_vec, A_vec] = getOptimalJetContributions(omega, B, B_Lie, q_dot_des, bias_jets);
    [k_vec, phi_vec] = getOptimalForceParameters(omega, num_of_links, c_vec, A_vec, max_thrust);
    input_func = InputOptimUtils.createInputFunction(omega, k_vec, phi_vec, c_vec);
    local_force_func = InputOptimUtils.getJetOnlyForceFunction(input_func, thrust_angle, num_of_links);
end


function [c_vec, A_vec] = getOptimalJetContributions(omega, B, B_Lie, q_dot_des, bias_jets)
    %% Get a linear combination of c's and A's that results in the desired
    % average velocity 

    unit_B_Lie = zeros(size(B_Lie));
    for idx = 1:size(B_Lie, 2)
        unit_B_Lie(:, idx) = B_Lie(:, idx)/norm(B_Lie(:, idx));
    end

    if bias_jets
        c_vec = lsqminnorm(B, q_dot_des);
        q_dot_res = q_dot_des - B*c_vec;
        A_vec = lsqminnorm((1/omega)*B_Lie, q_dot_res);
    else
        u = lsqminnorm([B, unit_B_Lie], q_dot_des);
        c_vec = u(1:size(B, 2));
        A_vec = u(size(B, 2) + 1: end); 
    end
end


function [k_vec, phi_vec] = getOptimalForceParameters(omega, num_of_links, c_vec, A_vec, max_thrust)
    %% Optimal jet input functions for the given the jet contributions

    % Set bounds for amplitude oscillations (these conditions lock the salp
    % jets)
    ub = zeros(num_of_links, 1);
    lb = zeros(num_of_links, 1);
    for idx = 1:num_of_links
        if c_vec(idx) > 0
            ub(idx) = min(max_thrust, c_vec(idx));
        elseif c_vec(idx) < 0
            lb(idx) = max(-max_thrust, c_vec(idx));
        end
    end

    % Set parameters for optimizer
    x0 = zeros([2*num_of_links - 1, 1]);
    ub = vertcat(ub, zeros(numel(A_vec), 1) + pi);
    lb = vertcat(lb, zeros(numel(A_vec), 1) - pi);

    % Create and run optimization problem
    ms = MultiStart('UseParallel', true);
    error_min = @(x) parameterErrorFunction(x, omega, A_vec, num_of_links);
    problem = createOptimProblem('lsqnonlin', 'x0', x0, 'objective', error_min, ...
        'ub', ub, 'lb', lb);
    parpool
    params = run(ms, problem, 100);
    delete(gcp('nocreate'))

    % Get input function from paramters
    [k_vec, phi_vec] = InputOptimUtils.parseParameterInputs(params, num_of_links);
end


function error = parameterErrorFunction(x, omega, A_vec, num_of_links)
    %% Error function to minimize for a given time

    [k_vec, phi_vec] = InputOptimUtils.parseParameterInputs(x, num_of_links);
    area_vect = InputOptimUtils.generateAreaVector(k_vec, phi_vec);
    
    error = zeros(size(A_vec));
    for idx = 1:numel(A_vec)
        if A_vec > 0 
            error(idx) = A_vec(idx)/omega + area_vect(idx);
        else
            error(idx) = A_vec(idx)/omega - area_vect(idx);
        end
    end
end
