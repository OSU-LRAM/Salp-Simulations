function [local_force_func, input_func] = getMinGaitForceModFunction(num_of_links, q_init, B_aug, domain, g_dot_des, thrust_angle, omega)
    %% Output force function for the jets 

    addpath('SalpUtils\');

    % Get gait related info
    r_init = q_init(4:end);
    shape_pos_func = GaitUtils.generateBasicGaitPos(omega, r_init);
    shape_vel_func = GaitUtils.generateBasicGaitVel(omega, r_init);
    gait_info = GaitUtils.evalGaitOverPeriod(shape_pos_func, shape_vel_func, (1/omega)*2*pi, 21);

    % Get inputs and forces
    input_func = getOptimalInputs(omega, g_dot_des, gait_info, B_aug, domain, num_of_links);
    local_force_func = InputOptimUtils.getJetOnlyForceFuncion(input_func, thrust_angle);
end


function input_func = getOptimalInputs(omega, g_dot_des, gait_info, B_aug, domain, num_of_links)
    %% Optimal jet inputs constants for the given desired generalized velocity


    % Set parameters for optimizer
    % TODO: Modify this later to use lsqnonlin with multistart instead of 
    % fmincon with globalsearch
    rng default
    gs = GlobalSearch;
    x0 = vertcat(zeros([num_of_links, 1]) + .5, randsample(100, ...
        num_of_links - 1)*(2*pi/100), zeros([num_of_links, 1]));
    error_min = @(x) parameterErrorFunction(x, omega, g_dot_des, gait_info, B_aug, domain, num_of_links);
    problem = createOptimProblem('fmincon', 'x0', x0, 'objective', error_min);
    params = run(gs, problem);

    % Get input function from paramters
    [k_vec, phi_vec, c_vec] = InputOptimUtils.parseParameterInputs(params, num_of_links);
    input_func = InputOptimUtils.getInputFunction(omega, k_vec, phi_vec, c_vec);
end


function error = parameterErrorFunction(x, omega, g_dot_des, gait_info, B_aug, domain, num_of_links)
    %% Error function to minimize for a given time

    [k_vec, phi_vec, c_vec] = InputOptimUtils.parseParameterInputs(x, num_of_links);
    area_vect = InputOptimUtils.generateAreaVector(k_vec, phi_vec);
    u = InputOptimUtils.getInputFunction(omega, k_vec, phi_vec, c_vec);

    % Get function to minimize
    % TODO: Find way to vectorize this later (vecnorm and passing the 
    % correct arguments to my functions
    error = 0;
    for idx = 1:numel(gait_info.t)
        
        % Get current evaluation paramters
        t = gait_info.t(idx);
        r = gait_info.shape_pos(:, idx);
        r_dot = gait_info.shape_vel(:, idx);

        % Get relevant parameters
        u_t = u(t);
        q_dot = vertcat(g_dot_des, r_dot);
        [B, B_Lie] = InputOptimUtils.interpAugB(B_aug, r, domain, num_of_links);

        % Get accumulated error
        error = error + norm(q_dot - (B*u_t + B_Lie*area_vect));
    end
end



