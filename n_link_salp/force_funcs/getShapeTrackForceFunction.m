function local_force_func = getShapeTrackForceFunction(num_of_joints, q_init, q_dot_des, omega, B_aug, domain, thrust_angle)
    %% Generate jet input force function

    [~, B_Lie] = interpAugmentedDistribution(B_aug, q_init, domain, num_of_joints);
    gait_vel_func = defineGaitVelocities(omega, B_Lie, q_dot_des, num_of_joints);
    local_force_func = @(t, q) getOptimalInputs(num_of_joints, t, q, B_aug, domain, gait_vel_func, thrust_angle);
end


function [local_forces, jet_forces] = getOptimalInputs(num_of_joints, t, q, B_aug, domain, gait_vel_func, thrust_angle)
    %% Optimal jet inputs constants for the given desired generalized velocity

    shape = q(4:end);
    q_dot_dim = size(B_aug{1}, 1);

    % Get magnitude of forces
    B_q = interpAugmentedDistribution(B_aug, shape, domain, num_of_joints);
    q_field_vel = cat(1, zeros([q_dot_dim - num_of_joints, 1]), gait_vel_func(t));
    force_norm = lsqminnorm(B_q, q_field_vel);

    % Get big rotation matrix that maps thrust input thrust vector to forces
    R = zeros([3*(num_of_joints + 1), num_of_joints + 1]);
    for idx = 1:num_of_joints + 1
        if mod(idx, 2) == 0
            angle = thrust_angle;
        else
            angle = -thrust_angle;
        end
        R(3*(idx-1)+1: 3*(idx-1)+3, idx) = [cos(angle); sin(angle); 0];
    end

    % Get corresponding forces
    local_forces = R*force_norm;
    jet_forces = local_forces;
end


function [B_base, B_Lie] = interpAugmentedDistribution(B_aug, config, domain, num_of_joints)
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

    % Separate both B_base (the base fields) and B_Lie (Lie bracket fields)
    B_base = B_aug_q(:, 1:num_of_joints + 1);
    B_Lie = B_aug_q(:, num_of_joints + 2: end);
end


function gait_vel_func = defineGaitVelocities(omega, B_Lie, q_dot_des, num_of_joints)
    %% Return shape velocity function that form a gait

    % Set parameters for optimizer
    rng default
    gs = GlobalSearch;
    x0 = cat(1, zeros([num_of_joints, 1]) + .5, randsample(100, ...
        num_of_joints - 1)*(2*pi/100));
    error_min = @(x) parameterErrorFunction(x, omega, B_Lie, q_dot_des, num_of_joints);
    problem = createOptimProblem('fmincon', 'x0', x0, 'objective', error_min, ...
        'ub', [15, 15, pi], 'lb', [-15, -15, -pi], 'nonlcon', @equalizeAmplitudes);
    params = run(gs, problem);

    % Get parameters from optimizer solution
    k_vec = params(1:num_of_joints);
    phi_vec = cat(1, 0, params(num_of_joints + 1: end));

    % Define gait velocities with the given parameters
    shape_vel = cell(num_of_joints, 1);
    for idx = 1:numel(shape_vel)
        shape_vel{idx} = @(t) -omega*k_vec(idx)*cos(omega*t + phi_vec(idx) - pi/2);
    end

    gait_vel_func = @(t) evalShapeVelocity(t, shape_vel);
end


function gait_vel = evalShapeVelocity(t, shape_vel)
    %% Evaluate the shape velocity for a given time t

    gait_vel = zeros(size(shape_vel));
    for idx = 1:numel(shape_vel)
        gait_vel(idx) = shape_vel{idx}(t);
    end
end


function [c, ceq] = equalizeAmplitudes(x)
    %% Constraint to equalize the resulting amplitudes in optimization

    c = [];
    ceq = norm(x(1) - x(2));
end


function error = parameterErrorFunction(x, omega, B_Lie, q_dot_des, num_of_joints)
    %% Find the parameters that form a gait and satisfy velocity requirements

    q_dot_dim = size(B_Lie, 1);

    % Assign variables
    k_vec = x(1:num_of_joints);
    phi_vec = cat(1, 0, x(num_of_joints + 1: 2*num_of_joints - 1));

    % Get relevant vectors to minimize along q_dot_des
    shape_vel_vect = generateShapeVelParameterVector(k_vec, phi_vec, ...
        num_of_joints, q_dot_dim);
    area_comp_vect = generateAreaVector(k_vec, phi_vec, B_Lie);

    error = norm(q_dot_des - omega * (shape_vel_vect + area_comp_vect));
end


function shape_vel_vect = generateShapeVelParameterVector(k_vec, phi_vec, num_of_joints, q_dot_dim)
    %% Generate desired shape velocities based on the shape space parameterization

    shape_vel_vect = zeros([num_of_joints, 1]);
    for idx = 1:num_of_joints
        shape_vel_vect(idx) = -k_vec(idx) * cos(phi_vec(idx) - pi/2);
    end
    shape_vel_vect = cat(1, zeros(q_dot_dim - num_of_joints, 1), shape_vel_vect);
end


function area_comp_vect = generateAreaVector(k_vec, phi_vec, B_Lie)
    %% Get the resulting velocities from the Lie bracket fields

    % TODO: I need to adjust this method (i.e., idx2 = idx + 1:numel(k_vec)

    area_vect = [];
    for idx = 1:numel(k_vec)
        for idx2 = idx:numel(k_vec)
            area_vect(end + 1) = pi * k_vec(idx) * k_vec(idx2) * sin(phi_vec(idx2) - ...
                phi_vec(idx));
        end
    end
    area_comp_vect = B_Lie * area_vect';
end
