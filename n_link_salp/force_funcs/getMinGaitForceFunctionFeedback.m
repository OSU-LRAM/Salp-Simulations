function period_force_func = getMinGaitForceFunctionFeedback(T, num_of_links, q_init, B_aug, domain, q_dot_des, thrust_angle, omega, max_thrust, bias_jets)
    %% Wrapper to generate force functions at every time period
    
    % Local global workspace variables to generate new force functions
    % everytime a period passes
    periods_passed = 0;  % Numbers of periods passed
    next_period_time = T;  % Variable to keep track if a period has passed
    force_func = getMinGaitForceFunction2(num_of_links, q_init, B_aug, domain, q_dot_des, thrust_angle, omega, max_thrust, bias_jets);

    function [local_forces, jet_forces] = periodicForceSpeedFeedback(t, T, q, num_of_links, q_init, B_aug, domain, q_dot_des, thrust_angle, omega, max_thrust, bias_jets)
        %% Actual force function executor

        if next_period_time < t
            next_period_time = (periods_passed + 2) * T;
            periods_passed = periods_passed + 1;
            force_func = getMinGaitForceFunction(num_of_links, q_init, B_aug, domain, q_dot_des, thrust_angle, omega, max_thrust, bias_jets);
        end

        [local_forces, jet_forces] = force_func(t, q);
    end

    period_force_func = @(t, q) periodicForceSpeedFeedback(t, T, q, num_of_links, q_init, B_aug, domain, q_dot_des, thrust_angle, omega, max_thrust, bias_jets);
end

