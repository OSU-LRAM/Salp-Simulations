classdef GaitUtils
    %% Static utility class to help with salp gait generation/manipulation,
    % etc.

    methods(Static)

        function gait_info = evalGaitOverPeriod(shape_pos_func, shape_vel_func, T, num_of_pts)
            %% Evaluate the functions over a period T

            shape_dim = numel(shape_pos_func(0));

            gait_info.t = linspace(0, T, num_of_pts);  % Get time domain
            
            shape_pos = zeros([shape_dim, numel(gait_info.t)]);
            shape_vel = zeros([shape_dim, numel(gait_info.t)]);
            for idx = 1:shape_dim
                t = gait_info.t(idx);
                gait_info.shape_pos(:, idx) = shape_pos_func(t);
                gait_info.shape_vel(:, idx) = shape_vel_func(t);
            end
            gait_info.shape_pos = shape_pos;
            gait_info.shape_vel = shape_vel;
        end
        
        function shape_pos_func = generateBasicGaitPos(omega, r_init)
            %% Generate basic shape position sinusoidal functions for the salp to follow

            % Define individual shape functions
            shape_funcs = cell(numel(r_init), 1);
            for idx = 1:numel(r_init)
                shape_funcs{idx} = @(t) r_init(idx)/2 + (r_init(idx)/2) * cos(omega * t);
            end

            shape_pos_func = @(t) GaitUtils.evalParameterFunctions(t, shape_funcs);
        end


        function shape_vel_func = generateBasicGaitVel(omega, r_init)
            %% Generate basic shape velocity sinusoidal functions for the salp to follow

            % Define individual shape functions
            shape_funcs = cell(numel(r_init), 1);
            for idx = 1:numel(r_init)
                shape_funcs{idx} = @(t) -omega * (r_init(idx)/2) * cos(omega * t - pi/2);
            end

            shape_vel_func = @(t) GaitUtils.evalParameterFunctions(t, shape_funcs);
        end


        function values = evalParameterFunctions(t, para_funcs)
            %% Evaluate a set t-parameterized functions
                
            values = zeros([numel(para_funcs), 1]);
            for idx = 1:numel(para_funcs)
                values(idx) = para_funcs{idx}(t);
            end

        end

    end
end

