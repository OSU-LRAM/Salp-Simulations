classdef InputOptimUtils
    %% Static utility class with methods to perform input optimization on
    % the salp
    
    methods(Static)


        function local_force_func = getJetOnlyForceFunction(input_func, thrust_angle, num_of_links)

            force_func = InputOptimUtils.getJetForceFunction(input_func, thrust_angle, num_of_links);

            % Generate jet force in the salp animation code format
            function [local_forces, jet_forces] = generateJetForce(t, X)
                jet_forces = force_func(t);
                local_forces = jet_forces;
            end

            local_force_func = @generateJetForce;
        end 


        function force_func = getJetForceFunction(input_funcs, thrust_angle, num_of_links)

            % Get matrix that places inputs in direction of the forces
            R = zeros([3*num_of_links, num_of_links]);
            for idx = 1:num_of_links
                if mod(idx, 2) == 0
                    angle = thrust_angle;
                else
                    angle = -thrust_angle;
                end
                R(3*(idx-1)+1: 3*(idx-1)+3, idx) = [cos(angle); sin(angle); 0];
            end

            % Get force function from the inputs
            force_func = @(t) R * InputOptimUtils.evalParameterFunctions(t, input_funcs);
        end

        

        function [k_vec, phi_vec, c_vec] = parseParameterInputs(x, num_of_links)
            %% Helper function to parse parameters for input function

            k_vec = x(1:num_of_links);
            phi_vec = vertcat(0, x(num_of_links + 1: 2*num_of_links - 1));
            if nargout == 3
                c_vec = x(2*num_of_links: end);
            end
            
        end


        function input_func_t = getInputFunction(omega, k_vec, phi_vec, c_vec)
            %% Evaluate input function at the given time t

            input_func = InputOptimUtils.createInputFunction(omega, k_vec, phi_vec, c_vec);
            input_func_t = @(t) InputOptimUtils.evalParameterFunctions(t, input_func);
        end


        function input_funcs = createInputFunction(omega, k_vec, phi_vec, c_vec)
            %% Create function that generates inputs for a given time t

            input_funcs = cell(numel(k_vec), 1);
            for idx = 1:numel(k_vec)
                input_funcs{idx} = @(t) k_vec(idx) * cos(omega*t + phi_vec(idx)) + c_vec(idx);
            end
        end


        function [B, B_Lie] = interpAugB(B_aug, config, domain, num_of_links)
            %% Interpolate augmented distribution for a given configuration
        
            config = num2cell(config);
            B_aug_q = zeros([size(B_aug{1}, 1), numel(B_aug)]);
            for row_idx = 1:size(B_aug{1}, 1)
        
                % Specify slicing for the nd-array
                S.type = '()';
                S.subs = horzcat({row_idx}, repmat({':'}, 1, ndims(B_aug{1}) - 1));
        
                % Iterate over vector fields and interpolate
                for field_idx = 1:numel(B_aug)
                    B_aug_q_row_field = subsref(B_aug{field_idx}, S);
                    size(B_aug_q_row_field)
                    B_aug_q(row_idx, field_idx) = interpn(domain{:}, squeeze(B_aug_q_row_field), config{:}, 'spline');
                end
            end
        
            % Separate both B_base (the base fields) and B_Lie (Lie bracket fields)
            B = B_aug_q(:, 1:num_of_links);
            B_Lie = B_aug_q(:, num_of_links + 1: end);
        end
    
    
        function area_vect = generateAreaVector(k_vec, phi_vec)
            %% Get the resulting velocities from the Lie bracket fields
        
            area_vect = [];
            for idx = 1:numel(k_vec)
                for idx2 = idx + 1:numel(k_vec)
                    area_vect(end + 1) = pi * k_vec(idx) * k_vec(idx2) * sin(phi_vec(idx2) - ...
                        phi_vec(idx));
                end
            end
            area_vect = area_vect.';
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

