function varargout = animateDistribution_discrete(salp, solver_config, anim_config, verify)
% ANIMATESALP Summary of this function goes here
%   Detailed explanation goes here

    if (nargin == 4 && verify) || nargin == 3
        varargout{1} = verifyConfigurations(salp, solver_config, anim_config);
    end
    
    simul_results = integrateSalpMotion(salp, solver_config);
    plotSalpStates(salp, simul_results);
    animateSalp(salp, simul_results, anim_config);
end


function plotSalpStates(salp, simul_results)

    % Get revelant variables
    salp.baseframe = 'com-mean';
    ts = simul_results.ts;
    states = simul_results.states;
    alphas = sym('alpha_', [numel(salp.linklengths) - 1, 1]);
    [~, ~, ~, base_to_frame_trans] = N_link_chain(salp, alphas);

    % Get com states
    com_states = zeros(3, size(states, 2));
    for idx = 1:size(states, 2)
        com_states(:, idx) = subs(base_to_frame_trans, alphas, states(4:5, idx))* states(1:3, idx);
    end
    
    figure(1337)
    for idx = 1:5
        subplot(5, 1, idx)
        if idx < 4
            plot(ts, com_states(idx, :))
        else
            plot(ts, states(idx, :))
        end
    end
end


function verified_configs = verifyConfigurations(salp, solver_config, anim_config)

    %% Verify salp fields

    fields = ["drag_ratio", "baseframe", "linklengths", "num_of_links", ...
        "num_of_joints"];
    verifyRequiredFields(salp, "salp", fields);

    % Check optional fields
 
    if isfield(salp, "thruster_angles")
        disp("System will use the angles specified in 'thruster_angles' to" + ...
            "adjust the thruster angles");
    else
        disp("system will determine thruster angles per timestep");
    end

    if isfield(salp, "max_thrust")
        disp("System will use 'max_thrust' to normalize the jet forces");
    else
        disp("system will determine max thrust from the simulation and use" + ...
            "that to normalize the jet forces");
    end

    %% Verify solver configuration fields

    fields = ["t_span", "dt", "init_state", "local_force_func", ...
        "base_force_func"];
    verifyRequiredFields(solver_config, "solver_config", fields);

    %% Verify animation configuration fields

    fields = ["dir_loc", "gif_name", "tick_colors", "link_colors", ...
        "thruster_colors", "in_jet_colors", "out_jet_colors"];
    verifyRequiredFields(anim_config, "anim_config", fields);

    verified_configs = true;
end


function verifyRequiredFields(config, config_name, fields)
% Notify user that there are missing required fields

    for idx = 1:numel(fields)
        if ~isfield(config, fields)
            error(fields(idx) + " is a required field for " + config_name);
        end
    end
end