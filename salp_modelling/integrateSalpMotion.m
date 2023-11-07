function simul_results = integrateSalpMotion(salp, solver_config)
%% Solve for the salp motion
%
%   Inputs:
%
%       salp - Structure with information about the salp
%       solver_config - Information about the solver
%
%   Outputs:
%
%       simul_results - Structure with simulation results

    % Simulation results
    simul_results = struct();

    % Generate time vector for output
    simul_results.dt = solver_config.dt;
    simul_results.ts = 0:solver_config.dt:solver_config.t_span(end);
    
    % Set ODE description for swimming system
    % sets (forcingFunction,L,d,K) for ODE, changes to them after
    % this point will not be seen in ODE results
    odefun = @(t,X) salpODE(t, X, salp, solver_config);

    % Integrate swimmer ODE over desired time period given initial
    % conditions
    sol = ode45(odefun,solver_config.t_span, solver_config.init_state);

    % If we want the limit cycle results
    % Evaluate solution to ODE at desired time steps for ease of animation
    simul_results.states = deval(sol, simul_results.ts);

    % Get normalized jet force magnitudes and angles
    simul_results.forces = getNormalizedForceJets(salp, solver_config, simul_results);

end



function forces = getNormalizedForceJets(salp, solver_config, simul_results)
%% Get normalized jet forces in the system
%
%   Inputs:
%
%       ts - Time instances used for the solver
%
%   Outputs:
%
%       forces - Normalized jet forces for the system

    %     % Create storage for force information
    %     forces = struct();
    %     
    %     % Populate fields
    %     forces.jets = zeros(salp.num_of_links, numel(ts));
    %     forces.thruster_angles = cell(numel(ts), 1);

    % Variables to create and store normalized forces
    ts = simul_results.ts;
    forces = createForceJetStruct(salp, ts);

    % Get maximum norm in matrix and thruster angles
    max_norm = 0;
    for col_idx = 1:numel(ts)
 
        t = ts(col_idx);
        state = simul_results.states(:, col_idx);

        % Get force jets for current time
        [~, jets] = solver_config.local_force_func(t, state);
        
        % Create storage for thruster angles at ts
        thruster_angles = zeros(salp.num_of_links, 1);

        for row_idx = 1:salp.num_of_links
            comp_idx = 3*(row_idx - 1) + 1;

            % Update thruster angle and norm
            [thruster_angle, jet_norm] = getThrusterAngle(salp, jets(comp_idx: comp_idx + 1), row_idx);
            thruster_angles(row_idx) = thruster_angle;
            forces.jets(row_idx, col_idx) = jet_norm;

            % Keep track of max norm
            if jet_norm > max_norm
                max_norm = jet_norm;
            end
        end

        % Save thruster angles for current time ts
        forces.thruster_angles{col_idx} = thruster_angles;
    end

    % If max norm was specified, pass to system
    if isfield(forces, "max_thrust")
        max_norm = solver_config.max_thrust;
    end

    % Normalize jet forces from simulation
    disp(forces.jets);
    forces.jets = forces.jets/max_norm;
    disp(forces.jets);

end


function forces = createForceJetStruct(salp, ts)
%% Create storage for force information

    forces = struct();
    
    % Populate fields
    forces.thruster_angles = cell(numel(ts), 1);
    forces.jets = zeros(salp.num_of_links, numel(ts));
end


function [thruster_angle, jet_norm] = getThrusterAngle(salp, force_jet, link_num)
% Calculate thruster angles for a given force jet

    jet_norm = norm(force_jet);
    if isfield(salp, "thruster_angles") && salp.thruster_angles
        thruster_angle = salp.thruster_angles(link_num);
    else
        thruster_angle = atan2(force_jet(2), force_jet(1));
    end

end


function X_dot = salpODE(t, X, salp, solver_config)
%% Describes ODE for salp swimmer
%
%   Inputs:
%
%       t - Current simulation time
%       X - Current state of the system
%       salp - Structure with information about the salp
%       solver_config - Information about the solver
%
%   Outputs:
%
%       X_dot - State derivative vector

    %Calculate swimmer body/shape velocity from thrust at time t
    speeds = thrustToSpeeds(t, X, salp, solver_config);

    %Decompose velocity into body velocity and shape velocity
    g_circ = speeds(1:3);
    r_dot = speeds(4:end);

    %Turn swimmer body velocity into swimmer world velocity
    g_dot = TeLg(X(1:3))*g_circ;

    %Compose state rate of change
    %Position rate of change is swimmer world velocity
    %Shape rate of change is swimmer joint velocity
    X_dot = [g_dot;r_dot];

end


function speeds = thrustToSpeeds(t, X, salp, solver_config)
%Calculate body and joint velocity given a swimmer shape and set of thrust
%vectors acting on each link. Requires access to sysplotter functions
%
%   Outputs:
%
%       speeds - Local base link and shape velocities 

    r = X(4:end);  % Get current salp shape

    %Calculate body jacobian for each link
    [~,~,J_full] = N_link_chain(salp, r);

    %Calculate distribution matrix and drag metric
    J_Dual_full = getDualJacobians_discrete(J_full);
    M = getDragMetric(salp.linklengths, salp.drag_ratio, J_full);

    %Collect forces acting on local link frames and base frames
    local_forces = solver_config.local_force_func(t, X);
    forces = J_Dual_full*local_forces + solver_config.base_force_func(t, X);

    %Calculate body frame and joint velocities that would result from these
    %forces
    speeds = M\forces;

end

