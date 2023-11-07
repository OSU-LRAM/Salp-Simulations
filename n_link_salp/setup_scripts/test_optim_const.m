
clc, close all

addpath('SalpUtils')
addpath("..\salp_modelling\")  % Add path of the salp modelling and animation
addpath("..\n_link_salp\SalpUtils\")

baseframe = 'com-mean';
load(['salp_3_link_basis_15_t_angle_30baseframe_' baseframe '.mat'])

%% Extract distribution

B_aug = cell(1, numel(poss_basis));
for idx = 1:numel(poss_basis)
    B_aug{idx} = real(poss_basis{idx}.field);
end

%% System properties for the salp

n = 3;  % Number of links
drag_ratio = 3;
linklengths = ones([n, 1]);
f_jets = cell(n, 1);
f_jets(:) = {zeros([1,3])};
f_jets{3} = [0 10 0];

salp = setupSalp(linklengths, drag_ratio, f_jets);

%% Jet force parameters

deg_of_freedom = 2;
T_vec_idxs = [2, 4];
T_vec = [.001 .01 .1 1 10 100 1000];
thrust_angle = pi/6;
max_force_thrust = 50;
bias_constant_jets = true;
q_dot_des = zeros([(n - 1) + 3, 1]);
q_dot_des(deg_of_freedom) = .4;

%% Setup solver configuration
% rng default
g_init = zeros(1, 3);  % Initial position in world frame
r_init = 2*pi*rand([1, n - 1]) - pi;  % Initial shape of the system
r_init = [-1, 1]*(pi/8);  % Initial shape of the system
q_init = [g_init r_init]';
fixed_f_ext = zeros([3 + (n - 1), 1]);  % Precomputed constant additional body 

%% Solver settings to plot resulting states
solver_config = struct();
solver_config.dt = 1/30;  % timestep for simulation (30 fps)
solver_config.t_span = [0, 2];
solver_config.init_state = [g_init r_init];
solver_config.base_force_func = @(t, X) fixed_f_ext;

%% Create plots
f_omegas = {@(omega) omega/(2*pi), @(omega) omega, @(omega) 1/omega, @(omega) 1, @(omega) 1/2, @(omega) 2, @(omega) 1/(omega * 16), @(omega) omega^2, @(omega) sqrt(omega), @(omega) exp(omega), @(omega) log(omega)};
% f_omega_inv_1 = [1, 1, nan, nan, 1/16, exp(1)];
f_omega_strs = {'\frac{\omega}{2\pi} = T^{-1}','\omega', '\frac{1}{\omega}', '1', '\frac{1}{2}', '2', '\frac{1}{16\omega}', '\omega^2', '\sqrt{\omega}', 'e^{\omega}', 'ln(\omega)'};
for idx = 1:numel(f_omegas)
    plotSinusoidParameters(T_vec, T_vec_idxs, salp, solver_config, f_omegas{idx}, f_omega_strs{idx}, n, thrust_angle, q_init, B_aug, domain, q_dot_des, max_force_thrust)
end

% Test functions

function plotSinusoidParameters(T_vec, T_vec_idxs, salp, solver_config, f_omega, f_omega_str, num_of_links, thrust_angle, q_init, B_aug, domain, q_dot_des, max_thrust)

    %% Get constants
    [~, A_const_vec] = extractOptimalJetContributions(num_of_links, q_init, B_aug, domain, q_dot_des);

%     % Check constants
%     r_init = q_init(4:end);
%     [B, B_Lie] = InputOptimUtils.interpAugB(B_aug,r_init,domain,num_of_links);
%      norm(q_dot_des - (B*c_const_vec2 + B_Lie*A_const_vec2)) > norm(q_dot_des - (B*c_const_vec1 + B_Lie*A_const_vec1))

    %% Generate parameters
    omega_vec = 2*pi*(T_vec).^-1;
    f_omega_vec = zeros(size(omega_vec));
    A_mat = zeros(numel(B_aug) - num_of_links, numel(T_vec));
    force_funcs = cell(numel(T_vec), 1);
    area_mat = zeros(size(A_mat));
    for idx = 1:numel(T_vec)

        % Get parameters
        omega = omega_vec(idx);
        [c_vec, A_vec] = extractOptimalJetContributions(num_of_links, q_init, B_aug, domain, q_dot_des);
        [k_vec, phi_vec] = getOptimalForceParameters(num_of_links, c_vec, A_vec/f_omega(omega), max_thrust);

        % Save parameter data
        f_omega_vec(idx) = f_omega(omega);
        A_mat(:, idx) = A_vec/f_omega(omega);
        area_mat(:, idx) = InputOptimUtils.generateAreaVector(k_vec, phi_vec);
        input_func = InputOptimUtils.createInputFunction(omega, k_vec, phi_vec, c_vec);
        force_funcs{idx} = InputOptimUtils.getJetOnlyForceFunction(input_func, thrust_angle, num_of_links);
    end


    %% Get resulting states for selected T's
    simul_results = cell(numel(T_vec_idxs), 1);
    for idx = 1:numel(T_vec_idxs)
        omega_idx = T_vec_idxs(idx);
        solver_config.local_force_func = force_funcs{omega_idx};
        simul_results{idx} = integrateSalpMotion(salp, solver_config);
    end


    %% Create plots

    % Plot constants
    plotConstantValues(T_vec, f_omega_vec, f_omega_str, area_mat, A_const_vec)

    % Plot resulting displacement
    for idx = 1:numel(T_vec_idxs)
        T = T_vec(T_vec_idxs(idx));
        plotSalpStates(T, salp, simul_results{idx}, f_omega_str)
    end
end


function [c_vec, A_vec] = searchGloballyMinimumJetContributions(num_of_links, q_init, B_aug, domain, q_dot_des)

    r_init = q_init(4:end);
    [B, B_Lie] = InputOptimUtils.interpAugB(B_aug,r_init,domain,num_of_links);

    % Build and solve optimization problem
    gs = GlobalSearch;
    x0 = zeros(num_of_links + size(B_Lie, 2), 1);
    lsq_norm = @(x) norm(q_dot_des - [B B_Lie]*x);
    problem = createOptimProblem('fmincon', ...
        'x0', x0, 'objective', lsq_norm);

    % Solve and separate parameters
    params = run(gs, problem);
    c_vec = params(1:num_of_links);
    A_vec = params(num_of_links + 1:end);
end


function plotConstantValues(T_vec, f_omega_vec, f_omega_str, area_mat, A_const_vec)
    %% Create plots highlighting constant values that preserved by the optimizer

    % Get values for plotting
    f_omega_pts = area_mat .* f_omega_vec;
    ylabel_sscript = {'{12}', '{13}', '{23}'};

    % Get subplots highlighting constant values being preserved
    figure
    sgtitle(['$f(\omega) = ' f_omega_str '$'], ...
        "Interpreter", 'latex')
    for idx = 1:size(area_mat, 1)
        subplot(size(area_mat, 1), 1, idx)
        scatter(log10(T_vec), f_omega_pts(idx, :), ...
            "MarkerFaceColor", 'w', ...
            "MarkerEdgeColor", 'k')
        yline(A_const_vec(idx), ...
            "Color", 'r', ...
            "Label", ['$C(r_{init}, \dot{q}_{des}) = ' num2str(A_const_vec(idx)) '$'], ...
            "Interpreter", "latex", ...
            "LineStyle", "--", ...
            "LabelHorizontalAlignment","left");
        xlim([-3 3])
        ylim([A_const_vec(idx) - 5, A_const_vec(idx) + 5])
        xlabel('$log(T)$', "Interpreter", 'latex', "FontWeight", 'bold');
        ylabel(['$f(\omega)A_' ylabel_sscript{idx} '$'], ...
            "FontWeight", 'bold', "Interpreter",'latex')
    end
end


function plotSalpStates(T, salp, simul_results, f_omega_str)

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
    
    % Create state plot
    figure
    ylabels = {'x', 'y', '\theta', '\alpha_1', '\alpha_2'};
    sgtitle(['$f(\omega) = ' f_omega_str ', T = ' num2str(T) '$'], ...
        "Interpreter", 'latex');
    for idx = 1:5
        subplot(5, 1, idx)
        if idx < 4
            plot(ts, com_states(idx, :))
        else
            plot(ts, states(idx, :))
        end
        xlabel('$t$', "Interpreter", 'latex', "FontWeight", 'bold');
        ylabel(['$' ylabels{idx} '$'], "FontWeight", 'bold', "Interpreter",'latex')
    end
end


function [c_vec, A_vec] = extractOptimalJetContributions(num_of_links, q_init, B_aug, domain, q_dot_des)
    %% Output force function for the jets

    r_init = q_init(4:end);
    [B, B_Lie] = InputOptimUtils.interpAugB(B_aug,r_init,domain,num_of_links);
    [c_vec, A_vec] = getOptimalJetContributions(B, B_Lie, q_dot_des);
end


function [c_vec, A_vec] = getOptimalJetContributions(B, B_Lie, q_dot_des)
    %% Get a linear combination of c's and A's that results in the desired
    % average velocity 

    u = lsqminnorm([B, B_Lie], q_dot_des);
    c_vec = u(1:size(B, 2));
    A_vec = u(size(B, 2) + 1: end); 
end


function [k_vec, phi_vec] = getOptimalForceParameters(num_of_links, c_vec, area_vec, max_thrust)
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
    ub = vertcat(ub, zeros(numel(area_vec), 1) + pi);
    lb = vertcat(lb, zeros(numel(area_vec), 1) - pi);

    % Create and run optimization problem
    ms = MultiStart('UseParallel', true);
    error_min = @(x) parameterErrorFunction(x, area_vec, num_of_links);
    problem = createOptimProblem('lsqnonlin', 'x0', x0, 'objective', error_min, ...
        'ub', ub, 'lb', lb);
    parpool
    params = run(ms, problem, 100);
    delete(gcp('nocreate'))

    % Get input function from paramters
    [k_vec, phi_vec] = InputOptimUtils.parseParameterInputs(params, num_of_links);
end


function error = parameterErrorFunction(x, area_vec, num_of_links)
    %% Error function to minimize for a given time

    [k_vec, phi_vec] = InputOptimUtils.parseParameterInputs(x, num_of_links);
    area_vect = InputOptimUtils.generateAreaVector(k_vec, phi_vec);
    
    error = zeros(size(area_vec));
    for idx = 1:numel(area_vec)
        if area_vec > 0 
            error(idx) = area_vec(idx) + area_vect(idx);
        else
            error(idx) = area_vec(idx) - area_vect(idx);
        end
    end
end