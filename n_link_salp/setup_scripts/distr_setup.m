
addpath("..\salp_modelling\")  % Add path of the salp modelling and animation


%% System properties for the salp

drag_ratio = 3;
linklengths = ones([3, 1]);
f_jets = {[0 0 0], [0 0 0], [1 0 0]};

salp = setupSalp(linklengths, drag_ratio, f_jets);

%% Setup graphing configuration

graph_config = struct();
graph_config.delta = .5;  % Spacing for graph
graph_config.range = [-1, 1];  % Range of axes
graph_config.skip_zeros = true;  % Skip graphing f jet components that are zero
graph_config.r = true;  % Specifies plotting shape distribution
graph_config.g = true;  % Specifies plotting displacement distribution

%% Graph distribution

% plotDistribution_discrete(salp, graph_config);

%% Setup solver configuration

g_init = zeros(1, 3);  % Initial position in world frame
r_init = [0.5, 0.5];  % Initial shape of the system
fixed_f_ext = zeros([5, 1]);  % Precomputed constant additional
fixed_f_jets = [salp.f_jets{:}]';  % Precomputed constant jet forces

solver_config = struct();
solver_config.dt = 1/30;  % timestep for simulation (30 fps)
solver_config.t_span = [0, 5];
solver_config.init_state = [g_init r_init];
% solver_config.jet_force_func = @(t, X) fixed_f_jets;  % This has to be consistent to how jet forces are generated in local_force_func 
solver_config.base_force_func = @(t, X) fixed_f_ext;
% solver_config.local_force_func = @(t, X) fixed_f_jets;  % Forces that are applied locally on the links
solver_config.local_force_func = @(t, X) local_fixed_force(fixed_f_jets);  % Forces that are applied locally on the links

%% Animation configuration

anim_config = struct();

anim_config.dir_loc = 'C:\Users\Yousef\Documents\University\Research\Misc scripts\salp_animation\animations';
anim_config.gif_name = 'test_name'; 
anim_config.tick_colors = [.8,.8,.8];  % Grey
anim_config.link_colors = [0.0745,0.6235,1.0000];  % Light Blue
anim_config.thruster_colors = [0.4941,0.1843,0.5569];  % Purple
anim_config.in_jet_colors = [1.0000,1.0000,0.0667];  % Yellow
anim_config.out_jet_colors = [1.0000,0.4118,0.1608];  % Orange

animateDistribution_discrete(salp, solver_config, anim_config);