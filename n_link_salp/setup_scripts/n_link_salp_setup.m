
addpath("..\salp_modelling\")  % Add path of the salp modelling and animation

%% System properties for the salp

num_of_links = 3;

drag_ratio = 3;
linklengths = ones([num_of_links, 1]);
f_jets = cell(num_of_links, 1);
f_jets(:) = {zeros([1,3])};
f_jets{3} = [0 10 0];

salp = setupSalp(linklengths, drag_ratio, f_jets);

%% Setup solver configuration

g_init = zeros(1, 3);  % Initial position in world frame
r_init = zeros([1, num_of_links - 1]) + 0.5;  % Initial shape of the system
fixed_f_ext = zeros([3 + (num_of_links - 1), 1]);  % Precomputed constant additional body forces
fixed_f_jets = [salp.f_jets{:}]';  % Precomputed constant jet forces

solver_config = struct();
solver_config.dt = 1/30;  % timestep for simulation (30 fps)
solver_config.t_span = [0, 4];
solver_config.init_state = [g_init r_init];
solver_config.jet_force_func = @(t, X) fixed_f_jets;  % This has to be consistent to how jet forces are generated in local_force_func
solver_config.base_force_func = @(t, X) fixed_f_ext;
solver_config.local_force_func = @(t, X) fixed_f_jets;  % Forces that are applied locally on the links

%% Animation configuration

anim_config = struct();

anim_config.dir_loc = 'C:\Users\Yousef\Documents\University\Research\Misc scripts\salp_animation\animations';
anim_config.gif_name = 'link_5_3_y_10_N'; 
anim_config.tick_colors = [.8,.8,.8];  % Grey
anim_config.link_colors = [0.0745,0.6235,1.0000];  % Light Blue
anim_config.thruster_colors = [0.4941,0.1843,0.5569];  % Purple
anim_config.in_jet_colors = [1.0000,1.0000,0.0667];  % Yellow
anim_config.out_jet_colors = [1.0000,0.4118,0.1608];  % Orange

animateDistribution_discrete(salp, solver_config, anim_config);
