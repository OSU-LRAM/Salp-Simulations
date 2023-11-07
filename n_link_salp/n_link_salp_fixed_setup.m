addpath("..\salp_modelling\")  % Add path of the salp modelling and animation

%% System properties for the salp

n = 3; % Number of links in salp
salp = struct();

salp.baseframe = 'tail';  % Setup by default to simplify animation generation
salp.drag_ratio = 3;
salp.linklengths = ones([n, 1]);
salp.num_of_links = size(salp.linklengths, 1);
salp.num_of_joints = salp.num_of_links - 1;

%% Solver parameters

T = .1;
omega = 2*pi/T;
num_of_cycles = 20;
thrust_angle = pi/6;

%% Setup solver configuration
g_init = zeros(1, 3);  % Initial position in world frame (not Lie algebra element g_circ_left)
r_init = zeros(1, n-1);  % Initial joint angles of salp (initial shape)
q_init = [g_init r_init]';

solver_config = struct();
solver_config.dt = 1/30;  % timestep for simulation (30 fps)
solver_config.t_span = [0, num_of_cycles*T];
solver_config.init_state = [g_init r_init];
solver_config.base_force_func = @(t, X) zeros([3 + (n - 1), 1]);
solver_config.local_force_func = getLocalFixedForceFunction([1, 2, 3], thrust_angle);

%% Animation configuration

anim_config = struct();

anim_config.baseframe = 'com-mean'; % What local frame to use for the animation

% Gif directory location and name
anim_config.dir_loc = '';  % You can put a directory path here to save the gif
anim_config.gif_name = 'example_animation'; 

% Colors for salp in the animation
anim_config.tick_colors = [.8,.8,.8];  % Grey
anim_config.link_colors = [0.0745,0.6235,1.0000];  % Light Blue
anim_config.thruster_colors = [0.4941,0.1843,0.5569];  % Purple
anim_config.in_jet_colors = [1.0000,1.0000,0.0667];  % Yellow
anim_config.out_jet_colors = [1.0000,0.4118,0.1608];  % Orange

animateDistribution_discrete(salp, solver_config, anim_config);
