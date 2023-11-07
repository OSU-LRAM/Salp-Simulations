
addpath("..\salp_modelling\")  % Add path of the salp modelling and animation
load("salp_3_link_basis_11_t_angle_30.mat")

%% Extract distribution

B_aug = cell(1, numel(poss_basis));
for idx = 1:numel(poss_basis)
    B_aug{idx} = poss_basis{idx}.field;
end

%% System properties for the salp

n = 3;  % Number of links
num_of_joints = n - 1;
drag_ratio = 3;
linklengths = ones([n, 1]);
f_jets = cell(n, 1);
f_jets(:) = {zeros([1,3])};
f_jets{3} = [0 10 0];

salp = setupSalp(linklengths, drag_ratio, f_jets);

%% Jet force parameters

deg_of_freedom = 2;

T = .1;
num_of_cycles = 20;
omega = 2*pi/T;
thrust_angle = pi/6;
max_force_thrust = 50;
bias_constant_jets = true;
q_dot_des = zeros([(n - 1) + 3, 1]);
q_dot_des(deg_of_freedom) = 1;

%% Setup solver configuration

g_init = zeros(1, 3);  % Initial position in world frame
r_init = zeros([1, n - 1]) + pi/8;  % Initial shape of the system
r_init = zeros(1, n-1);
q_init = [g_init r_init]';
fixed_f_ext = zeros([3 + (n - 1), 1]);  % Precomputed constant additional body 

solver_config = struct();
solver_config.dt = 1/30;  % timestep for simulation (30 fps)
solver_config.t_span = [0, num_of_cycles*T];
solver_config.init_state = [g_init r_init];
solver_config.base_force_func = @(t, X) fixed_f_ext;
solver_config.local_force_func = getShapeTrackForceFunction(num_of_joints, r_init, q_dot_des, omega, B_aug, domain, thrust_angle);

%% Animation configuration

anim_config = struct();

anim_config.dir_loc = 'C:\Users\Yousef\Documents\University\Research\Misc scripts\salp_animation\animations';

if bias_constant_jets
    bias_str = 'y';
else
    bias_str = 'n';
end

if deg_of_freedom == 1
    deg_of_freedom_str = 'x';
elseif deg_of_freedom == 2
    deg_of_freedom_str = 'y';
elseif deg_of_freedom == 3
    deg_of_freedom_str = 'theta';
end

anim_config.gif_name = ['link_' num2str(2) '_t_angle_' num2str((180/pi)*thrust_angle) 'period_' num2str(T) '_move_' deg_of_freedom_str '_bias_' bias_str]; 
anim_config.tick_colors = [.8,.8,.8];  % Grey
anim_config.link_colors = [0.0745,0.6235,1.0000];  % Light Blue
anim_config.thruster_colors = [0.4941,0.1843,0.5569];  % Purple
anim_config.in_jet_colors = [1.0000,1.0000,0.0667];  % Yellow
anim_config.out_jet_colors = [1.0000,0.4118,0.1608];  % Orange

animateDistribution_discrete(salp, solver_config, anim_config);
