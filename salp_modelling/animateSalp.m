
function animateSalp(salp, simul_results, anim_config)
%% General code to generate the salp animation

    % Initialize objects that will be used in loop 
    figure(10);
    f10 = gcf;
    base_anim_geom = generateBaseSalp(salp, simul_results, anim_config);  % Base geometry

    % Get geometry per time step
    for idx = 1:numel(simul_results.ts)

        % Make figure a blank slate
        clf;
        axis equal;
        hold on;

        % Get link transformations for current timestep
        link_trans_mats = getLinkTransformations(salp, simul_results.states(:, idx), anim_config.baseframe);

        % Generate geometry for current timestep
        plotGrid(base_anim_geom);
        plotCurrentSalpGeometry(salp, base_anim_geom, link_trans_mats, simul_results.forces, idx)

        % Center grid on salp
        centerGridOnSalp(salp, link_trans_mats);

        % Set figure position on screen
        f10.Position = [100,120,900,800];

        % Make sure everything has finished drawing
        drawnow;

        % Save frame to gif
        saveAnimationFrame(idx, simul_results, anim_config)

    end
end


function centerGridOnSalp(salp, link_trans_mats)

    % Calculate centroid of salp for choosing axis limits
    sum_position = zeros(2, 1);
    for idx = 1:salp.num_of_links
        sum_position = sum_position + link_trans_mats{idx}(1:2,3);
    end

    nLinks = salp.num_of_links;
    salp_centroid = sum_position/salp.num_of_links;

    % Provide buffer in all directions around c.g.

    L = mean(salp.linklengths);
 
    minX = salp_centroid(1)-nLinks*L*2/3;
    maxX = salp_centroid(1)+nLinks*L*2/3;
    minY = salp_centroid(2)-nLinks*L*2/3;
    maxY = salp_centroid(2)+nLinks*L*2/3;

    % Fit axis so we're always zoomed on the swimmer
    axis([minX,maxX,minY,maxY]);
end


function saveAnimationFrame(ts, simul_results, anim_config)
% Save the frame for the salp animation
    
    % If this is the first frame
    if ts == 1
        % Start the gif
        try
            if isempty(anim_config.dir_loc)
                file_delim = '';
            else
                file_delim = '\';
            end
            gif_loc = strcat(anim_config.dir_loc, file_delim, anim_config.gif_name, '.gif');
            gif(gif_loc, 'DelayTime', simul_results.dt);
        catch
            error('Animation requires "gif" package to be installed');
        end
    else
        % Add frame to gif
        gif;
    end
end


function link_trans_mats = getLinkTransformations(salp, state, anim_baseframe)
% Get link transformations to transform geometry (assumes link base is in
% tail to simplify how links are transformed)
    
    % Storage for transformations
    link_trans_mats = cell(salp.num_of_links, 1);

    % Initialize first transformation 
    world_pos = state(1:2);
    world_angle = state(3);
    link_trans_mats{1} = getAffineTransformation(world_angle, world_pos, [1 1]);
    
    % Get other link transformations
    for idx = 1:salp.num_of_links - 1

        % Get relevant angle and link lengths
        alpha = state(3 + idx);
        curr_len = salp.linklengths(idx);
        next_len = salp.linklengths(idx + 1);

        % Get respective SE2 elements
        M1 = getAffineTransformation(alpha, [curr_len/2 0], [1 1]);
        M2 = getAffineTransformation(0, [next_len/2 0], [1 1]);

        % Save link-to-link transformation in cell
        link_trans_mats{idx + 1} = link_trans_mats{idx}*M1*M2;     
    end

    % Get rebasing matrix and rebase matrices
    if ~strcmp(salp.baseframe, anim_baseframe)
        salp.baseframe = anim_baseframe;
        [~, ~, ~, frame_0] = N_link_chain(salp, state(4:end));
        anim_rebase_frame = inv(frame_0);
        for idx = 1:salp.num_of_links
            link_trans_mats{idx} = anim_rebase_frame * link_trans_mats{idx};
        end
    end
    


end


function trans_mat = getAffineTransformation(theta, shift, scale)
% Get matrix representation for the given transformation parameters
%
% Inputs:
%
%   points: A set of points to transform
%
%   theta: Angle in which to rotate points
%
%   trans: Translation vector
%
% Outputs:
%
%   new_points: New set of transformed points

    % Rotation matrix
    R = [cos(theta),-sin(theta), 0;
         sin(theta), cos(theta), 0;
         0         , 0         , 1];

    % Translation matrix
    T = eye(3);
    T(1:2, 3) = shift;
    
    % Scale matrix
    S = diag([scale(:)' 1]);

    trans_mat = T*R*S;
end


function plotGrid(base_anim_geom)
% Plot grid marks to distinguish salp motion

    ticks = base_anim_geom.ticks;

    for idx = 1:numel(ticks.points{1})
        x_ticks = ticks.points{1};
        plot([x_ticks(idx),x_ticks(idx)], ...
            [ticks.plot_min{2},ticks.plot_max{2}], ...
            'Color', ticks.color);
    end

    for idx = 1:numel(ticks.points{2})
        y_ticks = ticks.points{2};
        plot([ticks.plot_min{1},ticks.plot_max{1}], ...
            [y_ticks(idx),y_ticks(idx)], ...
            'Color', ticks.color);
    end
end


function plotCurrentSalpGeometry(salp, base_anim_geom, link_trans_mat, forces, ts)
    
    % Get base geometry components
    salp_comps = {
        base_anim_geom.links;
        base_anim_geom.thrusters;
        base_anim_geom.in_force_jets;
        base_anim_geom.out_force_jets
        };

    % Fill in each of the salp components into a frame
    for idx = 1:numel(salp_comps)
        salp_comp = salp_comps{idx};


        for idx2 = 1:salp.num_of_links

            % Get relevant info to shift points
            M = link_trans_mat{idx2};
            points = salp_comp.points{idx2};

            % If component is not links, rotate component by thruster angle
            if ~isequaln(salp_comp, base_anim_geom.links)
                forces.thruster_angles{ts}(idx2)
                M = M*getAffineTransformation(forces.thruster_angles{ts}(idx2), [0, 0], [1, 1]);
            end

            % If component is one of the forces, adjust to current force
            % output
            if isequaln(salp_comp, base_anim_geom.in_force_jets) || isequaln(salp_comp, base_anim_geom.out_force_jets)
                points = points * forces.jets(idx2, ts);
            end

            % Shift points to correct place
            shifted_points = shiftPoints(M, points);
            fill(shifted_points(1,:), shifted_points(2,:), salp_comp.color);
        end
    end
end


function shifted_points = shiftPoints(trans_mat, points)
    
    move = trans_mat(1:2, 3);
    rot_mat = trans_mat(1:2, 1:2);

    shifted_points = rot_mat*points + move;
end


function base_anim_geom = generateBaseSalp(salp, simul_results, anim_config)
% Create base link geometry for the salp
% It creates each segment of the salp based on the length of the link

    % Create base salp geometry storage
    base_anim_geom = struct();
    base_anim_geom.ticks = createBaseSalpComponent(salp, anim_config.tick_colors);
    base_anim_geom.links = createBaseSalpComponent(salp, anim_config.link_colors);
    base_anim_geom.thrusters = createBaseSalpComponent(salp, anim_config.thruster_colors);
    base_anim_geom.in_force_jets = createBaseSalpComponent(salp, anim_config.in_jet_colors);
    base_anim_geom.out_force_jets = createBaseSalpComponent(salp, anim_config.out_jet_colors);

    %% Set points for system
   
    % Base grid to draw background lines that will help judge swimmer 
    % motion
    ticks = cell(2, 1);
    plotMin = ticks;
    plotMax = plotMin;
    backgroundWidth = 1;
    for idx = 1:2
        plotMin{idx} = min(simul_results.states(idx,:)) - salp.num_of_links*mean(salp.linklengths);
        plotMax{idx} = max(simul_results.states(idx,:)) + salp.num_of_links*mean(salp.linklengths);
        ticks{idx} = [-1*fliplr(0:backgroundWidth:abs(plotMin{idx})),backgroundWidth:backgroundWidth:abs(plotMax{idx})];
        
    end

    base_anim_geom.ticks.points = ticks;
    base_anim_geom.ticks.plot_min = plotMin;
    base_anim_geom.ticks.plot_max = plotMax;

    % Create overall base n-link geometry
    jet_points = cell(numel(salp.linklengths), 1);
    for idx = 1:salp.num_of_links
        
        L = salp.linklengths(idx);  % Current link length
      
        % Build link geometry
        major_axis = L/2;
        minor_axis = L/10;
        ellipse_thetas = linspace(0, 2*pi, 100);
        base_anim_geom.links.points{idx} = [...
            major_axis*cos(ellipse_thetas); ...
            minor_axis*sin(ellipse_thetas)];

        % Build thruster geometry
        thruster_radius = L/6;
        thruster_thetas = linspace(-pi/2,pi/2,50);
        base_anim_geom.thrusters.points{idx} = [...
            thruster_radius*cos(thruster_thetas);...
            thruster_radius*sin(thruster_thetas)];

        % Build force jet geometry
        thrust_jet_radius = thruster_radius*.8;
        max_thrust_length = 2/3*L;
        thrust_jet_thetas = linspace(pi/2,3*pi/2,50);
        jet_points{idx} = [...
            max_thrust_length*cos(thrust_jet_thetas);...
            thrust_jet_radius*sin(thrust_jet_thetas)];

        base_anim_geom.in_force_jets.points{idx} = jet_points{idx};
        base_anim_geom.out_force_jets.points{idx} = .5*jet_points{idx};
    end
    
end 


function comp = createBaseSalpComponent(salp, comp_color)
% Create structure to store salp components
    
    comp = struct();
    comp.color = comp_color;
    comp.points = cell(salp.num_of_links, 1);
end