function plotDistribution_discrete(salp, graph_config)
    
    %% Create plotting  elements

    % Create shape domain to plot in
    shape_domain = createShapeDomain(salp, graph_config);

    % Check which force components are nonzero (if it was requested)
    f_idxs = {};
    if isfield(graph_config, "skip_zeros") && graph_config.skip_zeros
        f_idxs = checkForNonZeroForces(salp);
    end

    % Generate distribution vector fields to plot
    vect_fields = generateDistributionVectorFields(...
        salp, ...
        graph_config, ...
        shape_domain, ...
        f_idxs ...
        );

    %% Generate plots
    
    if isfield(vect_fields, "Bg")
        titles = createGraphTitles(salp, "g", f_idxs);
        generateDistributionPlots(317, salp, graph_config, vect_fields.Bg, shape_domain, titles);
    end
    if isfield(vect_fields, "Br")
        shape_size = numel(salp.f_jets) - 1;
        if shape_size == 2
            shape_domain = shape_domain(:, 1:2);
        end
        titles = createGraphTitles(salp, "r", f_idxs);
        generateDistributionPlots(318, salp, graph_config, vect_fields.Br, shape_domain, titles);
    end

end


function vect_fields = generateDistributionVectorFields(salp, graph_config, shape_domain, f_idxs)
% Generate distribution vector fields for a given shape domain

    %% Extract relevant information from salp system
    if isempty(f_idxs)
        num_of_fields = 3*numel(salp.f_jets);
    else
        num_of_fields = numel(f_idxs);
    end
    
    %% Create storage for vector fields
    vect_fields = struct();

    % Create empty Bg vector field if it was requested
    if isfield(graph_config, "g") && graph_config.g
        vect_fields.Bg = createEmptyVectorFields( ...
            num_of_fields, ...
            [3, size(shape_domain, 1)] ...
            );
    end

    % Create empty Br vector field if it was requested
    if isfield(graph_config, "r") && graph_config.r
        shape_size = numel(salp.f_jets) - 1;
        vect_fields.Br = createEmptyVectorFields( ...
            num_of_fields, ...
            [shape_size, size(shape_domain, 1)] ...
            );
    end

    %% Calculate distribution and populate vector fields
    for idx = 1:size(shape_domain, 1)
        shape = shape_domain(idx, 1:2);
        B = getDistribution_discrete(salp, shape, f_idxs);

        % Fill out respective vector fields
        for idx2 = 1:num_of_fields

            % Build Br vector fields if it was requested
            if isfield(vect_fields, "Br")
                vect_fields.Br{idx2}(:, idx) = B(4:end, idx2);
            end

            % Build Bg vector fields if it was requested
            if isfield(vect_fields, "Bg")
                vect_fields.Bg{idx2}(:, idx) = B(1:3, idx2);
            end 
        end
    end

end


function generateDistributionPlots(fignum, salp, graph_config, vect_fields, shape_domain, titles)

    is_3d = size(shape_domain, 2) == 3;

    % Generate axes to plot in

    if isfield(graph_config, "skip_zeros") && graph_config.skip_zeros
        m = 1;
        n = numel(vect_fields);
        num_of_fields = n;
    else
        m = 3;
        n = numel(salp.f_jets);
        num_of_fields = 3*n;
    end
    ax = create_subaxes(fignum, m, n, num_of_fields);

    % Create shape domain
    X = shape_domain(:, 1)';
    Y = shape_domain(:, 2)';
    if is_3d
        Z = shape_domain(:, 3)';
    end
    
    % Plot the vector fields
    for idx = 1:num_of_fields
        U = vect_fields{idx}(1, :);
        V = vect_fields{idx}(2, :);
        if is_3d
            W = vect_fields{idx}(3, :);
            quiver3(X, Y, Z, U, V, W, ...
                "Parent", ax{idx}, ...
                "Color", "black");
        else
            quiver(X, Y, U, V, ...
                "Parent", ax{idx}, ...
                "Color", "black");
        end

        % Set plot properties
        axis(ax{idx}, "equal");
        axis(ax{idx}, "square");
        box(ax{idx}, "on");
        ax{idx}.Title.String = titles{idx};
        ax{idx}.XLabel.String = "\alpha_{1}";
        ax{idx}.YLabel.String = "\alpha_{2}";
        
    end
    
end


function graph_titles = createGraphTitles(salp, var_str, f_idxs)
% Helper function to create graph titles for the distribution fields

    var_comp_str = ["x", "y", "\theta"];

    % Get the number of fields
    if isempty(f_idxs)
        num_of_fields = 3*numel(salp.f_jets);
    else
        num_of_fields = numel(f_idxs);
    end
    graph_titles = cell(num_of_fields, 1);


    % Create titles
    for idx = 1:num_of_fields

        % Choose what strings to use
        if isempty(f_idxs)
            col_idx = fix((idx - 1)/3) + 1;
            comp_idx = mod((idx - 1), 3) + 1;
        else
            col_idx = f_idxs{idx}(1);
            comp_idx = f_idxs{idx}(2);
        end

        % Save title
        graph_titles{idx} = var_str + "_{" + ...
                + col_idx + var_comp_str(comp_idx) + ... 
                "}";
    end
    
end


function B = getDistribution_discrete(salp, shape, f_idxs)
% Get distribution matrix for a set of forces
%
% Inputs:
%
%   geom: Geometry struct 
%
% Outputs:
%
%   B: Distribution matrix for a given set of Jacobians

    [~,~,J_full] = N_link_chain(salp, shape);
    J_dual_full = getDualJacobians_discrete(J_full);
    M = getDragMetric(salp.linklengths, salp.drag_ratio, J_full);

    B = -M\J_dual_full;

    % Choose rows if f_idxs is not empty
    if ~isempty(f_idxs)
        B_new = zeros(size(B, 1), numel(f_idxs));
        f_jet_size = numel(salp.f_jets{1});
        for idx = 1:numel(f_idxs)
            col_idx = f_jet_size*(f_idxs{idx}(1) - 1) + f_idxs{idx}(2);
            B_new(:, idx) = B(:, col_idx);
        end
        B = B_new;
    end

end

function f_idxs = checkForNonZeroForces(salp)
% Save force jet components that are not zero
    
    f_idxs = cell(0);

    for idx = 1:numel(salp.f_jets)
        for idx2 = 1:numel(salp.f_jets{idx})

            % Save index if force component is not zero
            if salp.f_jets{idx}(idx2) ~= 0
                f_idxs{end + 1} = [idx, idx2];
            end

        end
    end

end


function vect_fields = createEmptyVectorFields(num_of_fields, field_size)
% Generate empty matrices in cell array to store multiple vector fields

    vect_fields = cell(1, num_of_fields);
    for idx = 1:num_of_fields
        vect_fields{idx} = zeros(field_size);
    end
end


function [domain, domain_size] = createShapeDomain(salp, graph_config)
% Create shape domain based on the system and graph configuration

    shape_size = salp.num_of_links;

    % Extract relevant configuration
    delta = graph_config.delta;
    range = graph_config.range;
    
    % Create XY plane
    points = range(1):delta:range(2);
    [X, Y] = meshgrid(points);
    domain_size = numel(X);

    % Create domain array
    if isfield(graph_config, "g") || (isfield(graph_config, "r") && shape_size == 3 && salp.r)
        domain = [X(:) Y(:) zeros(numel(X), 1)];
    else
        domain = [X(:) Y(:)];
    end

end


function [ax,f] = create_subaxes(fignum, m, n, p)
% Clear out a specified figure and create a clean set of axes in that
% figure with equal-axis aspect ratio
%
% Inputs:
%
%   fignum: The number of the figure (or a figure handle) in which to
%       create the axes 
%   m: the number of rows of subplots to create
%   n: the number of columns of subplots to create
%   p: the number of subplots to create (should be less than or equal to
%       m*n)
%
% Outputs:
%
%   ax: A cell array of handles to the created subplot axes
%   f: A handle to the figure that was created

    f = figure(fignum);   
    clf(f, "reset");
    ax = cell(1, p);
    
    %%%%%%%
    % Loop over the number of subplots requested
    for idx = 1: p
        ax{idx} = subplot(m, n, idx, "Parent", f);
    end
end
