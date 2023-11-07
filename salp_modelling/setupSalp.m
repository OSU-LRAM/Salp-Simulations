function salp = setupSalp(linklengths, drag_ratio, f_jets)
% Shorthand to setup the salp object
%
% Inputs:
%
%   linklengths: Array of link lengths for the salp system
%
%   drag_ratio: A number representing the drag ratio for the salp system
%
%   f_jets: Cell array with arrays that represent the force components of 
%       the jet forces. This is used to specify what distribution vector
%       fields you want to plot.
%
% Outputs:
%
%   salp: Structure representing the salp system

    salp = struct();

    salp.baseframe = 'tail';  % Setup by default to simplify animation generation
    salp.drag_ratio = drag_ratio;
    salp.linklengths = linklengths;
    salp.num_of_links = size(salp.linklengths, 1);
    salp.num_of_joints = salp.num_of_links - 1;

end

