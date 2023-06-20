%Generates 3 individual joint stability plots (one for thrust from each link in the salp)
%given a salp definition (i.e. joint stiffnesses, drag ratio) where the 
%salp is continuously thrusting.  Also generate four joint trajectories 
%from different joint initial conditions on this stability plot

function plotStruct = makeStabilityDiagram(salp,startPoints)

    %Grab salp definitions
    L = salp.L;
    dragRatio = salp.dragRatio;
    jointStiffness = salp.jointStiffness;
    nJoints = 2;

    %Make drag matrix for the links
    d = diag([L,dragRatio*L,dragRatio*L^3/12]);
    %Make stiffness matrix for system
    K = [zeros(3,nJoints);diag(-jointStiffness*ones(1,nJoints))];

    %Generate grid of joint values
    dens = 11;
    r1 = linspace(-pi,pi,dens);
    r2 = linspace(-pi,pi,dens);
    [R1,R2] = meshgrid(r1,r2);

    %Hold plot structure for thrust coming from each of the 3 links
    plotStruct = cell(1,3);

    %For thrust coming from each link
    for link = 1:3
        %Get thrust force model
        thrusterForces = zeros(9,1);
        thrusterForces((link-1)*3+1:(link-1)*3+3) = [1;0;0];
        thrusterForces = @(t) thrusterForces;

        %Set up holder for resultant joint velocities
        dR1 = zeros(size(R1));
        dR2 = zeros(size(R2));

        %At each joint angle for this thrust/stiffness/drag ratio, produce
        %resultant speeds
        for i = 1:numel(R1)
            speeds = thrustToSpeeds(0,[R1(i),R2(i)],thrusterForces,L,d,K);
            dR1(i) = speeds(4);
            dR2(i) = speeds(5);
        end

        %Store vector field plot data
        plotStruct{link} = struct();
        plotStruct{link}.R1 = R1;
        plotStruct{link}.R2 = R2;
        plotStruct{link}.dR1 = dR1;
        plotStruct{link}.dR2 = dR2;

        %Get the path evolution for each joint initial condition
        setvals = 0;
        for sp = 1:size(startPoints,1)
            path = getFullPath(salp,thrusterForces,startPoints(sp,:));
            if ~setvals
                setvals = 1;
                plotStruct{link}.pathx = zeros(size(startPoints,1),numel(path(1,:)));
                plotStruct{link}.pathy = zeros(size(startPoints,1),numel(path(2,:)));
            end
            plotStruct{link}.pathx(sp,:) = path(1,:);
            plotStruct{link}.pathy(sp,:) = path(2,:);
        end

end