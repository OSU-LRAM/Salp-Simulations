tic;

%Remember to initialize sysplotter and have gif Add-On installed

%Declare salp characteristic values
salp = struct();
%Length of individual zooid units
salp.L = .25;
%Ratio of Y drag to X drag of a link - Y drag probably higher
salp.dragRatio = 3;
%Stiffness of springs connecting links
salp.jointStiffness = .15;
%Value of maximum thrust available from zooids
salp.thrustAmp = 2;
%Angle of thrusters relative to back of link
salp.thrustAngle = pi/12;
%Phase lag between link and its neighbors for synchronized swimming
salp.phaseDelay = pi/4;
%Frequency at which each thruster operates
salp.thrustFrequency = 1;
%Number of periods to simulate
salp.numPeriods = 8;
%Percent of time thrusters are activated
salp.percentThrusting = .75;
%Type of thruster activation: can be 'synch', 'asynch', or 'random'
salp.forceType = 'synch';
%Frequency at which to animate salp
salp.FPS = 60;

%Run salp simulation
[ts,states,forces] = integrateSalpMotion(12,salp);
toc

%animate the results
animateSalp(ts,states,forces,salp);