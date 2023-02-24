%integrates motion of a salp-like swimmer over time

%out:
%ts - time vector over which output takes place
%states - matrix of states (body positions/orientations and joint angles)
%over time
%forces - matrix of input force vectors over time

%in:
%nLinks - number of links composing the salp
%getSteadyState - boolean determining whether only the end steady state
%results are determinend vs the evolution over time

function [ts,states,forces] = integrateSalpMotion(nLinks,salp,getSteadyState)

    %Set defaults for variables in case function is called with few inputs
    if ~exist('getSteadyState','var')
        %Default return only the limit cycle evaluation
        getSteadyState = 0;
    end

    %Can be 'synch', 'asynch', or 'random'
    %synch refers to synchronized thrusters with phase offset
    %asynch refers to asynchronized thrusters firing constantly
    %random refers to asynchronized thrusters with individual stochastic
    %phase delays
    forceType = salp.forceType;
    %Get thruster frequency
    frequency = salp.thrustFrequency;
    %Number of joints in system
    nJoints = nLinks - 1;
    %Zero vector initial condition for swimmer ODE
    initialState = zeros(3+nJoints,1);
    %Evaluate 5 time periods before finding limit cycle
    numPeriods = salp.numPeriods;
    %Time period of thruster actuation
    T = 1/frequency;
    %Time span over which swimmer ODE will be evaluated
    tSpan = [0,numPeriods*T];
    %Sample output at animation frequency for ease of animation
    FPS = salp.FPS;
    dt = 1/FPS;
    %Generate time vector for output
    ts = [0:dt:tSpan(end)];

    %Set constants for swimmer
    %Length is a quarter-meter
    L = salp.L;
    %Y-drag coefficient is dragRatio times X-drag coefficient
    dragRatio = salp.dragRatio;
    %Stiffness of joint springs
    jointStiffness = salp.jointStiffness;

    %Make drag matrix for the links
    d = diag([L,dragRatio*L,dragRatio*L^3/12]);
    %Make stiffness matrix for system
    K = [zeros(3,nJoints);diag(-jointStiffness*ones(1,nJoints))];

    %Force function descriptions
    %Magnitude of max thrust in N
    thrustAmp = salp.thrustAmp;
    %Angle of thrust with respect to rear of link
    thrustAngle = salp.thrustAngle;
    %How much current link lags behind distal link in thruster phase
    phaseDelay = salp.phaseDelay;

    %Generate forcing function fpr evaluation at generic time
    %If thrusters are synchronized
    if strcmp(forceType,'synch')
        %Phase delays start at 0 for head link and stack down to tail
        phaseDelays = phaseDelay*[0:nLinks-1]';
        forcingFunction = generateForceVectorFunctions(nLinks,thrustAmp,frequency,thrustAngle,phaseDelays);
    %If thrusters are not synchronized
    elseif strcmp(forceType,'asynch')
        %Phase delays are all random
        phaseDelays = 2*pi*[rand(nLinks,1)];
        forcingFunction = generateForceVectorFunctions(nLinks,thrustAmp,frequency,thrustAngle,phaseDelays);
    %If there is a randomized delay between thruster fires
    elseif strcmp(forceType,'random')
        forcingFunction = generateDispersedForceVectorFunctions(nLinks,thrustAmp,frequency,thrustAngle,numPeriods,salp.percentThrusting);
    end
    
    %Set ODE description for swimming system
    %sets (forcingFunction,L,d,K) for ODE, changes to them after
    %this point will not be seen in ODE results
    odefun = @(t,X) salpODE(t,X,forcingFunction,L,d,K);
    %integrate swimmer ODE over desired time period given initial
    %conditions
    sol = ode45(odefun,tSpan,initialState);
    %If we want the limit cycle results
    if getSteadyState
        %Take the swimmer state at the end of the initial cycles and set it
        %as the initial condition for the limit cycle
        final_loop_start = sol.y(:,end);
        %Move swimmer back to origin
        final_loop_start(1:3) = [0;0;0];
        %Calculate results from limit cycle initial condition
        sol = ode45(odefun,[0,T],final_loop_start);

        %Replace time vector with time vector on limit cycle
        ts = [0:dt:T];
    end

    %Evaluate solution to ODE at desired time steps for ease of animation
    states = deval(sol,ts);
    %Get sample value of forcing function so we know dimensionality
    f0 = forcingFunction(0);
    %Initialize forces as zero matrix
    forces = zeros(numel(f0),numel(ts));
    for i = 1:numel(ts)
        %Evaluate forcing function across desired time interval
        forces(:,i) = forcingFunction(ts(i));
    end

end

%Describes ODE for salp swimmer

%out:
%dX - state (X) rate of change

%in:
%t - time at which to evaluate system
%X - system state, X = [bodyPosition;jointAngles]
%Ft - forcing function (anonymous function that can be evaluated at time t)
%describing forces acting on each link due to exerted thrust
%L - length of each link
%d - drag matrix for each link
%K - stiffness matrix such that spring contribution to body forces
% [Fbody;Fshape] = K*r, for body shape r
function dX = salpODE(t,X,Ft,L,d,K)

    %Decompose system state into system position and shape variables
    g_vec = X(1:3);
    r = X(4:end);

    %Calculate swimmer body/shape velocity from thrust at time t
    speeds = thrustToSpeeds(t,r,Ft,L,d,K);

    %Decompose velocity into body velocity and shape velocity
    g_circ = speeds(1:3);
    r_dot = speeds(4:end);

    %Find left action that takes body velocities to world velocities
    theta = X(3);
    g_mat = [cos(theta),-sin(theta),0;...
    sin(theta),cos(theta),0;...
    0,0,1];
    %Turn swimmer body velocity into swimmer world velocity
    g_dot = g_mat*g_circ;

    %Compose state rate of change
    %Position rate of change is swimmer world velocity
    %Shape rate of change is swimmer joint velocity
    dX = [g_dot;r_dot];

end
