%Gets a joint path evolution given a certain initial condition and salp
%definition
function path = getFullPath(salp,Ft,initPoint)
    %Grab salp definitions
    L = salp.L;
    dragRatio = salp.dragRatio;
    jointStiffness = salp.jointStiffness;
    nJoints = 2;

    %Make drag matrix for the links
    d = diag([L,dragRatio*L,dragRatio*L^3/12]);
    %Make stiffness matrix for system
    K = [zeros(3,nJoints);diag(-jointStiffness*ones(1,nJoints))];

    %Set up simulation parameters
    dt = 1/50;
    T = 5;
    tSpan = [0,T];
    ts = [0:dt:T];

    %Simulate joint path evolution
    odefun = @(t,X) jointODE(t,X,Ft,L,d,K);
    sol = ode45(odefun,tSpan,initPoint);
    path = deval(sol,ts);

end
    
%ODE describing joint path evolution for a certain salp definition
function dX = jointODE(t,X,Ft,L,d,K)

    r = X;

    %Calculate swimmer body/shape velocity from thrust at time t
    speeds = thrustToSpeeds(0,r,Ft,L,d,K);

    %Grab shape velocities
    dX = speeds(4:end);

end