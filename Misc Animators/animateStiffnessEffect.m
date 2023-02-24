gifTime = 10;
FPS = 30;
stiffPhase = linspace(0,2*pi,gifTime*FPS);
maxCurv = 1;
bump = .05;
curvs = maxCurv/2 + maxCurv/2*cos(stiffPhase) + bump;

numLinks = 30;
blu = [0.0745    0.6235    1.0000];

lw = .1;
phaseDelay = pi/3;

L = .25;
a = L/2;
AR = .2;
ellipseThetas = linspace(0,2*pi,100);
ellipsePoints = [a*cos(ellipseThetas);a*AR*sin(ellipseThetas)];

fig10 = figure(10);
set(fig10,'Color',[1,1,1]);
for i = 1:numel(curvs)
    clf;
    hold on;

    maxCurvature = curvs(i);
    thisTheta = 0;
    curvPhase = 0;
    orientation = 0;
    linkPosition = [0;0];

    stiffness = maxCurv + 2*bump - maxCurvature;
    normStiffness = 6*stiffness/(maxCurv + 2*bump);
    box = [0,normStiffness,normStiffness,0;6.5,6.5,7.5,7.5];
    boxOutline = [0,6,6,0,0;6.5,6.5,7.5,7.5,6.5];
    fill(box(1,:),box(2,:),[1,0,0]);
    plot(boxOutline(1,:),boxOutline(2,:),'k','LineWidth',2);
    text(3,8,'Joint Stiffness','FontSize',12,'FontWeight','bold','HorizontalAlignment','center')

    thisEllipse = ellipsePoints;

    for j = 1:numLinks
        fill(thisEllipse(1,:),thisEllipse(2,:),blu);

        linkPosition = linkPosition + [a*cos(orientation);a*sin(orientation)];

        curvPhase = curvPhase + phaseDelay;
        thisTheta = -maxCurvature/2*cos(curvPhase) + maxCurvature/2;
        orientation = orientation + thisTheta;

        linkPosition = linkPosition + [a*cos(orientation);a*sin(orientation)];
        thisEllipse = shiftPoints(ellipsePoints,orientation,linkPosition);
    end

    axis([-2.5,numLinks*L+.5,-.5,9]);

    if i == 1
        gif('StiffnessIllustration.gif','DelayTime',1/FPS);
    else
        gif;
    end

end

        


%Generates rotation matrix to rotate 2D points
function R = getRotMatrix(theta)

    R = [cos(theta),-sin(theta);sin(theta),cos(theta)];

end


%Rotates and translates datapoints
%points - [xpoints;ypoints]
%theta - angle to rotate by
%dpoint - [dx;dy]
function newPoints = shiftPoints(points,theta,dpoints)

    R = getRotMatrix(theta);
    
    rotatedPoints = R*points;
    
    newPoints = rotatedPoints + dpoints;

end
