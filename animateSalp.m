%Makes a gif illustrating salp motion from ODE simulation

%ts - time vector
%states - [world positions of tail link;joint angles] for each time in ts
%forces - vector of forces for each link at each time in ts
%salp - structure describing salp characteristic values
%doOtherPlot - boolean, decides whether to plot salp data synchronized with
%animations

function animateSalp(ts,states,forces,salp,doOtherPlot)

    if ~exist('doOtherPlot','var')
        doOtherPlot = 0;
    end

    %Get link length from system and store it
    L = salp.L;

    %Get thruster angles and calculate angle relative to front of each link
    %for left and right-handed thrusters
    sameSide = 0;

    thrustAngle = salp.thrustAngle;
    theta_left = pi-thrustAngle;
    if sameSide
        theta_right = theta_left;
    else
        theta_right = pi+thrustAngle;
    end

    %Major axis radius for salp link
    a = L/2;
    
    %Time step for animation
    dt = ts(2)-ts(1);
    %Get number of links in swimmer
    nLinks = size(states,1)-2;

    %Prep values at which to draw background lines that will help judge
    %swimmer motion
    backgroundWidth = 1;
    plotMinX = min(states(1,:)) - 1.5*nLinks*L;
    plotMaxX = max(states(1,:)) + 1.5*nLinks*L;
    plotMinY = min(states(2,:)) - 1.5*nLinks*L;
    plotMaxY = max(states(2,:)) + 1.5*nLinks*L;
    xTicks = [-1*fliplr([0:backgroundWidth:abs(plotMinX)]),[backgroundWidth:backgroundWidth:abs(plotMaxX)]];
    yTicks = [-1*fliplr([0:backgroundWidth:abs(plotMinY)]),[backgroundWidth:backgroundWidth:abs(plotMaxY)]];
    
    %Draw a single swimmer link that we'll place at each position in the
    %chain
    aspectRatio = 1/5;
    ellipseThetas = linspace(0,2*pi,100);
    ellipsePoints = [a*cos(ellipseThetas);a*aspectRatio*sin(ellipseThetas)];

    %Draw a thruster that we'll place on each link
    thrusterRadius = L/6;
    thrusterThetas = linspace(pi/2,3*pi/2,50);
    thrusterPoints = [thrusterRadius*cos(thrusterThetas);thrusterRadius*sin(thrusterThetas)];

    %Draw the force jet that we'll place on each thruster
    thrustForceRadius = thrusterRadius*.8;
    maxThrustLength = 2/3*L;
    thrustForceThetas = linspace(3*pi/2,5*pi/2,50);
    outForce = [maxThrustLength*cos(thrustForceThetas);thrustForceRadius*sin(thrustForceThetas)];
    
    %Useful colors for testing animation
    black = [0,0,0]; %Black
    purple = [0.4941,0.1843,0.5569]; %Purple
    darkBlue = [0,0.4471,0.7412]; %Dark Blue
    lightBlue = [0.0745,0.6235,1.0000]; %Light Blue
    orange = [1.0000,0.4118,0.1608]; %Orange
    yellow = [1.0000,1.0000,0.0667]; %Yellow
    grey = [.8,.8,.8]; %Grey

    %Set animation colors
    linkColor = lightBlue;
    thrusterColor = purple;
    outForceColor = orange;
    inForceColor = yellow;
    backgroundTicksColor = grey;
    
    maxCurv = max(max(abs(states(4:end,:))));

    %Initialize figure
    fig10 = figure(10);
    set(fig10,'Color',[1,1,1]);
    %For every time instance
    for i = 1:numel(ts)
    
        %Make figure a blank slate
        clf;
        
        %If we want to plot synchronous with animation, divide into
        %subplots and do that
        if doOtherPlot
            subplot(1,2,2);
            curvs = flipud(states(4:end,i));
            plot(curvs);
            axis([0,nLinks+1,-maxCurv,maxCurv]);
            title('Curvature Over Body');
            xlabel('Links');
            xticks([1,nLinks]);
            xticklabels({'Head','Tail'});
            ylabel('Curvature (rad)');
            subplot(1,2,1);
            
        end
            
        axis equal;
        hold on;

        title('Salp Animation');
        
        %Plot all the background lines for judging swimmer motion
        for j = 1:numel(xTicks)
            plot([xTicks(j),xTicks(j)],[plotMinY,plotMaxY],'Color',backgroundTicksColor);
        end
        for j = 1:numel(yTicks)
            plot([plotMinX,plotMaxX],[yTicks(j),yTicks(j)],'Color',backgroundTicksColor);
        end
    
        %Store swimmer data at this time
        thisState = states(:,i);
        theseForces = forces(:,i);
    
        %Get position and orientation of tail link
        linkPosition = thisState(1:2);
        linkTheta = thisState(3);

        %Vector for holding summed link positions - used for finding
        %swimmer c.g. when we choose the axis limits for each frame
        positions = zeros(2,1);
    
        %For every link in the swimmer
        for link = 1:nLinks

            %Skip first three values in state since they represent world
            %frame positions. This is index of joint shape value
            linkIndex = link+3;

            %Get the force magnitude acting on this link
            forceIndices = 3*(link-1)+[1:3];
            thisForceVector = theseForces(forceIndices);
            forceMagnitude = norm(thisForceVector(1:2));
            normalizedForceMagnitude = forceMagnitude/salp.thrustAmp;

            %Scale thruster animation by normalized force magnitude
            thisThrusterOut = outForce;
            thisThrusterOut(1,:) = thisThrusterOut(1,:)*normalizedForceMagnitude;

            %Make the inside thrust color half the length of the outside color
            thisThrusterIn = outForce;
            thisThrusterIn(1,:) = .5*thisThrusterIn(1,:)*normalizedForceMagnitude;
            thisThrusterIn(2,:) = .5*thisThrusterIn(2,:);

            %Choose whether this thruster is left-sided or right-sided
            %based on link index
            if mod(link,2)
                thrusterTheta = theta_left;
            else
                thrusterTheta = theta_right;
            end
    
            %Move ellipse, thruster, and thrust to the correct
            %position/orientation
            linkPoints = shiftPoints(ellipsePoints,linkTheta,linkPosition);
            linkThrusterPoints = shiftPoints(thrusterPoints,linkTheta+thrusterTheta,linkPosition);
            linkThrustOutPoints = shiftPoints(thisThrusterOut,linkTheta+thrusterTheta,linkPosition);
            linkThrustInPoints = shiftPoints(thisThrusterIn,linkTheta+thrusterTheta,linkPosition);

            %Add ellipse, thruster, and thrust to the plot
            fill(linkPoints(1,:),linkPoints(2,:),linkColor);
            fill(linkThrusterPoints(1,:),linkThrusterPoints(2,:),thrusterColor);
            fill(linkThrustOutPoints(1,:),linkThrustOutPoints(2,:),outForceColor);
            fill(linkThrustInPoints(1,:),linkThrustInPoints(2,:),inForceColor);

            %Keep track of summed link positions for c.g. calculations
            positions = positions + linkPosition;

            %If the link isn't the head
            if link < nLinks
                %Walk forward half the length on the current length
                dX = getRotMatrix(linkTheta)*[a;0];
                %Rotate by the joint angle
                linkTheta = linkTheta + thisState(linkIndex);
                %Walk forward half the length on the next link
                linkPosition = linkPosition + dX + getRotMatrix(linkTheta)*[a;0];
            end
        end

        %Calculate swimmer c.g. for choosing axis limits
        positions = positions/nLinks;
        %Provide buffer in all directions around c.g.
        minX = positions(1)-nLinks*L*2/3;
        maxX = positions(1)+nLinks*L*2/3;
        minY = positions(2)-nLinks*L*2/3;
        maxY = positions(2)+nLinks*L*2/3;

        %Fit axis so we're always zoomed on the swimmer
        axis([minX,maxX,minY,maxY]);
        %Set figure position on screen
        f10.Position = [100,120,1800,800];
        %Make sure everything has finished drawing
        drawnow;
    
        %If this is the first frame
        if i == 1
            try
                %Start the gif
                gif('SalpSwimming.gif','DelayTime',dt);
            catch
                error('Animation requires "gif" package to be installed');
            end 
        else
            %Add frame to gif
            gif;
        end
    
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