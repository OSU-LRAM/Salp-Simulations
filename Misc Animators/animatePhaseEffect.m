%Animates effect of phase delay between links in salp chain

gifTime = 10;
FPS = 30;
phases = linspace(0,2*pi,gifTime*FPS);

numLinks = 20;

lw = .1;

L = .25;
a = L/2;
AR = .5;
ellipseThetas = linspace(0,2*pi,100);
ellipsePoints = [a*cos(ellipseThetas);a*AR*sin(ellipseThetas)];

clockR = .75*numLinks*L/2;
clockX = numLinks*L/2;
clockY = numLinks*L/2 + .25;
clockPoints = [clockR*cos(ellipseThetas) + clockX;clockR*sin(ellipseThetas) + clockY];
clockWidth = 3;

handWidth = 3;
handR = .75*clockR;

fig10 = figure(10);
set(fig10,'Color',[1,1,1]);
for i = 1:numel(phases)
    clf;
    hold on;

    thisPhase = phases(i);
    phaseDelay = -phases(i);
    linkX = a;
    thisEllipse = ellipsePoints;

    plot(clockPoints(1,:),clockPoints(2,:),'k','LineWidth',clockWidth);

    handPoints = [clockX,clockX + handR*cos(thisPhase);clockY,clockY + handR*sin(thisPhase)];
    plot(handPoints(1,:),handPoints(2,:),'k','LineWidth',handWidth);
    plot([clockX,clockX+clockR],[clockY,clockY],'k--','LineWidth',.5)

    text(clockX,5,'Phase Delay','FontSize',12,'FontWeight','bold','HorizontalAlignment','center')

    axis equal;
    axis([-.5,numLinks*L+.5,-.5,5.5]);
    drawnow;

    for j = 1:numLinks
        thisPhase = thisPhase + phaseDelay;
        thisEllipse = thisEllipse + [L;0];
        thisMag = .5*cos(thisPhase) + .5;

        color = [1,1-thisMag,1-thisMag];
        fill(thisEllipse(1,:),thisEllipse(2,:),color);
        plot(thisEllipse(1,:),thisEllipse(2,:),'k','LineWidth',lw);
    end


    
    if i==1
        gif('PhaseIllustration.gif','DelayTime',1/FPS);
    else
        gif;
    end
end
