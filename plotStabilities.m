%Generates gif of joint stability vector field and joint flow from 4
%initial conditions over varying stiffnesses

%Remember to initialize sysplotter and have gif Add-On installed
ks = [0:.001/3:.075];
dt = 1/30;

%Declare salp characteristic values
salp = struct();
%Length of individual zooid units
salp.L = .25;
%Ratio of Y drag to X drag of a link - Y drag probably higher
salp.dragRatio = 3;

%Initialize figure and figure size
fig2 = figure(2);
fig2.Position = [97 145 1648 518];
set(gcf,'color','w');

%For each stiffness we want to visualize
for stiffness_index = 1:numel(ks)

    %Set salp stiffness
    salp.jointStiffness = ks(stiffness_index);
    
    %Set initial joint conditions using intermediate polar coordinates
    mag = pi/2;
    thetas = [pi/4,3*pi/4,5*pi/4,7*pi/4];
    startPoints = zeros(numel(thetas),2);
    for i = 1:numel(thetas)
        startPoints(i,:) = [mag*cos(thetas(i)),mag*sin(thetas(i))];
    end
    
    %Grab the plot data
    results = makeStabilityDiagram(salp,startPoints);
    
    clf;
    tcl = tiledlayout(1,3);
    %Plot each stability vector field and the four joint flows
    for i = 1:3
        nexttile;
        quiver(results{i}.R1,results{i}.R2,results{i}.dR1,results{i}.dR2);
        hold on;
        for j = 1:size(results{i}.pathx,1)
            plot(results{i}.pathx(j,:),results{i}.pathy(j,:),'LineWidth',2);
        end
        axis square
        xticks([-pi,0,pi]);
        yticks([-pi,0,pi]);
        xticklabels({'-\pi','0','\pi'});
        yticklabels({'-\pi','0','\pi'});
        title(['Link ',num2str(i),' Firing']);
        xlabel('Tail Joint');
        ylabel('Head Joint');
    end
    
    title(tcl,{['Stiffness: ',num2str(round(salp.jointStiffness,3))],''});
    drawnow;

    if stiffness_index == 1
        gif('Stability vs Stiffness.gif','DelayTime',dt);
    else
        gif;
    end
end