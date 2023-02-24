%Calculate body and joint velocity given a swimmer shape and set of thrust
%vectors acting on each link. Requires access to sysplotter functions

%out:
%speeds - column vector of body velocities and shape velocities
%in:
%shape - shape configuration of joints
%Ft - all thrust vectors e.g. [F1;F2;F3] where F1 = [Fx1;Fy1;Ftheta1]
%as a function of time
%L - length of one link
%d - drag metric on each link s.t. Fi = d*vi for v in link frame
%K - matrix of joint stiffnesses s.t. [Fb;Fr] = K*r
function speeds = thrustToSpeeds(t,shape,Ft,L,d,K)

    %Make shape into column vector, establish number of links in the system
    shape = shape(:);
    numLinks = numel(shape)+1;

    %Make geometry structure used to calculate body jacobians
    geom = struct();
    geom.linklengths = L*ones(1,numLinks);
    %Make center of tail base frame for ease of animation
    geom.baseframe = 'tail';

    %Calculate body jacobian for each link
    [~,~,J_full] = N_link_chain(geom,shape);

    %Calculate distribution matrix and drag metric
    B = getDistribution_discrete(J_full);
    D = getDragMetric(d,J_full);

    %Collect forces acting on body frame and joints from thrusters and
    %springs
    Forces = B*Ft(t) + K*shape;
    %Calculate body frame and joint velocities that would result from these
    %forces
    speeds = inv(D)*Forces;

end

%Calculate drag metric from individual link drag matrices
function D = getDragMetric(d,J_full)

    %Calculate size of output
    %n = 3 + number of joints
    n = size(J_full{1},2);
    %initialize empty drag metric
    D = zeros(n);
    %For every link
    for i = 1:numel(J_full)
        %Pull back forces on link into chosen body frame using jacobians
        D = D + J_full{i}'*d*J_full{i};
    end

end