%generates vector of thruster forces for each link as a function of time,
%where each link has an uncertain delay after firing

%out:
%forcingFunctions - vector of thrust functions for each link

%in:
%nLinks - number of links in swimming system
%thrustAmp - maximum amplitude of thrust (thrust bound is [0,thrustAmp]
%frequency - frequency of thrust function
%thrustAngle - angle of thrust with respect to tail of link
%numPeriods - number of 'Thrust' periods to simulate
%percentThrusting - percentage of time a zooid is thrusting on average,
%scalar from 0 (never firing) to 1 (always firing)

function forcingFunctions = generateDispersedForceVectorFunctions(nLinks,thrustAmp,frequency,thrustAngle,numPeriods,percentThrusting)

    %Cell vector to hold individual time activations of thrust for each link
    functionCells = cell(1);

    sameSide = 0;
    
    %Calculate average time between thrust finishing and new ignition
    T = 1/frequency;
    if percentThrusting == 0
        deltaT = 1e6;
    else
        deltaT = T/percentThrusting - T;
    end
    
    %For every link, generate a vector of thrust ignition times based on
    %the average thruster uptime
    for i = 1:nLinks
        ignitionVector = [-rand()*(T+deltaT)];
        while ignitionVector(end) < T*numPeriods;
            ignitionVector(end+1) = ignitionVector(end) + T + mean(rand(1,3))*2*deltaT;
        end
        functionCells{i} = ignitionVector;
    end

    %Make anonymous function that evaluates all force functions at the same
    %time and returns them as a giant column vector
    forcingFunctions = @(t) combineAllForces(t,functionCells,sameSide,thrustAmp,frequency,thrustAngle);

end

%Evaluates force function for each link at a given time t and returns them
%all as a giant column vector F
function F = combineAllForces(t,functionCells,sameSide,thrustAmp,frequency,thrustAngle)

    %Matrix that rotates +X vector to be in direction of thrust for a jet
    %pointing out the 'left' side of the swimmer
    R_left = [cos(pi-thrustAngle),-sin(pi-thrustAngle),0;sin(pi-thrustAngle),cos(pi-thrustAngle),0;0,0,0];
    %Matrix that rotates +X vector to right-side jet
    R_right = [cos(pi-thrustAngle),sin(pi-thrustAngle),0;-sin(pi-thrustAngle),cos(pi-thrustAngle),0;0,0,0];
    
    T = 1/frequency;
    
    %Initialize F as a zero-vector
    F = zeros(3*numel(functionCells),1);

    index = 1;
    %For every link (and associated forcing function)
    for i = 1:numel(functionCells)
        
        %Determine if jets will point left or right, or if alternating
        if mod(i,2)
            R = R_left;
        else
            if sameSide
                R = R_left;
            else
                R = R_right;
            end
        end
        
        %Calculate thrust amplitude at this time
        ignitions = functionCells{i};
        relevantIgnition = ignitions(find(ignitions > t,1)-1);
        thrustT = t-relevantIgnition;
        if thrustT > T
            thrust = zeros(3,1);
        else
            thrust = -R*[-thrustAmp/2*cos(2*pi*frequency*thrustT)+thrustAmp/2;0;0];
        end
        
        %Add evaluated force function to complete force vector
        F(index:index+2) = thrust;
        index = index + 3;
    end

end
