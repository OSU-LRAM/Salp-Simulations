%generates vector of thruster forces for each link as a function of time -
%sunusoidal thrust of magnitude 0 to thrustAmp

%out:
%forcingFunctions - vector of thrust functions for each link

%in:
%nLinks - number of links in swimming system
%thrustAmp - maximum amplitude of thrust (thrust bound is [0,thrustAmp]
%frequency - frequency of thrust function
%thrustAngle - angle of thrust with respect to tail of link
%phaseDelays - vector of phase offsets for each link with respect to the
%distal link (head link should have phase offset of zero)

function forcingFunctions = generateForceVectorFunctions(nLinks,thrustAmp,frequency,thrustAngle,phaseDelays)

    %Matrix that rotates +X vector to be in direction of thrust for a jet
    %pointing out the 'left' side of the swimmer
    R_left = [cos(pi-thrustAngle),-sin(pi-thrustAngle),0;sin(pi-thrustAngle),cos(pi-thrustAngle),0;0,0,0];
    %Matrix that rotates +X vector to right-side jet
    R_right = [cos(pi-thrustAngle),sin(pi-thrustAngle),0;-sin(pi-thrustAngle),cos(pi-thrustAngle),0;0,0,0];

    %Cell vector to hold individual time functions of thrust for each link
    functionCells = cell(1);

    sameSide = 0;
    %For every link
    for i = 1:nLinks
        %Determine if jet will point left or right, alternate them
        if mod(i,2)
            R = R_left;
        else
            if sameSide
                R = R_left;
            else
                R = R_right;
            end
        end
        %Make force vector in link frame of sinusoidal jet action
        %Jet magnitude goes from 0 to thrustAmp at desired thrustAngle
        theseForces = @(t) -R*[thrustAmp/2*sin(2*pi*frequency*t + phaseDelays(i))+thrustAmp/2;0;0];
        %Save link forcing function into cell array
        functionCells{i} = theseForces;
    end

    %Make anonymous function that evaluates all force functions at the same
    %time and returns them as a giant column vector
    forcingFunctions = @(t) combineAllForces(t,functionCells);

end

%Evaluates force function for each link at a given time t and returns them
%all as a giant column vector F
function F = combineAllForces(t,functionCells)

    %Initialize F as a zero-vector
    F = zeros(3*numel(functionCells),1);

    index = 1;
    %For every link (and associated forcing function)
    for i = 1:numel(functionCells)
        %Add evaluated force function to complete force vector
        F(index:index+2) = functionCells{i}(t);
        index = index + 3;
    end

end
