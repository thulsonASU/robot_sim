function dqdt = systemDynamics(t, q, Tao, NumberOfJoints)
    % Initialize dqdt vector
    dqdt = zeros(NumberOfJoints, 1);
    
    % Compute the dynamics for each joint
    for i = 1:NumberOfJoints
        % Assume a simple model where the joint angle rate of change (dqdt) is directly proportional to the torque
        % This would be the case if the moment of inertia and the angular acceleration were constant and equal to 1
        dqdt(i) = Tao(i);
    end
end