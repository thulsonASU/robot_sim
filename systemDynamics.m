function dydt = systemDynamics(t, y, eqs, q_all, NumberOfJoints)

    for i = 1:NumberOfJoints
        % Extract the joint angles, velocities, and accelerations from the input vector
        Y(i) = y(i);
        Y(i + NumberOfJoints) = y(i + NumberOfJoints);
        Y(i + 2*NumberOfJoints) = y(i + 2*NumberOfJoints);
    end

    % Old Example Equations
    % % Calculate torques Tao_1, Tao_2, Tao_3
    % Tao_1 = t3_dd*((3*cos(t2 + t3))/2 + cos(t3) + 9/4) - t1_d*((t3_d*(3*sin(t2 + t3) + 2*sin(t3)))/2 + 3*t2_d*(sin(t2 + t3)/2 + 5*sin(t2))) + t2_dd*((3*cos(t2 + t3))/2 + 21*cos(t2) + 2*cos(t3) + 89/4) + t1_dd*(3*cos(t2 + t3) + 30*cos(t2) + 2*cos(t3) + 107/2) - t2_d*(3*t1_d*(sin(t2 + t3)/2 + 5*sin(t2)) + (t3_d*(3*sin(t2 + t3) + 2*sin(t3)))/2 + 3*t2_d*(sin(t2 + t3)/2 + 7*sin(t2))) - (t3_d*(3*sin(t2 + t3) + 2*sin(t3))*(t1_d + t2_d + t3_d))/2;
    % Tao_2 = t3_dd*(cos(t3) + 9/4) + t1_d*(t1_d*((3*sin(t2 + t3))/2 + 15*sin(t2)) - t3_d*sin(t3)) + t2_dd*(12*cos(t2) + 2*cos(t3) + 89/4) - t2_d*(6*t2_d*sin(t2) + t3_d*sin(t3)) + t1_dd*((3*cos(t2 + t3))/2 + 21*cos(t2) + 2*cos(t3) + 89/4) - t3_d*sin(t3)*(t1_d + t2_d + t3_d);
    % Tao_3 = (9*t3_dd)/4 + t2_dd*(cos(t3) + 9/4) + t1_dd*((3*cos(t2 + t3))/2 + cos(t3) + 9/4) + t1_d*(t1_d*((3*sin(t2 + t3))/2 + sin(t3)) + t2_d*sin(t3)) + t2_d*sin(t3)*(t1_d + t2_d);
    
    % tao is calculated with the rhs from the subbed equations
    for i = 1:NumberOfJoints
        Tao(i) = subs(eqs(i), q_all, Y);
    end
    
    % Calculate the derivatives of the joint angles
    dydt = [transpose(Y(NumberOfJoints+1:end)); transpose(Tao)];
    
    % convert dydt to double
    dydt = double(dydt);
end