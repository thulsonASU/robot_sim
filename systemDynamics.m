function dydt = systemDynamics(t, y)
    % Extract joint angles and their derivatives from y
    t1 = y(1);
    t2 = y(2);
    t3 = y(3);
    t1_d = y(4);
    t2_d = y(5);
    t3_d = y(6);
    t1_dd = y(7);
    t2_dd = y(8);
    t3_dd = y(9);
    
    % Calculate torques Tao_1, Tao_2, Tao_3
    Tao_1 = t3_dd*((3*cos(t2 + t3))/2 + cos(t3) + 9/4) - t1_d*((t3_d*(3*sin(t2 + t3) + 2*sin(t3)))/2 + 3*t2_d*(sin(t2 + t3)/2 + 5*sin(t2))) + t2_dd*((3*cos(t2 + t3))/2 + 21*cos(t2) + 2*cos(t3) + 89/4) + t1_dd*(3*cos(t2 + t3) + 30*cos(t2) + 2*cos(t3) + 107/2) - t2_d*(3*t1_d*(sin(t2 + t3)/2 + 5*sin(t2)) + (t3_d*(3*sin(t2 + t3) + 2*sin(t3)))/2 + 3*t2_d*(sin(t2 + t3)/2 + 7*sin(t2))) - (t3_d*(3*sin(t2 + t3) + 2*sin(t3))*(t1_d + t2_d + t3_d))/2;
    Tao_2 = t3_dd*(cos(t3) + 9/4) + t1_d*(t1_d*((3*sin(t2 + t3))/2 + 15*sin(t2)) - t3_d*sin(t3)) + t2_dd*(12*cos(t2) + 2*cos(t3) + 89/4) - t2_d*(6*t2_d*sin(t2) + t3_d*sin(t3)) + t1_dd*((3*cos(t2 + t3))/2 + 21*cos(t2) + 2*cos(t3) + 89/4) - t3_d*sin(t3)*(t1_d + t2_d + t3_d);
    Tao_3 = (9*t3_dd)/4 + t2_dd*(cos(t3) + 9/4) + t1_dd*((3*cos(t2 + t3))/2 + cos(t3) + 9/4) + t1_d*(t1_d*((3*sin(t2 + t3))/2 + sin(t3)) + t2_d*sin(t3)) + t2_d*sin(t3)*(t1_d + t2_d);
    
    % Calculate the derivatives of the joint angles
    dydt = [t1_d; t2_d; t3_d; t1_dd; t2_dd; t3_dd; Tao_1; Tao_2; Tao_3];

end