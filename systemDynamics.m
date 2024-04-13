function dXdt = systemDynamics(t, X, tao)
    % Extract state variables from X
    t1 = X(1);
    t2 = X(2);
    t3 = X(3);
    t1_d = X(4);
    t2_d = X(5);
    t3_d = X(6);

    % Define the symbolic variables
    syms t1_dd_sym t2_dd_sym t3_dd_sym
    
    % Define your equations for joint torques (Tau) here
    % Replace the example equations with your actual equations
    eq1 = t3_dd_sym*((3*cos(t2 + t3))/2 + cos(t3) + 9/4) - t1_d*((t3_dd_sym*(3*sin(t2 + t3) + 2*sin(t3)))/2 + 3*t2_d*(sin(t2 + t3)/2 + 5*sin(t2))) + t2_dd_sym*((3*cos(t2 + t3))/2 + 21*cos(t2) + 2*cos(t3) + 89/4) + t1_dd_sym*(3*cos(t2 + t3) + 30*cos(t2) + 2*cos(t3) + 107/2) - t2_d*(3*t1_d*(sin(t2 + t3)/2 + 5*sin(t2)) + (t3_dd_sym*(3*sin(t2 + t3) + 2*sin(t3)))/2 + 3*t2_d*(sin(t2 + t3)/2 + 7*sin(t2))) - (t3_dd_sym*(3*sin(t2 + t3) + 2*sin(t3))*(t1_d + t2_d + t3_dd_sym))/2;
    eq2 = t3_dd_sym*(cos(t3) + 9/4) + t1_d*(t1_d*((3*sin(t2 + t3))/2 + 15*sin(t2)) - t3_dd_sym*sin(t3)) + t2_dd_sym*(12*cos(t2) + 2*cos(t3) + 89/4) - t2_d*(6*t2_d*sin(t2) + t3_dd_sym*sin(t3)) + t1_dd_sym*((3*cos(t2 + t3))/2 + 21*cos(t2) + 2*cos(t3) + 89/4) - t3_dd_sym*sin(t3)*(t1_d + t2_d + t3_dd_sym);
    eq3 = (9*t3_dd_sym)/4 + t2_dd_sym*(cos(t3) + 9/4) + t1_dd_sym*((3*cos(t2 + t3))/2 + cos(t3) + 9/4) + t1_d*(t1_d*((3*sin(t2 + t3))/2 + sin(t3)) + t2_d*sin(t3)) + t2_d*sin(t3)*(t1_d + t2_d);
    
    % Define the equations in terms of the symbolic variables
    eq1_sym = eq1 == t1_dd_sym;
    eq2_sym = eq2 == t2_dd_sym;
    eq3_sym = eq3 == t3_dd_sym;
    
    % Solve the system of equations
    sol = solve([eq1_sym, eq2_sym, eq3_sym], [t1_dd_sym, t2_dd_sym, t3_dd_sym]);
    
    % Extract the solutions
    t1_dd = double(sol.t1_dd_sym);
    t2_dd = double(sol.t2_dd_sym);
    t3_dd = double(sol.t3_dd_sym);
    
    % Pack the derivatives into a column vector
    dXdt = [t1_d; t2_d; t3_d; t1_dd; t2_dd; t3_dd];
end