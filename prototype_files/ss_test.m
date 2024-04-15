% Define the symbolic variables and equations
syms q1 q2 q3 q1_dot q2_dot q3_dot u1 u2 u3 t1 t2 t3 t1_d t2_d t3_d t1_dd t2_dd t3_dd w_a1 w_a2 w_a3 Tao_1 Tao_2 Tao_3

% % Calculate torques Tao_1, Tao_2, Tao_3
% Tao_1 = t3_dd*((3*cos(t2 + t3))/2 + cos(t3) + 9/4) - t1_d*((t3_d*(3*sin(t2 + t3) + 2*sin(t3)))/2 + 3*t2_d*(sin(t2 + t3)/2 + 5*sin(t2))) + t2_dd*((3*cos(t2 + t3))/2 + 21*cos(t2) + 2*cos(t3) + 89/4) + t1_dd*(3*cos(t2 + t3) + 30*cos(t2) + 2*cos(t3) + 107/2) - t2_d*(3*t1_d*(sin(t2 + t3)/2 + 5*sin(t2)) + (t3_d*(3*sin(t2 + t3) + 2*sin(t3)))/2 + 3*t2_d*(sin(t2 + t3)/2 + 7*sin(t2))) - (t3_d*(3*sin(t2 + t3) + 2*sin(t3))*(t1_d + t2_d + t3_d))/2;
% Tao_2 = t3_dd*(cos(t3) + 9/4) + t1_d*(t1_d*((3*sin(t2 + t3))/2 + 15*sin(t2)) - t3_d*sin(t3)) + t2_dd*(12*cos(t2) + 2*cos(t3) + 89/4) - t2_d*(6*t2_d*sin(t2) + t3_d*sin(t3)) + t1_dd*((3*cos(t2 + t3))/2 + 21*cos(t2) + 2*cos(t3) + 89/4) - t3_d*sin(t3)*(t1_d + t2_d + t3_d);
% Tao_3 = (9*t3_dd)/4 + t2_dd*(cos(t3) + 9/4) + t1_dd*((3*cos(t2 + t3))/2 + cos(t3) + 9/4) + t1_d*(t1_d*((3*sin(t2 + t3))/2 + sin(t3)) + t2_d*sin(t3)) + t2_d*sin(t3)*(t1_d + t2_d);

% convert Tao to a system of equations
eq1 = u1 == t3_dd*((3*cos(t2 + t3))/2 + cos(t3) + 9/4) - t1_d*((t3_d*(3*sin(t2 + t3) + 2*sin(t3)))/2 + 3*t2_d*(sin(t2 + t3)/2 + 5*sin(t2))) + t2_dd*((3*cos(t2 + t3))/2 + 21*cos(t2) + 2*cos(t3) + 89/4) + t1_dd*(3*cos(t2 + t3) + 30*cos(t2) + 2*cos(t3) + 107/2) - t2_d*(3*t1_d*(sin(t2 + t3)/2 + 5*sin(t2)) + (t3_d*(3*sin(t2 + t3) + 2*sin(t3)))/2 + 3*t2_d*(sin(t2 + t3)/2 + 7*sin(t2))) - (t3_d*(3*sin(t2 + t3) + 2*sin(t3))*(t1_d + t2_d + t3_d))/2;
eq2 = u2 == t3_dd*(cos(t3) + 9/4) + t1_d*(t1_d*((3*sin(t2 + t3))/2 + 15*sin(t2)) - t3_d*sin(t3)) + t2_dd*(12*cos(t2) + 2*cos(t3) + 89/4) - t2_d*(6*t2_d*sin(t2) + t3_d*sin(t3)) + t1_dd*((3*cos(t2 + t3))/2 + 21*cos(t2) + 2*cos(t3) + 89/4) - t3_d*sin(t3)*(t1_d + t2_d + t3_d);
eq3 = u3 == (9*t3_dd)/4 + t2_dd*(cos(t3) + 9/4) + t1_dd*((3*cos(t2 + t3))/2 + cos(t3) + 9/4) + t1_d*(t1_d*((3*sin(t2 + t3))/2 + sin(t3)) + t2_d*sin(t3)) + t2_d*sin(t3)*(t1_d + t2_d);

% Define the state variables
x1 = t1;
x2 = t2;
x3 = t3;
x4 = t1_d;
x5 = t2_d;
x6 = t3_d;
% Define the input variables
u = [u1; u2; u3];
% Define the state equations
dx1 = x4;
dx2 = x5;
dx3 = x6;
dx4 = simplify(solve(eq1, t1_dd));
dx5 = simplify(solve(eq2, t2_dd));
dx6 = simplify(solve(eq3, t3_dd));

dx = [dx1; dx2; dx3; dx4; dx5; dx6];                                                                                                                                                                                                                             (4*u3)/9 - (4*t2_dd*(cos(t3) + 9/4))/9 - (4*t1_dd*((3*cos(t2 + t3))/2 + cos(t3) + 9/4))/9 - (4*t1_d*(t1_d*((3*sin(t2 + t3))/2 + sin(t3)) + t2_d*sin(t3)))/9 - (4*t2_d*sin(t3)*(t1_d + t2_d))/9

% Define the state vector
x = [x1; x2; x3; x4; x5; x6];
% Define the state matrix (A)
A = jacobian(dx, x)
% Define the input matrix (B)
B = jacobian([dx1; dx2; dx3; dx4; dx5; dx6], u)

% Define the initial conditions
x0 = [0; 0; 0; 0; 0; 0]; % Example initial conditions
% Define the time span for the simulation
tspan = [0 10]; % Simulate from t=0 to t=10
% Define the input function (if the input is time-varying)
u = @(t) [0; 0; 0]; % Example input function

q = [t1; t2; t3; t1_d; t2_d; t3_d];

disp(class(x))
% Simulate the system
[t, x] = ode45(@(t, x) systemDynamics_ss_test(t, x, u(t), A, B), tspan, x0);

% Plot the results
figure;
plot(t, x);
xlabel('Time (s)');
ylabel('State Variables');
legend('q1', 'q2', 'q3', 'q1_dot', 'q2_dot', 'q3_dot');
title('State Variables vs. Time');