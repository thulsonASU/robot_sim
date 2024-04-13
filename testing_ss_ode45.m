% Define the time span
tspan = [0, 10]; % Time span from 0 to 10 seconds

% Define initial conditions
tao0 = 0; % Initial torque

% Solve the differential equation for motor torque using ode45
[t, tao] = ode45(@dynamics_test, tspan, tao0);

% Plot the motor torque vs. time
plot(t, tao, 'b');
xlabel('Time');
ylabel('Motor Torque (\tau)');
title('Motor Torque vs Time');
