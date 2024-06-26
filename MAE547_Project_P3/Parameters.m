% URDF Parameters
param.m1 = 1.3565;
param.m2 = 1.1178;
param.m3 = 0.0405;
param.pi = 3.14;

% From the urdf provided
param.l0 = 0.4;
param.l3 = 0.3/2;
param.l1 = 0.4;
param.l2 = 0.46;
param.r1 = 0.02;
param.b2 = 0.03;
param.b3 = 0.01;
param.g = 9.81;

l1 = param.l1;
l2 = param.l2;
l3 = param.l3;

t1 = 0.1;
d2 = -0.1;
d3 = 0.1;

xd = Kinematics(q);

% PD Controller parameters
Kp = diag([50 80 100 1 1 1]);
Kd = diag([2 10 30 1 1 1]);

K = 200*diag([1 1 1 1 1 1]);
param.Kp = Kp;
param.K = K;

% Environment only one component for a simplified force control
param.wall = 0.6;
param.axis = 1; % 1-x 2-y 3-z wrt of frame 0