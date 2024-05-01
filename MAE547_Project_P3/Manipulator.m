function [sys,x0,str,ts] = Manipulator(t,x,u,flag)

switch flag
  case 0
    [sys,x0,str,ts]=mdlInitializeSizes();
  case 1
    sys=mdlDerivatives(t,x,u);
  case 3
    sys=mdlOutputs(t,x,u);
  case { 2, 4, 9 }
    sys = [];
  otherwise
    DAStudio.error('Simulink:blocks:unhandledFlag', num2str(flag));
end

function [sys,x0,str,ts]=mdlInitializeSizes()
    dof = 3;
    sizes = simsizes;
    sizes.NumContStates  = 2*dof;
    sizes.NumDiscStates  = 0;
    sizes.NumOutputs     = 2*dof;
    sizes.NumInputs      = dof; %3*dof
    sizes.DirFeedthrough = 0;
    sizes.NumSampleTimes = 1;

    sys = simsizes(sizes);
    x0  = zeros(6,1);
    str = [];
    ts  = [0 0];


function sys=mdlDerivatives(t, x, u)
    q = x(1:3);
    d_q = x(4:6);
    tau = u(1:3);
    % he = u(4:9);
    r1 = 0.02;
    b2 = 0.03;
    b3 = 0.01;
    m1 = 1;
    m2 = 1;
    m3 = 1;
    l1 = 0.4;
    l2 = 0.46;
    l3 = 0.15;
    g = 9.81;

    B = double(B_Lagrangian(q,b2,b3,m3,l1,m1,m2,l2,l3,r1));
    C = double(C_Lagrangian(q, d_q,m3,l1));
    G = double(G_Lagrangian(q,m3,l1,m1,m2));
    % J = double(Jacobian(q, param));
    
    dd_q = inv(B)*(tau - C*d_q - G);
    % dd_q = inv(B)*(tau - C*d_q - G - J'*he);

    sys = [d_q; dd_q];

function sys=mdlOutputs(t, x, u)
    sys = x;