% global a k_r1 k_r2 pi_m pi_l


impedanceGUI = uifigure('Name', 'Impedance Control Simulink Model', 'NumberTitle', 'off');

% load manipulator dynamic parameters without load mass
%PARAM  Augmented link parameters.

% L. Villani, G. Oriolo, B. Siciliano
% February 2009

% link parameters
  % lengths
    a  = [1;1];

  % masses
    m_l1  = 50;
    m_l2  = 50;

  % distances of centers of mass from joint axes
    l_1  =  0.5;
    l_2  =  0.5;

  % inertia moments relative to the centers of mass
    I_l1 = 10;
    I_l2 = I_l1;

% motor parameters
  % masses
    m_m1 = 5;
    m_m2 = 5;

  % inertia moments
    I_m1 = 0.01;
    I_m2 = 0.01;

  % transmission ratios
    k_r1 = 100;
    k_r2 = 100;

% augmented link parameters
  % masses
    m_1 = m_l1 + m_m2;
    m_2 = m_l2;

  % first moments of inertia
    m1_lC1 = m_l1*(l_1 - a(1));
    m2_lC2 = m_l2*(l_2 - a(2));

  % inertia moments relative to origins of link frames
    I_1 = I_l1 + m_l1*(l_1 - a(1))^2 + I_m2;
    I_2 = I_l2 + m_l2*(l_2 - a(2))^2;

% minimum set of dynamic parameters
  pi_m(1) = a(1)*m_1 + m1_lC1 + a(1)*m_2;
  pi_m(2) = a(1)*m1_lC1 + I_1 + k_r1^2*I_m1;
  pi_m(3) = a(2)*m_2 + m2_lC2;
  pi_m(4) = a(2)*m2_lC2 + I_2;
  pi_m(5) = I_m2;

  pi_l = pi_m;

% gravity acceleration
  g = 9.81;

% friction matrix
  K_r = diag([k_r1 k_r2]);
  F_v = K_r*diag([0.01 0.01])*K_r;

% sample time of controller
  Tc = 0.001;

% impedance controller gains
  K_d = 1600*diag([1 1]);
  K_p = 5000*diag([1 1]);
  M_d = 100*diag([1 1]);
  iM_d = inv(M_d);

% constraint frame matrix
  R_c = [cos(pi/4)  -sin(pi/4);
        sin(pi/4)  cos(pi/4)];

% stiffness matrix in base frame
  K = R_c*[0 0;0 5e3]*R_c';

% position of undeformed plane
  o_r = [1;0];

% initial position
  p_i = [1+0.2*cos(pi/4); 0];

% final position
  p_f = [1.2+0.2*cos(pi/4);0.2];

% initial joint configuration
  q_i = inv_k2u(a,p_i);

% duration of simulation
  t_d = 2.5;

% trajectory in base frame
  tra_5;

% sample time for plots
  Ts = Tc;