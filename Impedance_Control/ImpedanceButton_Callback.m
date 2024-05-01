function ImpedanceButton_Callback()
    % Create a new figure window
    ImpedanceFig = uifigure('Name', 'Impedance Simulation', 'NumberTitle', 'off');

    closeButton = uibutton(ImpedanceFig, 'Position', [10, 10, 100, 22], 'Text', 'Close', 'ButtonPushedFcn', @(btn,event) close(ImpedanceFig));
    closeButton.BackgroundColor = [0.7 0.2 0.2]; % Red color

    % Add UI controls to the new figure window as needed
    % Add labels for instructions
    uilabel(ImpedanceFig, 'Text', 'Each comma seperated value corresponds to each joint in series. ', 'Position', [20, 395, 400, 22]);
    uilabel(ImpedanceFig, 'Text', 'Link parameters', 'Position', [40, 360, 500, 22]);
    
    % Add labels for each field
    uilabel(ImpedanceFig, 'Text', 'Link Lengths:', 'Position', [10, 340, 300, 20]);
    uilabel(ImpedanceFig, 'Text', 'Link Masses:', 'Position', [10, 320, 300, 20]);
    uilabel(ImpedanceFig, 'Text', 'Distance of CM from Joint Axis:', 'Position', [10, 300, 300, 20]);
    uilabel(ImpedanceFig, 'Text', 'Inertia Moments:', 'Position', [10, 280, 300, 20]);

    %Motor parameters
    uilabel(ImpedanceFig, 'Text', 'Motor Masses:', 'Position', [10, 260, 300, 20]);
    uilabel(ImpedanceFig, 'Text', 'Motor Inertia Moments:', 'Position', [10, 240, 300, 20]);
    uilabel(ImpedanceFig, 'Text', 'Transmission Ratios:', 'Position', [10, 220, 300, 20]);
    uilabel(ImpedanceFig, 'Text', 'Gravity:', 'Position', [10, 200, 300, 20]);

    %Impedance control gains
    uilabel(ImpedanceFig, 'Text', 'KD:', 'Position', [10, 180, 300, 20]);
    uilabel(ImpedanceFig, 'Text', 'KP:', 'Position', [10, 160, 300, 20]);
    uilabel(ImpedanceFig, 'Text', 'MD:', 'Position', [10, 140, 300, 20]);
    uilabel(ImpedanceFig, 'Text', 'Duration of Simulation:', 'Position', [10, 120, 300, 20]);
    uilabel(ImpedanceFig, 'Text', 'Sample Time of Controller:', 'Position', [10, 100, 300, 20]);

    % Make a ui field to fill in for Kr MassLI MassMI IntertiaLI IntertiaMI center of mass lengths
    a_Field = uieditfield(ImpedanceFig, 'text', 'Position', [200, 340, 100, 20], 'Value', '1,1');
    ml_Field= uieditfield(ImpedanceFig, 'text', 'Position', [200, 320, 100, 20], 'Value', '50,50');
    l_Field = uieditfield(ImpedanceFig, 'text', 'Position', [200, 300, 100, 20], 'Value', '0.5,0.5');
    Il_Field = uieditfield(ImpedanceFig, 'text', 'Position', [200, 280, 100, 20], 'Value', '10,10');
    mm_Field = uieditfield(ImpedanceFig, 'text', 'Position', [200, 260, 100, 20], 'Value', '5,5');
    Im_Field = uieditfield(ImpedanceFig, 'text', 'Position', [200, 240, 100, 20], 'Value', '0.01,0.01');
    kr_Field = uieditfield(ImpedanceFig, 'text', 'Position', [200, 220, 100, 20], 'Value', '100,100');
    g_Field = uieditfield(ImpedanceFig, 'text', 'Position', [200, 200, 100, 20], 'Value', '9.81');
    Kd_Field = uieditfield(ImpedanceFig, 'text', 'Position', [200, 180, 100, 20], 'Value', '1600');
    Kp_Field = uieditfield(ImpedanceFig, 'text', 'Position', [200, 160, 100, 20], 'Value', '5000');
    Md_Field = uieditfield(ImpedanceFig, 'text', 'Position', [200, 140, 100, 20], 'Value', '100');
    t_Field = uieditfield(ImpedanceFig, 'text', 'Position', [200, 120, 100, 20], 'Value', '2.5');
    tc_Field = uieditfield(ImpedanceFig, 'text', 'Position', [200, 100, 100, 20], 'Value', '0.001');

    % button to update the equations with the new values
    runSimButton = uibutton(ImpedanceFig, 'Position', [120, 10, 100, 22], 'Text', 'Execute Values', 'ButtonPushedFcn', @(btn,event) runSimButton_Callback(...
        a_Field, ...
        ml_Field, ...
        l_Field, ...
        Il_Field, ...
        mm_Field, ...
        Im_Field, ...
        kr_Field, ...
        g_Field, ...
        Kd_Field, ...
        Kp_Field, ...
        Md_Field, ...
        t_Field, ...
        tc_Field));
    runSimButton.BackgroundColor = [0.2 0.7 0.2]; % Green color

    % show plots button
    % plotButton = uibutton(dynamicsFig, 'Position', [230, 10, 100, 22], 'Text', 'Plot Simulation', 'ButtonPushedFcn', @(btn,event) plotButton_Callback(dynamicsFig));

    % % Label for Simulink Buttons
    % uilabel(dynamicsFig, 'Text', 'Simulink Control Models', 'Position', [40, 100, 400, 22]);
    % % Buttons for Simulink Models by Tatwik, Danis, and Zane
    % % Indirect Force Control via Compliance
    % indirectForceButton = uibutton(dynamicsFig, 'Position', [10, 80, 200, 22], 'Text', 'Compliance Force Control', 'ButtonPushedFcn', @(btn,event) indirectForceButton_Callback( ...
    %     dynamicsFig ... % Put Params here or do something else to get them passed to the model
    %     ));
    % % Indirect Force Control via Impedance
    % indirectImpedanceButton = uibutton(dynamicsFig, 'Position', [10, 50, 200, 22], 'Text', 'Impedance Force Control', 'ButtonPushedFcn', @(btn,event) indirectImpedanceButton_Callback( ...
    %     dynamicsFig ... % Put Params here or do something else to get them passed to the model
    %     ));

end

function runSimButton_Callback( ...
    a_Field, ...
        ml_Field, ...
        l_Field, ...
        Il_Field, ...
        mm_Field, ...
        Im_Field, ...
        kr_Field, ...
        g_Field, ...
        Kd_Field, ...
        Kp_Field, ...
        Md_Field, ...
        t_Field, ...
        tc_Field)

    % Get the values from the user input fields
    % convert the EditField values to numbers str2num

    % convert chars to array
    a = str2num(a_Field.Value);
    m_l = str2num(ml_Field.Value);
    l = str2num(l_Field.Value);
    I_l = str2num(Il_Field.Value);
    m_m = str2num(mm_Field.Value);
    I_m = str2num(Im_Field.Value);
    k_r = str2num(kr_Field.Value);
    g = str2num(g_Field.Value);
    Kd = str2num(Kd_Field.Value);
    Kp = str2num(Kp_Field.Value);
    Md = str2num(Md_Field.Value);
    t_d = str2num(t_Field.Value);
    Tc = str2num(tc_Field.Value);
    
    % user input checks to make sure the dimensions are correct for the number of joints in the system
    NumberOfJoints = 2;
    if length(a) ~= NumberOfJoints
        errordlg('The number of gear ratios does not match the number of joints in the system.', 'Error', 'modal');
        return;
    elseif length(m_l) ~= NumberOfJoints
        errordlg('The number of link masses does not match the number of joints in the system.', 'Error', 'modal');
        return;
    elseif length(l) ~= NumberOfJoints
        errordlg('The number of motor masses does not match the number of joints in the system.', 'Error', 'modal');
        return;
    elseif length(I_l) ~= NumberOfJoints
        errordlg('The number of link inertias does not match the number of joints in the system.', 'Error', 'modal');
        return;
    elseif length(m_m) ~= NumberOfJoints
        errordlg('The number of motor inertias does not match the number of joints in the system.', 'Error', 'modal');
        return;
    elseif length(I_m) ~= NumberOfJoints
        errordlg('The number of static frictions does not match the number of joints in the system.', 'Error', 'modal');
        return;
    elseif length(k_r) ~= NumberOfJoints
        errordlg('The number of viscosity frictions does not match the number of joints in the system.', 'Error', 'modal');
        return;
    elseif length(g) ~= 1
        errordlg('The number of joint lengths does not match the number of joints in the system.', 'Error', 'modal');
        return;
    elseif length(Kd) ~= 1
        errordlg('The number of center of mass lengths does not match the number of joints in the system.', 'Error', 'modal');
        return;
    elseif length(Kp) ~= 1
        errordlg('The number of torques does not match the number of joints in the system.', 'Error', 'modal');
        return;
    elseif length(Md) ~= 1
        errordlg('The number of initial conditions does not match the number of joints in the system. x0 is defined as [position, velocity] for each joint.', 'Error', 'modal');
        return;
    elseif length(t_d) ~= 1
        errordlg('The number of initial conditions does not match the number of joints in the system. x0 is defined as [position, velocity] for each joint.', 'Error', 'modal');
        return;
    elseif length(Tc) ~= 1
        errordlg('The number of initial conditions does not match the number of joints in the system. x0 is defined as [position, velocity] for each joint.', 'Error', 'modal');
        return;
    end

% Impedance parameters code starts here
global a k_r1 k_r2 pi_m pi_l

% load manipulator dynamic parameters without load mass
  % link parameters
  % lengths
    % a  = [1;1];

  % masses
    % m_l1  = 50;
    % m_l2  = 50;
    m_l1 = m_l(1);
    m_l2 = m_l(2);


  % distances of centers of mass from joint axes
    % l_1  =  0.5;
    % l_2  =  0.5;
    l_1 = l(1);
    l_2 = l(2);

  % inertia moments relative to the centers of mass
    % I_l1 = 10;
    % I_l2 = 10;
    I_l1 = I_l(1);
    I_l2 = I_l(2);

% motor parameters
  % masses
    % m_m1 = 5;
    % m_m2 = 5;
    m_m1 = m_m(1);
    m_m2 = m_m(2);

  % inertia moments
    % I_m1 = 0.01;
    % I_m2 = 0.01;
    I_m1 = I_m(1);
    I_m2 = I_m(2);

  % transmission ratios
    % k_r1 = 100;
    % k_r2 = 100;
    k_r1 = k_r(1);
    k_r2 = k_r(2);

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
  % g = 9.81;

% friction matrix
  K_r = diag([k_r1 k_r2]);
  F_v = K_r*diag([0.01 0.01])*K_r;

% sample time of controller
  % Tc = 0.001;

% impedance controller gains
  % K_d = 1600*diag([1 1]);
  % K_p = 5000*diag([1 1]);
  % M_d = 100*diag([1 1]);
  % iM_d = inv(M_d);
    K_d = Kd*diag([1 1]);
    K_p = Kp*diag([1 1]);
    M_d = Md*diag([1 1]);
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
% cd("Impedance_Control\")
% initial joint configuration
  q_i = inv_k2u(a,p_i);

% duration of simulation
  % t_d = 2.5;

% trajectory in base frame
% time base vector
T = (0:Tc:t_d)';

% segment length
D_s = norm(p_f - p_i);

% trapezoidal velocity profile trajectory for path coordinate from 0 to 1
ds_c = 0.5;                % maximum velocity
t_f  = 1;                  % final time
[T_1,s,ds,dds,err] = trapez(0,1,ds_c/D_s,t_f,Tc);

n = size(T,1);
m = size(T_1,1);
o_d = zeros(n,2);
do_d = o_d;
ddo_d = o_d;

% tip trajectory
o_d(1:m,:) = ones(m,1)*p_i' + s*(p_f - p_i)';
o_d(m+1:n,:) = ones(n-m,1)*p_f';
do_d(1:m,:) = ds*(p_f - p_i)';
ddo_d(1:m,:) = dds*(p_f - p_i)';


% sample time for plots
Ts = Tc;

filename = "impedanceControl_ws.mat";
save(filename)

% sim('s9_3.mdl')
sim('Impedance_Control.slx')
    % plotButton_Callback(dynamicsFig,TransMats_Joint2Joint, q, JointSymbolic, JointLengths);


disp('Simulation complete');

e = ans.e;
f_e = ans.f_e;
time = ans.time;
pos_desired = ans.pos_desired;
vel_desired = ans.vel_desired;
pos_actual = ans.pos_actual;
vel_actual = ans.vel_actual;

hold off
clf

% position error in the base frame
  subplot(2,4,1)
  plot(time, e);
  axis([0 t_d -6e-2 6e-2]);
  %set(gca,'fontname','Times','fontsize',12,'fontweight','normal');
  xlabel('[s]');
  ylabel('[m]');
  title('Position error in the base frame');
  legend({'x', 'y'});

% contact force in the base frame
  subplot(2,4,2)
  plot(time, f_e);
  axis([0 t_d -550 550]);
  %set(gca,'fontname','Times','fontsize',12,'fontweight','normal');
  xlabel('[s]');
  ylabel('[N]');
  title('contact force in the base frame');
  legend({'x', 'y'});


% position error in the rotated base frame
  ec = e*R_c;
  subplot(2,4,3)
  plot(time, ec);
  axis([0 t_d -6e-2 6e-2]);
  set(gca,'fontname','Times','fontsize',12,'fontweight','normal');
  xlabel('[s]');
  ylabel('[m]');
  title('position error in the rotated base frame');
  legend({'x_c', 'y_c'});

% contact force in the rotated base frame
  f_ec = f_e*R_c;
  subplot(2,4,4)
  plot(time, f_ec);
  axis([0 t_d -550 550]);
  %set(gca,'fontname','Times','fontsize',12,'fontweight','normal');
  xlabel('[s]');
  ylabel('[N]');
  title('contact force in the rotated base frame');
  legend({'x_c', 'y_c'});


% desired position
  subplot(2,4,5)
  plot(time, pos_desired);
  %set(gca,'fontname','Times','fontsize',12,'fontweight','normal');
  xlabel('[s]');
  ylabel('[m]');
  title('Desired Position');
  legend({'x_c', 'y_c'});

% desired velocity
  subplot(2,4,6)
  plot(time, vel_desired);
  %set(gca,'fontname','Times','fontsize',12,'fontweight','normal');
  xlabel('[s]');
  ylabel('[m/s]');
  title('Desired Velocity');
  legend({'x_c', 'y_c'});

% Actual position
  subplot(2,4,7)
  plot(time, pos_actual);
  %set(gca,'fontname','Times','fontsize',12,'fontweight','normal');
  xlabel('[s]');
  ylabel('[m]');
  title('Actual Position');
  legend({'x_c', 'y_c'});

% Actual velocity
  subplot(2,4,8)
  plot(time, vel_actual);
  %set(gca,'fontname','Times','fontsize',12,'fontweight','normal');
  xlabel('[s]');
  ylabel('[m/s]');
  title('Actual Velocity');
  legend({'x_c', 'y_c'});
end

