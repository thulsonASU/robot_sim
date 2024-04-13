%Inverse Kine
%Dynamics

%%Assumtions include
%All joints are revolute or Prismatic
%All links are straight

%% GUI
% Create a figure window
fig = uifigure('Name', 'DH Table GUI', 'Position', [100, 100, 550, 300]);

% Add text above the table
uilabel(fig, 'Text', 'For joint variables please enter d1...dn or t1...tn', 'Position', [20, 250, 400, 20]);

% Create a table for DH parameters
dhTable = uitable(fig, 'Position', [20, 20, 400, 220], 'ColumnName', {'a', 'alpha', 'd', 'theta'}, 'ColumnEditable', [true true true true]);

% Create buttons to add or remove rows
addButton = uibutton(fig, 'Position', [430, 220, 100, 22], 'Text', 'Add Row', 'ButtonPushedFcn', @(btn,event) addButton_Callback(dhTable));
removeButton = uibutton(fig, 'Position', [430, 190, 100, 22], 'Text', 'Remove Row', 'ButtonPushedFcn', @(btn,event) removeButton_Callback(dhTable));



%Create Gravity Dropdown Menu
uilabel(fig, 'Text', 'Select Gravity Direction', 'Position', [430, 160, 100, 22]);
dropdown = uidropdown(fig, 'Position', [430, 140, 100, 22], 'Items', {'x', 'y', 'z', '-x', '-y', '-z'}, 'Value', 'x', 'ValueChangedFcn', @(dropdown,event) dropdown_Callback(dropdown));
selectedOption = 'x';


% Create a "Run" button
runButton = uibutton(fig, 'Position', [430, 110, 100, 22], 'Text', 'Run', 'ButtonPushedFcn', @(btn,event) runButton_Callback(dhTable,dropdown));
runButton.BackgroundColor = [0.2 0.7 0.2]; % Green color
data = {'','','',''};
dhTable.Data = data;



function dropdown_Callback(dropdown)
selectedOption = dropdown.Value;
disp(['Selected Gravity Direction: ' selectedOption]);
end

function addButton_Callback(dhTable)
% Add a row to the DH table
data = dhTable.Data;
newRow = {'', '', '', ''};
data = [data; newRow];
dhTable.Data = data;
end

function removeButton_Callback(dhTable)
% Remove the last row from the DH table

data = dhTable.Data;
if size(data, 1) > 1
    data(end, :) = [];
    dhTable.Data = data;
else
    errordlg('Cannot remove more rows.', 'Error', 'modal');
end
end
%% Run Function
function runButton_Callback(dhTable,dropdown)
clc;
disp('Running...');
GravityDirection = dropdown.Value;
DHTable = dhTable.Data;
%disp('DH Parameters:');
%disp(DHTable);
%Find Number of Joints
SizeOfData = size(DHTable);
NumberOfJoints = SizeOfData(1);

%Find type of joint
JointLengths = [];
JointSymbolic = string([]);
for i = 1:NumberOfJoints
    
    if contains(DHTable{i,4},"t")
        JointTypes(i) = "r";
        q(i) = "t"+i;
        DHTable{i,4} = str2sym(convertCharsToStrings(DHTable{i,4}));
        
        if contains(DHTable{i,2},"pi")
            DHTable{i,2} = str2sym(convertCharsToStrings(DHTable{i,2}));
        else
            DHTable{i,2} = str2num(DHTable{i,2});
        end
        DHTable{i,3} = str2num(DHTable{i,3});
        DHTable{i,1} = str2num(DHTable{i,1});
    elseif contains(DHTable{i,3},"d")
        JointTypes(i) = "p";
        q(i) = "d"+i;
        DHTable{i,3} = str2sym(convertCharsToStrings(DHTable{i,3}));
        
        if contains(DHTable{i,4},"pi")
            DHTable{i,4} = str2sym(convertCharsToStrings(DHTable{i,4}));
        else
            DHTable{i,4} = str2num(DHTable{i,4});
        end
        if contains(DHTable{i,2},"pi")
            DHTable{i,2} = str2sym(convertCharsToStrings(DHTable{i,2}));
        else
            DHTable{i,2} = str2num(DHTable{i,2});
        end
        DHTable{i,1} = str2num(DHTable{i,1});
    end
    
    %Fixing some joint calculation things, I think....
    if DHTable{i,1} == 0
        if JointTypes(i) == "p"
        else
            JointLengths(i) = DHTable{i,3};
            si = string(i);
            JointSymbolic(i) = "a"+si;
            
            DHTable{i,3} = str2sym(JointSymbolic(i));
        end
    elseif DHTable{i,3} == 0
        JointLengths(i) = DHTable{i,1};
        si = string(i);
        JointSymbolic(i) = "a"+si;
        
        DHTable{i,1} = str2sym(JointSymbolic(i));
    end
end
fprintf("Robot has " + NumberOfJoints + " Joints, in the configuration of:")
disp(JointTypes)

% send joint lengths to workspace :)
assignin('base', 'JointLengths', JointLengths);

%q
%Find all the transformation Matricies:
for i = 1:NumberOfJoints
    %Please note: This is just solving from joint to joint not base to
    %joint
    TransMats_Joint2Joint(:,:,i) = ...
        DHSolver(DHTable{i,1},DHTable{i,2},DHTable{i,3},DHTable{i,4});
end
for i = 1:NumberOfJoints
    % Multiply elements up to index i
    product = TransMats_Joint2Joint(:,:,1);
    for j = 2:i
        product = product * TransMats_Joint2Joint(:,:,j);
    end
    % Store the result at index i
    TransMats(:,:,i) = simplify(product);
end
%TransMats is a 3d matrix that contains T01 all the way to T0E - E being
%the number of the end effector.
%When you call TransMats and want, for example Z1, you would input it as
%such.
%TransMats(1:3,3,1) as 1 is the first matrix (T01) and we want Z1
fprintf("The translational Matrix for T0"+NumberOfJoints+" symbolically:\n\n")
disp(TransMats(:,:,NumberOfJoints))


fprintf("The translational Matrix for T0"+NumberOfJoints+" numerically:\n\n")
disp(subs(TransMats(:,:,NumberOfJoints),str2sym(JointSymbolic),JointLengths))

%% Inverse Jacs
for i = 1:NumberOfJoints
    %Definging My Kr's
    Kr_Val(i) = 1;
    Kr_Sym(i) = "Kr"+i;
    %Assume Kr is positive :)
    assume(str2sym(Kr_Sym(i)),"positive")
    
    
    MassLI_Val(i) = 1;
    MassLI_Sym(i) = "M_l"+i;
    MassMI_Val(i) = 1;
    MassMI_Sym(i) = "M_m"+i;
    
    IntertiaLI_Val(i) = 1;
    IntertiaLI_Sym(i) = "I_l"+i;
    IntertiaMI_Val(i) = 1;
    IntertiaMI_Sym(i) = "I_m"+i;
    %This will change with time lol it will be GUI I promise
    
    %Assume joints are real.
    assume(str2sym(q(i)), "real");
    %Assume joints lengths are positive lol
    if numel(JointSymbolic) > 0
        assume(str2sym(JointSymbolic(i)),"positive")
    end
end

for i = 1:NumberOfJoints
    if JointTypes(i) == "p"
        CenterOfMassLengths(i) = 0;
        CenterOfMassSymbolic(i) = "";
    else
        CenterOfMassLengths(i) = JointLengths(i)/2;
        CenterOfMassSymbolic(i) = "L"+i;
    end
    %assume center off masses are always positive :)
    assume(str2sym(CenterOfMassSymbolic(i)),"positive")
end

for i = 1:NumberOfJoints
    JPLi(:,:,i) = Jpj_Li(i,NumberOfJoints,JointTypes,TransMats,str2sym(JointSymbolic),str2sym(CenterOfMassSymbolic));
    JOLi(:,:,i) = Joj_Li(i,NumberOfJoints,JointTypes,TransMats);
    
    JPMi(:,:,i) = Jpj_Mi(i,NumberOfJoints,JointTypes,TransMats);
    JOMi(:,:,i) = Joj_Mi(i,NumberOfJoints,JointTypes,TransMats,str2sym(Kr_Sym(i)));
end

%JPLi
%JOLi
%JPMi
%JOMi

%% Calculating B
% (M_LI,M_MI,JP_LI,JO_LI,I_LI,I_MI,Trans,JP_MI,JO_MI,Num)
Bq = FindB(MassLI_Sym,MassMI_Sym,JPLi,JOLi,IntertiaLI_Sym,...
    IntertiaMI_Sym,TransMats,JPMi,JOMi,NumberOfJoints);

%Bq
%% Calculating B
% (Bq,q,Num)
C = FindC(Bq,q,NumberOfJoints);
%C
%% Calculating Gravity
% (GDir, M_L,M_M,JP_L,JP_M,Num,g)
G = FindG(GravityDirection,MassLI_Sym,MassMI_Sym,JPLi,JPMi,NumberOfJoints,9.81);
%G

%% Friction 
% Assuming 0 friction  :)

%% Adding them all together :)


%Make qd and qdd and taos 
for i = 1:NumberOfJoints
    qd = q+"_d";
    qdd = qd+"d";
    tao(i,1) = str2sym("Tao_"+i);
end
Bq_x_qdd = Bq*str2sym(qdd');
C_x_qd = C*str2sym(qd');
%Fv*qd
%Fs*sgn(q)

Tao = tao == Bq_x_qdd+C_x_qd+G;
fprintf("\nEquations Of Motion:\n")
for i = 1:NumberOfJoints
    fprintf("Equation #"+i+"\n")
    disp(Tao(i))

    % send to workspace
    assignin('base', 'Tao', Tao);
    assignin('base', 'NumberOfJoints', NumberOfJoints);
end

%% Solve it plz

% Make these user input later :)
syms Kr [1 NumberOfJoints] real
syms I_l [1 NumberOfJoints] real
syms I_m [1 NumberOfJoints] real
syms M_l [1 NumberOfJoints] real
syms M_m [1 NumberOfJoints] real
syms L [1 NumberOfJoints] real
syms tao [1 NumberOfJoints] real

% Simulate the system dynamics for 10 seconds
% Dynamic equations for the system are Tao(i) where i is the joint number
% Initial conditions for the system are q(0) and qd(0)
% The system is simulated using the ode45 function
% The simulation time is 10 seconds
% The initial conditions are q(0) = [0, 0, 0, 0, 0] and qd(0) = [0, 0, 0, 0, 0]
% The simulation time is 10 seconds

% Initial conditions
q0 = zeros(1, NumberOfJoints);
q_dot0 = zeros(1, NumberOfJoints);
q_ddot0 = zeros(1, NumberOfJoints);

temp_input = ones(1, NumberOfJoints);
% find all user inputs not yet defined and assign them to 1 :)
for i = 1:NumberOfJoints
    tao(i) = rhs(Tao(i)); % Get right hand side of the symbolic equation
    tao(i) = subs(tao(i), [Kr(i), I_l(i), I_m(i), M_l(i), M_m(i), L(i)], [temp_input(i),temp_input(i),temp_input(i),temp_input(i),1,JointLengths(i)]);
    disp(tao(i))
end

% Equations Of Motion:
% Equation #1
% Tao_1 == t3_dd*(I_l3 + I_m3*Kr3 + L3*M_l3*(L3 + a1*cos(t2 + t3) + a2*cos(t3))) - t1_d*(a1*t2_d*(L3*M_l3*sin(t2 + t3) + L2*M_l2*sin(t2) + M_l3*a2*sin(t2) + M_m3*a2*sin(t2)) + L3*M_l3*t3_d*(a1*sin(t2 + t3) + a2*sin(t3))) - t2_d*(a1*t1_d*(L3*M_l3*sin(t2 + t3) + L2*M_l2*sin(t2) + M_l3*a2*sin(t2) + M_m3*a2*sin(t2)) + a1*t2_d*(L3*M_l3*sin(t2 + t3) + L2*M_l2*sin(t2) + M_l3*a2*sin(t2) + 2*M_m3*a2*sin(t2)) + L3*M_l3*t3_d*(a1*sin(t2 + t3) + a2*sin(t3))) + t1_dd*(I_l1 + I_l2 + I_l3 + I_m2 + I_m3 + I_m1*Kr1^2 + L1^2*M_l1 + L2^2*M_l2 + L3^2*M_l3 + M_l2*a1^2 + M_l3*a1^2 + M_l3*a2^2 + M_m2*a1^2 + M_m3*a1^2 + M_m3*a2^2 + 2*L3*M_l3*a1*cos(t2 + t3) + 2*L2*M_l2*a1*cos(t2) + 2*L3*M_l3*a2*cos(t3) + 2*M_l3*a1*a2*cos(t2) + 2*M_m3*a1*a2*cos(t2)) + t2_dd*(I_l2 + I_l3 + I_m3 + L2^2*M_l2 + L3^2*M_l3 + M_l3*a2^2 + M_m3*a1^2 + M_m3*a2^2 + I_m2*Kr2 + L3*M_l3*a1*cos(t2 + t3) + L2*M_l2*a1*cos(t2) + 2*L3*M_l3*a2*cos(t3) + M_l3*a1*a2*cos(t2) + 2*M_m3*a1*a2*cos(t2)) - L3*M_l3*t3_d*(a1*sin(t2 + t3) + a2*sin(t3))*(t1_d + t2_d + t3_d)
 
% Equation #2
% Tao_2 == t1_d*(t1_d*(L3*M_l3*a1*sin(t2 + t3) + L2*M_l2*a1*sin(t2) + M_l3*a1*a2*sin(t2) + M_m3*a1*a2*sin(t2)) - L3*M_l3*a2*t3_d*sin(t3)) + t2_dd*(I_l2 + I_l3 + I_m3 + I_m2*Kr2^2 + L2^2*M_l2 + M_l3*(L3^2 + 2*cos(t3)*L3*a2 + a2^2) + M_m3*(a1^2 + 2*cos(t2)*a1*a2 + a2^2)) + t3_dd*(I_l3 + I_m3*Kr3 + L3*M_l3*(L3 + a2*cos(t3))) - t2_d*(L3*M_l3*a2*t3_d*sin(t3) + M_m3*a1*a2*t2_d*sin(t2)) + t1_dd*(I_l2 + I_l3 + I_m3 + L2^2*M_l2 + L3^2*M_l3 + M_l3*a2^2 + M_m3*a1^2 + M_m3*a2^2 + I_m2*Kr2 + L3*M_l3*a1*cos(t2 + t3) + L2*M_l2*a1*cos(t2) + 2*L3*M_l3*a2*cos(t3) + M_l3*a1*a2*cos(t2) + 2*M_m3*a1*a2*cos(t2)) - L3*M_l3*a2*t3_d*sin(t3)*(t1_d + t2_d + t3_d)
 
% Equation #3
% Tao_3 == t1_dd*(I_l3 + I_m3*Kr3 + L3*M_l3*(L3 + a1*cos(t2 + t3) + a2*cos(t3))) + t3_dd*(I_m3*Kr3^2 + M_l3*L3^2 + I_l3) + t1_d*(t1_d*(L3*M_l3*a1*sin(t2 + t3) + L3*M_l3*a2*sin(t3)) + L3*M_l3*a2*t2_d*sin(t3)) + t2_dd*(I_l3 + I_m3*Kr3 + L3*M_l3*(L3 + a2*cos(t3))) + L3*M_l3*a2*t2_d*sin(t3)*(t1_d + t2_d)
 
%   [TOUT,YOUT] = ODE45(ODEFUN,TSPAN,Y0) integrates the system of
%   differential equations y' = f(t,y) from time TSPAN(1) to TSPAN(end)
%   with initial conditions Y0. Each row in the solution array YOUT
%   corresponds to a time in the column vector TOUT. 
%     * ODEFUN is a function handle. For a scalar T and a vector Y,
%       ODEFUN(T,Y) must return a column vector corresponding to f(t,y).
%     * TSPAN is a two-element vector [T0 TFINAL] or a vector with
%       several time points [T0 T1 ... TFINAL]. If you specify more than
%       two time points, ODE45 returns interpolated solutions at the
%       requested times.
%     * YO is a column vector of initial conditions, one for each equation.

% Simulate the system dynamics
% [t, q] = ode45(@(t, q) systemDynamics(t, q, Tao, NumberOfJoints), [0, 10], [q0, q_dot0]);


end



