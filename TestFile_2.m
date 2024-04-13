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

    tao(i) = subs(rhs(Tao(i)),[str2sym(Kr_Sym),str2sym(MassLI_Sym),str2sym(MassMI_Sym),str2sym(IntertiaLI_Sym),...
        str2sym(IntertiaMI_Sym),str2sym(JointSymbolic),str2sym(CenterOfMassSymbolic)],...
        [Kr_Val,MassLI_Val,MassMI_Val,IntertiaLI_Val,IntertiaMI_Val,JointLengths,CenterOfMassLengths]);
end

% send to workspace for ode45 stuff (Will remove later when no longer needed for debugging) :)
assignin('base', 'tao', tao);
assignin('base', 'NumberOfJoints', NumberOfJoints);
assignin('base', 'JointLengths', JointLengths);
assignin('base', 'q', q);
assignin('base', 'qd', qd);
assignin('base', 'qdd', qdd);

%% Solve it plz

% Maybe Make User Inputs :)
initial_t1 = 1; % Initial joint position for joint 1
initial_t2 = 1; % Initial joint position for joint 2
initial_t3 = 1; % Initial joint position for joint 3
initial_t1_d = 0; % Initial joint velocity for joint 1
initial_t2_d = 0; % Initial joint velocity for joint 2
initial_t3_d = 0; % Initial joint velocity for joint 3

% Define initial conditions for the state variables (joint positions, velocities, accelerations)
X0 = [initial_t1; initial_t2; initial_t3; initial_t1_d; initial_t2_d; initial_t3_d];

% Run the simulation up to the time of the step
[t1, X1] = ode45(@systemDynamics, [0,5], X0);

% Apply the step input to the joint variables
X_step = X1(end, :);  % get the final state from the first simulation
X_step(1:3) = X_step(1:3) + [1 1 1];  % apply the step input

% Continue the simulation from the time of the step with the new initial conditions
[t2, X2] = ode45(@systemDynamics, [5, 10], X_step);

% Combine the results from the two simulations
t = [t1; t2];
X = [X1; X2];

% Extract the state variables (joint positions, velocities, accelerations) from the solution
t1 = X(:, 1);
t2 = X(:, 2);
t3 = X(:, 3);
t1_d = X(:, 4);
t2_d = X(:, 5);
t3_d = X(:, 6);

% Plot or analyze the results as needed
% Plot the joint positions over time
figure;
plot(t, t1, 'r', t, t2, 'g', t, t3, 'b');
xlabel('Time (s)');
ylabel('Joint Position (rad)');
legend('Joint 1', 'Joint 2', 'Joint 3');
title('Joint Positions vs Time');

% Plot the joint velocities over time
figure;
plot(t, t1_d, 'r', t, t2_d, 'g', t, t3_d, 'b');
xlabel('Time (s)');
ylabel('Joint Velocity (rad/s)');
legend('Joint 1', 'Joint 2', 'Joint 3');
title('Joint Velocities vs Time');
end


