% Inverse Kine
% Dynamics

%%Assumtions include
% All joints are revolute or Prismatic
% All links are straight

%% GUI
% Create a figure window
fig = uifigure('Name', 'Denavit-Hartenberg Table GUI', 'Position', [100, 100, 550, 300]);
% Add text above the table
uilabel(fig, 'Text', 'For joint variables please enter d1...dn or t1...tn', 'Position', [20, 250, 400, 20]);

% Create a button to open a figure window that references the DH parameters
DHfigButton = uibutton(fig, 'Position', [430, 80, 100, 22], 'Text', 'DH Reference', 'ButtonPushedFcn', @(btn,event) DHfigButton_Callback());

% Create a table for DH parameters
dhTable = uitable(fig, 'Position', [20, 20, 400, 220], 'ColumnName', {'a (m)', 'alpha (rad)', 'd (m)', 'theta (rad)'}, 'ColumnEditable', [true true true true]);

% Create buttons to add or remove rows
addButton = uibutton(fig, 'Position', [430, 220, 100, 22], 'Text', 'Add Row', 'ButtonPushedFcn', @(btn,event) addButton_Callback(dhTable));
removeButton = uibutton(fig, 'Position', [430, 190, 100, 22], 'Text', 'Remove Row', 'ButtonPushedFcn', @(btn,event) removeButton_Callback(dhTable));

%Create Gravity Dropdown Menu
uilabel(fig, 'Text', 'Select Gravity Direction', 'Position', [430, 160, 100, 22]);
dropdown = uidropdown(fig, 'Position', [430, 140, 100, 22], 'Items', {'x', 'y', 'z', '-x', '-y', '-z'}, 'Value', '-z', 'ValueChangedFcn', @(dropdown,event) dropdown_Callback(dropdown));
% selectedOption = 'x';

% Create a "Run" button
runButton = uibutton(fig, 'Position', [430, 110, 100, 22], 'Text', 'Run', 'ButtonPushedFcn', @(btn,event) runButton_Callback(dhTable,dropdown));
runButton.BackgroundColor = [0.2 0.7 0.2]; % Green color
data = {'','','',''};
dhTable.Data = data;

closeButton = uibutton(fig, 'Position', [430, 10, 100, 22], 'Text', 'Close', 'ButtonPushedFcn', @(btn,event) close(fig));
closeButton.BackgroundColor = [0.7 0.2 0.2]; % Red color

% % BUTTON FOR DYNAMICS SIMULATION
% dynamicsButton = uibutton(fig, 'Position', [430, 80, 100, 22], 'Text', 'Dynamics Sim', 'ButtonPushedFcn', ...
%     @(btn,event) dynamicsButton_Callback( ...
%     evalin('base', 'NumberOfJoints'), ...
%     evalin('base', 'Tao'), ...
%     evalin('base', 'Kr_Sym'), ...
%     evalin('base', 'MassLI_Sym'), ...
%     evalin('base', 'MassMI_Sym'), ...
%     evalin('base', 'IntertiaLI_Sym'), ...
%     evalin('base', 'IntertiaMI_Sym'), ...
%     evalin('base', 'FrictionS_Sym'), ...
%     evalin('base', 'FrictionV_Sym'), ...
%     evalin('base', 'JointSymbolic'), ...
%     evalin('base', 'CenterOfMassSymbolic'), ...
%     evalin('base', 'JointLengths'), ...
%     evalin('base', 'CenterOfMassLengths'),...
%     evalin('base', 'q'), ...
%     evalin('base', 'qd'), ...
%     evalin('base', 'qdd')));

function DHfigButton_Callback()
    % Create a figure window
    DHfig = uifigure('Name', 'Denavit-Hartenberg Reference', 'Position', [100, 100, 568, 417]);
    % insert png of DH parameters into fig MATLAB
    % Create an axes object
    ax = axes(DHfig, 'Units', 'pixels', 'Position', [0, 0, 568, 417]); % replace x, y, width, height with the desired position and size
    % Load the image
    img = imread('./images/DH.png'); % replace 'image.png' with the path to your image
    % Display the image
    imshow(img, 'Parent', ax);
    % close
    closeButton = uibutton(DHfig, 'Position', [10, 10, 100, 22], 'Text', 'Close', 'ButtonPushedFcn', @(btn,event) close(DHfig));
    closeButton.BackgroundColor = [0.7 0.2 0.2]; % Red color
end

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
                JointSymbolic(i) = "d"+si;
                
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
        if numel(JointSymbolic) >= i
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
     %simplify(G)

    %% Friction 
    for i = 1:NumberOfJoints
        FrictionS_Sym(i) = "Fs"+i;
        FrictionS_Val(i) = 1;
        
        FrictionV_Sym(i) = "Fv"+i;
        FrictionV_Val(i) = 1;
        
    end
    %% Adding them all together :)

    %Make qd and qdd and taos 
    for i = 1:NumberOfJoints
        qd = q+"_d";
        qdd = qd+"d";
        tao(i,1) = str2sym("Tao_"+i);
    end
    Bq_x_qdd = Bq*str2sym(qdd');
    C_x_qd = C*str2sym(qd');
    for i = 1:NumberOfJoints
        Fs(i,:) = str2sym(FrictionS_Sym(i))*str2sym(qd(i));
        Fv(i,:) = str2sym(FrictionV_Sym(i)); % * sgn(qd) % Inside the simulation a piecewise is added to get sign of qd
    end
    Tao = tao == Bq_x_qdd+C_x_qd+Fs+Fv+G;
    fprintf("\nEquations Of Motion:\n")
    for i = 1:NumberOfJoints
        fprintf("Equation #"+i+"\n")
        disp(simplify(Tao(i)))

    end

    % preallocate Q_Val and R_Val in workspace as empty
    Q_Val = [];
    R_Val = [];

    % send to workspace for ode45 stuff
    assignin('base', 'Tao', Tao);
    assignin('base', 'NumberOfJoints', NumberOfJoints);
    assignin('base', 'JointSymbolic', JointSymbolic);
    assignin('base', 'CenterOfMassSymbolic', CenterOfMassSymbolic);
    assignin('base', 'Kr_Sym', Kr_Sym);
    assignin('base', 'MassLI_Sym', MassLI_Sym);
    assignin('base', 'MassMI_Sym', MassMI_Sym);
    assignin('base', 'IntertiaLI_Sym', IntertiaLI_Sym);
    assignin('base', 'IntertiaMI_Sym', IntertiaMI_Sym);
    assignin('base', 'FrictionS_Sym', FrictionS_Sym);
    assignin('base', 'FrictionV_Sym', FrictionV_Sym);
    assignin('base', 'JointLengths', JointLengths);
    assignin('base', 'CenterOfMassLengths', CenterOfMassLengths);
    assignin('base', 'q', q);
    assignin('base', 'qd', qd);
    assignin('base', 'qdd', qdd);
    assignin('base', 'Q_Val', Q_Val);
    assignin('base', 'R_Val', R_Val);
    assignin('base', 'TransMats_Joint2Joint', TransMats_Joint2Joint);

    % call the dynamics button callback
    dynamicsButton_Callback(...
    evalin('base', 'TransMats_Joint2Joint'), ...
    evalin('base', 'NumberOfJoints'), ...
    evalin('base', 'Tao'), ...
    evalin('base', 'Kr_Sym'), ...
    evalin('base', 'MassLI_Sym'), ...
    evalin('base', 'MassMI_Sym'), ...
    evalin('base', 'IntertiaLI_Sym'), ...
    evalin('base', 'IntertiaMI_Sym'), ...
    evalin('base', 'FrictionS_Sym'), ...
    evalin('base', 'FrictionV_Sym'), ...
    evalin('base', 'JointSymbolic'), ...
    evalin('base', 'CenterOfMassSymbolic'), ...
    evalin('base', 'JointLengths'), ...
    evalin('base', 'CenterOfMassLengths'),...
    evalin('base', 'q'), ...
    evalin('base', 'qd'), ...
    evalin('base', 'qdd'));

end

% ============================================
% = TYLERS DYNAMICAL WONDERLAND STARTS BELOW =
% ============================================

function dynamicsButton_Callback(TransMats_Joint2Joint, NumberOfJoints, Tao, Kr_Sym, MassLI_Sym, MassMI_Sym, IntertiaLI_Sym, IntertiaMI_Sym, FrictionS_Sym, FrictionV_Sym, JointSymbolic, CenterOfMassSymbolic, JointLengths, CenterOfMassLengths, q, qd, qdd)
    % Create a new figure window
    dynamicsFig = uifigure('Name', 'Dynamics Simulation', 'NumberTitle', 'off');

    closeButton = uibutton(dynamicsFig, 'Position', [10, 10, 100, 22], 'Text', 'Close', 'ButtonPushedFcn', @(btn,event) close(dynamicsFig));
    closeButton.BackgroundColor = [0.7 0.2 0.2]; % Red color

    % Add UI controls to the new figure window as needed
    % Add labels for instructions
    uilabel(dynamicsFig, 'Text', 'Each comma seperated value corresponds to each joint in series. ', 'Position', [20, 395, 400, 22]);
    uilabel(dynamicsFig, 'Text', 'Provide the Values for the Following. Note: x0 position or velocity must be non zero ', 'Position', [40, 360, 500, 22]);
    
    % Add labels for each field
    uilabel(dynamicsFig, 'Text', 'Gear Ratios:', 'Position', [10, 340, 300, 20]);
    uilabel(dynamicsFig, 'Text', 'Mass of Links (kg):', 'Position', [10, 320, 300, 20]);
    uilabel(dynamicsFig, 'Text', 'Mass of Motors (kg):', 'Position', [10, 300, 300, 20]);
    uilabel(dynamicsFig, 'Text', 'Link Intertias (kg·m^2):', 'Position', [10, 280, 300, 20]);
    uilabel(dynamicsFig, 'Text', 'Motor Inertias (kg·m^2):', 'Position', [10, 260, 300, 20]);
    uilabel(dynamicsFig, 'Text', 'Static Frictions (N):', 'Position', [10, 240, 300, 20]);
    uilabel(dynamicsFig, 'Text', 'Viscosity Frictions (N·s/m^2):', 'Position', [10, 220, 300, 20]);

    uilabel(dynamicsFig, 'Text', 'Define u (initial torque):', 'Position', [280, 340, 200, 20]);
    uilabel(dynamicsFig, 'Text', 'Define x0 (pos and vel in rads):', 'Position', [280, 320, 200, 20]);
    uilabel(dynamicsFig, 'Text', 'Final simulation time (seconds):','Position', [280, 300, 300, 20]);

    % ODE45 Options
    uilabel(dynamicsFig, 'Text', 'Options for ODE45 Solver','Position', [40, 190, 200, 20]);
    uilabel(dynamicsFig, 'Text', 'Max Time Step:','Position', [10, 170, 100, 20]);
    uilabel(dynamicsFig, 'Text', 'Relative Tolerance:','Position', [10, 150, 120, 20]);
    uilabel(dynamicsFig, 'Text', 'Absolute Tolerance:','Position', [10, 130, 120, 20]);

    % Make a ui field to fill in for Kr MassLI MassMI IntertiaLI IntertiaMI center of mass lengths
    KrEditField = uieditfield(dynamicsFig, 'text', 'Position', [170, 340, 100, 20], 'Value', '1,1,1');
    MassLIEditField = uieditfield(dynamicsFig, 'text', 'Position', [170, 320, 100, 20], 'Value', '2,2,2');
    MassMIEditField = uieditfield(dynamicsFig, 'text', 'Position', [170, 300, 100, 20], 'Value', '2,2,2');
    IntertiaLIEditField = uieditfield(dynamicsFig, 'text', 'Position', [170, 280, 100, 20], 'Value', '1,1,1');
    IntertiaMIEditField = uieditfield(dynamicsFig, 'text', 'Position', [170, 260, 100, 20], 'Value', '1,1,1');
    FrictionSEditField = uieditfield(dynamicsFig, 'text', 'Position', [170, 240, 100, 20], 'Value', '0.05,0.05,0.05');
    FrictionVEditField = uieditfield(dynamicsFig, 'text', 'Position', [170, 220, 100, 20], 'Value', '0.05,0.05,0.05');

    uEditField = uieditfield(dynamicsFig, 'text', 'Position', [455, 340, 80, 20], 'Value', '0,0,0');
    x0EditField = uieditfield(dynamicsFig, 'text', 'Position', [455, 320, 80, 20], 'Value', '0,0,0,0,0,0');
    tf = uieditfield(dynamicsFig, 'text', 'Position', [455, 300, 50, 20], 'Value', '4');

    maxStepField = uieditfield(dynamicsFig, 'text', 'Position', [120, 170, 100, 20], 'Value', '0.1');
    relTolField = uieditfield(dynamicsFig, 'text', 'Position', [120, 150, 100, 20], 'Value', '1e-4');
    absTolField = uieditfield(dynamicsFig, 'text', 'Position', [120, 130, 100, 20], 'Value', '1e-5');

    % Checkbox so user can choose to use LQR control law for ODE45 solver
    lqrCheckBox = uicheckbox(dynamicsFig, 'Text', 'Use LQR Control Law', 'Position', [280, 260, 150, 20]);
    lqrEditButton = uibutton(dynamicsFig, 'Position', [440, 260, 100, 20], 'Text', 'Edit LQR Gains', 'ButtonPushedFcn', @(btn,event) lqrEditButton_Callback(dynamicsFig));

    % button to update the equations with the new values
    runSimButton = uibutton(dynamicsFig, 'Position', [120, 10, 100, 22], 'Text', 'Run Simulation', 'ButtonPushedFcn', @(btn,event) runSimButton_Callback(...
        TransMats_Joint2Joint,...
        dynamicsFig, ...    
        NumberOfJoints, ...
        Tao, ...
        Kr_Sym, ...
        MassLI_Sym, ...
        MassMI_Sym, ...
        IntertiaLI_Sym, ...
        IntertiaMI_Sym, ...
        FrictionS_Sym, ...
        FrictionV_Sym, ...
        JointSymbolic, ...
        CenterOfMassSymbolic, ...
        KrEditField, ...
        MassLIEditField, ...
        MassMIEditField, ...
        IntertiaLIEditField, ...
        IntertiaMIEditField, ...
        FrictionSEditField, ...
        FrictionVEditField, ...
        JointLengths, ...
        CenterOfMassLengths,...
        q, ...
        qd, ...
        qdd, ...
        uEditField, ...
        x0EditField, ...
        tf, ...
        maxStepField, ...
        relTolField, ...
        absTolField, ...
        lqrCheckBox, ...
        evalin('base', 'Q_Val'), ...
        evalin('base', 'R_Val')));
    runSimButton.BackgroundColor = [0.2 0.7 0.2]; % Green color

    % show plots button
    % plotButton = uibutton(dynamicsFig, 'Position', [230, 10, 100, 22], 'Text', 'Plot Simulation', 'ButtonPushedFcn', @(btn,event) plotButton_Callback(dynamicsFig));

    % Label for Simulink Buttons
    uilabel(dynamicsFig, 'Text', 'Simulink Control Models', 'Position', [40, 100, 400, 22]);
    % Buttons for Simulink Models by Tatwik, Danis, and Zane
    % Indirect Force Control via Compliance
    indirectForceButton = uibutton(dynamicsFig, 'Position', [10, 80, 200, 22], 'Text', 'Compliance Force Control', 'ButtonPushedFcn', @(btn,event) indirectForceButton_Callback( ...
        dynamicsFig ... % Put Params here or do something else to get them passed to the model
        ));
    % Indirect Force Control via Impedance
    indirectImpedanceButton = uibutton(dynamicsFig, 'Position', [10, 50, 200, 22], 'Text', 'Impedance Force Control', 'ButtonPushedFcn', @(btn,event) indirectImpedanceButton_Callback( ...
        dynamicsFig ... % Put Params here or do something else to get them passed to the model
        ));
end

% Add Button Callback for Indirect Force Control via Compliance in Simulink
function indirectForceButton_Callback(dynamicsFig)
    disp('WIP') % Zane and Gang
    % PUT CODE HERE FOR SIMULINK MODEL
end

% Add Button Callback for Indirect Force Control via Impedance in Simulink
function indirectImpedanceButton_Callback(dynamicsFig)
    disp('WIP') % Tatwik and Danis
    % PUT CODE HERE FOR SIMULINK MODEL
    % Open Simulink model
    open_system('Impedance_Control.slx');
    
    % Run Simulink simulation
    sim('Impedance_Control.slx');
%     M_d = 100*diag([]);
%     iM_d = inv(M_d)
% imp_Kd
% imp_Kp

        for j = 1:length(vars_in_eqs)
            %subs(Tao(i), )
        end
    sim('Impedance_Control',tspan)
end

% lqrEditButton callback function
function lqrEditButton_Callback(dynamicsFig)
    % Create a new figure window for editing LQR gains
    lqrFig = uifigure('Name', 'Edit LQR Gains', 'Position', [200, 200, 400, 200]);

    % Add text above the table
    uilabel(lqrFig, 'Text', 'Enter the LQR Gains for the Control Law', 'Position', [20, 170, 300, 20]);
    uilabel(lqrFig, 'Text', 'Assume velocity is five times more important than torque', 'Position', [20, 150, 400, 20]);

    uilabel(lqrFig, 'Text', 'How important are the states? (Q):', 'Position', [20, 120, 300, 20]);
    uilabel(lqrFig, 'Text', 'How important is energy consumption? (R):', 'Position', [20, 100, 300, 20]);
    QEditField = uieditfield(lqrFig, 'text', 'Position', [255, 120, 100, 20], 'Value', '500');
    REditField = uieditfield(lqrFig, 'text', 'Position', [255, 100, 100, 20], 'Value', '100');

    % Add a button to close the LQR gains figure
    closeButton = uibutton(lqrFig, 'Position', [10, 10, 100, 22], 'Text', 'Close', 'ButtonPushedFcn', @(btn,event) close(lqrFig));
    closeButton.BackgroundColor = [0.7 0.2 0.2]; % Red color

    % Add a button to update the LQR gains
    updateButton = uibutton(lqrFig, 'Position', [120, 10, 100, 22], 'Text', 'Update Gains', 'ButtonPushedFcn', @(btn,event) updateButton_Callback(QEditField,REditField));
    updateButton.BackgroundColor = [0.2 0.7 0.2]; % Green color

end

% callback function for updating the LQR gains
function updateButton_Callback(QEditField,REditField)
    Q_Val = str2num(QEditField.Value);
    R_Val = str2num(REditField.Value);

    % store the values in the dynamicsFig app data
    assignin('base', 'Q_Val', Q_Val);
    assignin('base', 'R_Val', R_Val);
end

% % substitute values
% eqs(i) = subs(rhs(Tao(i)),[str2sym(Kr_Sym),str2sym(MassLI_Sym),str2sym(MassMI_Sym),str2sym(IntertiaLI_Sym),...
%     str2sym(IntertiaMI_Sym),str2sym(JointSymbolic),str2sym(CenterOfMassSymbolic)],...
%     [Kr_Val,MassLI_Val,MassMI_Val,IntertiaLI_Val,IntertiaMI_Val,JointLengths,CenterOfMassLengths]);
% function callback for substituting the values into the equations using a button
function runSimButton_Callback( ...
    TransMats_Joint2Joint, ...
    dynamicsFig, ...
    NumberOfJoints, ...
    Tao, ...
    Kr_Sym, ...
    MassLI_Sym, ...
    MassMI_Sym, ...
    IntertiaLI_Sym, ...
    IntertiaMI_Sym, ...
    FrictionS_Sym, ...
    FrictionV_Sym, ...
    JointSymbolic, ...
    CenterOfMassSymbolic, ...
    KrEditField, ...
    MassLIEditField, ...
    MassMIEditField, ...
    IntertiaLIEditField, ...
    IntertiaMIEditField, ...
    FrictionSEditField, ...
    FrictionVEditField, ...
    JointLengths, ...
    CenterOfMassLengths,...
    q, ...
    qd, ...
    qdd, ...
    uEditField, ...
    x0EditField, ...
    tf, ...
    maxStepField, ...
    relTolField, ...
    absTolField, ...
    lqrCheckBox, ...
    Q_Val, ...
    R_Val ...
    )

    % Get the values from the user input fields
    % convert the EditField values to numbers str2num

    % convert chars to array
    Kr_Val = str2num(KrEditField.Value);
    MassLI_Val = str2num(MassLIEditField.Value);
    MassMI_Val = str2num(MassMIEditField.Value);
    IntertiaLI_Val = str2num(IntertiaLIEditField.Value);
    IntertiaMI_Val = str2num(IntertiaMIEditField.Value);
    FrictionS_Val = str2num(FrictionSEditField.Value);
    FrictionV_Val = str2num(FrictionVEditField.Value);
    maxStep_Val = str2num(maxStepField.Value);
    relTol_Val = str2num(relTolField.Value);
    absTol_Val = str2num(absTolField.Value);
    % User defined initial conditions
    x0 = str2num(x0EditField.Value);
    x0 = transpose(x0);
    % Torques
    u = str2num(uEditField.Value);
    u = transpose(u);
    
    % user input checks to make sure the dimensions are correct for the number of joints in the system
    if length(Kr_Val) ~= NumberOfJoints
        errordlg('The number of gear ratios does not match the number of joints in the system.', 'Error', 'modal');
        return;
    elseif length(MassLI_Val) ~= NumberOfJoints
        errordlg('The number of link masses does not match the number of joints in the system.', 'Error', 'modal');
        return;
    elseif length(MassMI_Val) ~= NumberOfJoints
        errordlg('The number of motor masses does not match the number of joints in the system.', 'Error', 'modal');
        return;
    elseif length(IntertiaLI_Val) ~= NumberOfJoints
        errordlg('The number of link inertias does not match the number of joints in the system.', 'Error', 'modal');
        return;
    elseif length(IntertiaMI_Val) ~= NumberOfJoints
        errordlg('The number of motor inertias does not match the number of joints in the system.', 'Error', 'modal');
        return;
    elseif length(FrictionS_Val) ~= NumberOfJoints
        errordlg('The number of static frictions does not match the number of joints in the system.', 'Error', 'modal');
        return;
    elseif length(FrictionV_Val) ~= NumberOfJoints
        errordlg('The number of viscosity frictions does not match the number of joints in the system.', 'Error', 'modal');
        return;
    elseif length(JointLengths) ~= NumberOfJoints
        errordlg('The number of joint lengths does not match the number of joints in the system.', 'Error', 'modal');
        return;
    elseif length(CenterOfMassLengths) ~= NumberOfJoints
        errordlg('The number of center of mass lengths does not match the number of joints in the system.', 'Error', 'modal');
        return;
    elseif length(u) ~= NumberOfJoints
        errordlg('The number of torques does not match the number of joints in the system.', 'Error', 'modal');
        return;
    elseif length(x0) ~= NumberOfJoints*2
        errordlg('The number of initial conditions does not match the number of joints in the system. x0 is defined as [position, velocity] for each joint.', 'Error', 'modal');
        return;
    end

    % This is the dumbest thing I have ever had to do in order to properly 
    % substitute in all symbolic variables with the correct user input values

    % define sym for center of mass symbolic as L1, L2, L3, etc. for number of joints so the lengths are equivalent to the number of user inputs
    for i = 1:NumberOfJoints
        CenterOfMassSymbolic(i) = "L"+i;
    end

    TaoSyms = sym('Tao_', [NumberOfJoints, 1]);
    u_sym = sym('u', [NumberOfJoints, 1]);
    % user input array
    user_input = [Kr_Val, MassLI_Val, MassMI_Val, IntertiaLI_Val, IntertiaMI_Val, FrictionS_Val, FrictionV_Val, JointLengths, CenterOfMassLengths];
    % corresponding symbolic variables
    user_input_sym = [str2sym(Kr_Sym), str2sym(MassLI_Sym), str2sym(MassMI_Sym), str2sym(IntertiaLI_Sym), str2sym(IntertiaMI_Sym), str2sym(FrictionS_Sym), str2sym(FrictionV_Sym), str2sym(JointSymbolic), str2sym(CenterOfMassSymbolic)];
    % append u_sym to user input and TaoSyms to user input syms for number of joints
    for i = 1:NumberOfJoints
        user_input = [user_input, u_sym(i)];
        user_input_sym = [user_input_sym, TaoSyms(i)];
    end

    for i = 1:NumberOfJoints
        vars_in_eqs = symvar(Tao(i));
        % disp('Variables in Equation ' + i + ':');
        % disp(vars_in_eqs);

        % Update the equations with the user input values based on the vars_in_eqs array referencing the user_input array
        for j = 1:length(vars_in_eqs)
            % if vars_in_eqs(j) contains user_input_sym, find the index of the variable in the user_input_sym array
            if ismember(vars_in_eqs(j), user_input_sym)
                % find the index of the variable in the user_input_sym array
                idx = find(vars_in_eqs(j) == user_input_sym);
                % update the equation with the user input value
                Tao(i) = subs(Tao(i), vars_in_eqs(j), user_input(idx));
            end
        end
    end
    
    disp('Updated Equations from User Input:')
    % Display the updated equations
    for i = 1:NumberOfJoints
        disp(['Equation ' num2str(i) ':']);
        disp(Tao(i));
    end

    % Solve it plz
    % okay :)
    % solves it with malicious intent

    q = str2sym(q);
    qd = str2sym(qd);
    qdd = str2sym(qdd);

    % preallocate the state variables
    x = sym('x', [NumberOfJoints*2, 1]);
    dx = sym('dx', [NumberOfJoints*2, 1]);

    % Define the state variables
    for i = 1:NumberOfJoints
        x(i) = q(i);
        x(i+NumberOfJoints) = qd(i);

        dx(i) = qd(i);
        dx(i+NumberOfJoints) = simplify(solve(Tao(i), qdd(i)));

    end

    % disp('dx:')
    % disp(dx)

    % Define the state matrix (A)
    A = jacobian(dx, x);
    % Define the input matrix (B)
    B = jacobian(dx, u_sym);

    % disp('State Variables:')
    % disp(x)
    % disp('Initial Conditions:')
    % disp(x0)
    % disp('State Matrix (A):')
    % disp(A)
    % disp('Input Matrix (B):')
    % disp(B)

    % LQR Stuff
    % check if the LQR checkbox is checked
    if lqrCheckBox.Value == 1  
        % check if QEditField and REditField are empty
        if isempty(Q_Val) || isempty(R_Val)
            disp('No User Input for Q and R')
            disp('Using Default Values, 500 for Q and 100 for R')
            % default values
            Q = 500*eye(NumberOfJoints*2);
            R = 100*eye(NumberOfJoints);
        else
            Q = Q_Val*eye(NumberOfJoints*2);
            R = R_Val*eye(NumberOfJoints);
        end

        % change the first half of Q to 5 times its value (Velocities)
        for i = 1:NumberOfJoints
            Q(i,i) = 5*Q(i,i);
        end
    else
        % Q and R are zero
        Q = 0;
        R = 0;
    end

    % Make some of these options user inputs :) (Done did it)
    options = odeset('MaxStep',maxStep_Val,'RelTol',relTol_Val,'AbsTol',absTol_Val*ones(1,(NumberOfJoints*2)));

    tf = str2num(tf.Value);
    tspan = [0 tf]; % Final simulation time
    u = @(t) u; % Make u a function of time

    % Solve the differential equations ode45 for state space
    [t, x] = ode45(@(t, x) systemDynamics(t, tf, x, u(t), A, B, FrictionS_Sym, FrictionS_Val, Q, R, lqrCheckBox.Value, dynamicsFig), tspan, x0, options);

    % pos of the joints over time
    pos = zeros(NumberOfJoints, length(t));
    for i = 1:NumberOfJoints
        pos(i,:) = x(:,i);
    end
    
    % Velocity of the joints over time
    vel = zeros(NumberOfJoints, length(t));
    for i = 1:NumberOfJoints
        vel(i,:) = x(:,i+NumberOfJoints);
    end

    tSpace = linspace(0, tf, size(vel, 2));
    % Calculate acceleration
    acc = zeros(NumberOfJoints, length(t));
    for i = 1:NumberOfJoints
        vel_diff = [0, diff(vel(i, :))];
        time_diff = diff(tSpace);
        acc_diff = vel_diff ./ [time_diff, 1];  % Append 1 to make time_diff the same size as vel_diff
        acc(i, :) = acc_diff;
    end

    setappdata(dynamicsFig, 'tSpace', tSpace);
    setappdata(dynamicsFig, 'pos', pos);
    setappdata(dynamicsFig, 'vel', vel);
    setappdata(dynamicsFig, 'acc', acc);

    plotButton_Callback(dynamicsFig,TransMats_Joint2Joint, q, JointSymbolic, JointLengths);
    disp('Simulation complete');
end

% plotting function callback for pos vel acc
function plotButton_Callback(dynamicsFig,TransMats_Joint2Joint, q, JointSymbolic, JointLengths)
    % Retrieve the results
    t = getappdata(dynamicsFig, 'tSpace');
    pos = getappdata(dynamicsFig, 'pos');
    vel = getappdata(dynamicsFig, 'vel');
    acc = getappdata(dynamicsFig, 'acc');
    
    % pos is the positions of each joint over time
    % vel is the velocities of each joint over time
    % acc is the accelerations of each joint for each time step

    % plotting needs to be robust for the number of joints
    % plot the positions of each joint over time

    % Create a new figure window with a larger size
    plotFig = uifigure('Name', 'Dynamics Simulation Plots', 'NumberTitle', 'off', 'Position', [100 100 600 500]);

    % Add a button to close the plot figure
    closeButton = uibutton(plotFig, 'Position', [10, 10, 100, 22], 'Text', 'Close', 'ButtonPushedFcn', @(btn,event) close(plotFig));
    closeButton.BackgroundColor = [0.7 0.2 0.2]; % Red color

    % Create axes for the position plot
    posAxes = uiaxes(plotFig, 'Position', [50 345 500 125]);
    % Plot the positions of each joint over time
    legendEntries = cell(1, size(pos, 1)); % Initialize cell array for legend entries
    for i = 1:size(pos, 1)
        plot(posAxes, t, pos(i,:));
        hold(posAxes, 'on');
        legendEntries{i} = sprintf('Joint %d', i); % Add entry to legend
    end
    legend(posAxes, legendEntries); % Add legend to plot
    title(posAxes, 'Joint Positions Over Time');
    xlabel(posAxes, 'Time (s)');
    ylabel(posAxes, 'Position (rad)');

    % Create axes for the velocity plot
    velAxes = uiaxes(plotFig, 'Position', [50 205 500 125]);
    % Plot the velocities of each joint over time
    legendEntries = cell(1, size(vel, 1)); % Initialize cell array for legend entries
    for i = 1:size(vel, 1)
        plot(velAxes, t, vel(i,:));
        hold(velAxes, 'on');
        legendEntries{i} = sprintf('Joint %d', i); % Add entry to legend
    end
    legend(velAxes, legendEntries); % Add legend to plot
    title(velAxes, 'Joint Velocities Over Time');
    xlabel(velAxes, 'Time (s)');
    ylabel(velAxes, 'Velocity (rad/s)');

    % Create axes for the acceleration plot
    accAxes = uiaxes(plotFig, 'Position', [50 65 500 125]);
    % Plot the accelerations of each joint over time
    legendEntries = cell(1, size(acc, 1)); % Initialize cell array for legend entries
    for i = 1:size(acc, 1)
        plot(accAxes, t, acc(i,:));
        hold(accAxes, 'on');
        legendEntries{i} = sprintf('Joint %d', i); % Add entry to legend
    end
    legend(accAxes, legendEntries); % Add legend to plot
    title(accAxes, 'Joint Accelerations Over Time');
    xlabel(accAxes, 'Time (s)');
    ylabel(accAxes, 'Acceleration (rad/s^2)');

    % Initialize link position vectors
    x = zeros(length(q), size(pos, 2));
    y = zeros(length(q), size(pos, 2));
    z = zeros(length(q), size(pos, 2));

    % Calculate position of each link for every end effector position
    for k = 1:size(pos, 2)
        % joint positions
        pos_current = pos(:, k);

        TransSyms = [q, str2sym(JointSymbolic)];
        TransVals = [transpose(pos_current), JointLengths];

        % Calculate position of each link
        for i = 1:length(q)
            % Extract the position from the transformation matrix
            positions = TransMats_Joint2Joint(1:3, 4, i);

            % Substitute symbolic variables with their corresponding numeric values
            positions = subs(positions, TransSyms, TransVals);

            x(i, k) = positions(1);
            y(i, k) = positions(2);
            z(i, k) = positions(3);
        end
    end

    % currentFig = gcf;
    % Create a UI alert dialog box
    uialert(plotFig, 'Press okay to playback dynamics', 'Start Animation', 'Icon', 'info', 'CloseFcn', @(~,~) startAnimation(q,x,y,z));
end

% Define the function to start the animation
function startAnimation(q,x,y,z)
    % Plot robot
    figure;
    hold on;  % Hold the current plot
    h = gobjects(1, length(q));  % Initialize array of graphics objects
    l = gobjects(1, length(q));  % Initialize array of line objects

    % Frame 0
    plot3(0, 0, 0, 'o-');
    % Plot a line from the origin to the first link and store the line object
    originLine = line([0, x(1, 1)], [0, y(1, 1)], [0, z(1, 1)], 'Color', 'k');

    for i = 1:length(q)
        h(i) = plot3(x(i, 1), y(i, 1), z(i, 1), 'o-');  % Plot each link
    end

    % Plot lines connecting the joints
    for i = 1:length(q)-1
        l(i) = line([x(i, 1), x(i+1, 1)], [y(i, 1), y(i+1, 1)], [z(i, 1), z(i+1, 1)], 'Color', 'k');
    end

    % Set x, y, and z limits
    xlim([min(x(:))-1, max(x(:))+1]);
    ylim([min(y(:))-1, max(y(:))+1]);
    zlim([min(z(:))-1, max(z(:))+1]);

    % Set the view to top down
    view(3);

    % Generate legend labels
    labels = cell(1, length(q) + 1);
    labels{1} = 'Origin';
    for i = 1:length(q)
        labels{i + 1} = sprintf('Link %d', i);
    end

    % Add a legend
    legend(labels, 'Location', 'best');

    hold off;  % Release the current plot
    grid on;
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title('Robot Operational Space Animation');

    % Animation loop
    for k = 2:size(x, 2)
        % Update the line from the origin to the first link
        originLine.XData = [0, x(1, k)];
        originLine.YData = [0, y(1, k)];
        originLine.ZData = [0, z(1, k)];

        for i = 1:length(q)
            % Update plot data for each link
            h(i).XData = x(i, k);
            h(i).YData = y(i, k);
            h(i).ZData = z(i, k);
        end

        % Update the lines connecting the joints
        for i = 1:length(q)-1
            l(i).XData = [x(i, k), x(i+1, k)];
            l(i).YData = [y(i, k), y(i+1, k)];
            l(i).ZData = [z(i, k), z(i+1, k)];
        end

        % Update the plot
        drawnow;

    end
end

% END OF TYLERS DYNAMICAL WONDERLAND