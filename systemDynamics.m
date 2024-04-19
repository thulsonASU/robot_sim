function dxdt = systemDynamics(t, tf, x, u, A, B, FrictionS_Sym, FrictionS_Val, Q, R, lqr_control_law)

    % FrictionS_Sym is the symbolic representation of the Fv friction from the governing equations the value is the user defined input
    % We need to update all FrictionS_Sym to FrictionS_Val*sign(qdot) in the A matrix if FricitonV_Sym exists in the A matrix
    % We need to update all FrictionS_Sym to FrictionS_Val*sign(qdot) in the B matrix if FricitonV_Sym exists in the B matrix
    % qdot is the x y z of the velocity of the system
    NumberOfJoints = size(B, 2); % number of joints in the system
    qdot = x(((NumberOfJoints*2)/2)+1:NumberOfJoints*2); % velocity of the system

    % Capture the acceleration of the system

    % Substitute the variables in A
    for i = 1:size(A, 1)
        for j = 1:size(A, 2)
            % Find the symbolic variables in A(i, j)
            vars_in_Aij = symvar(A(i, j));
            
            % If there are any symbolic variables, substitute them
            if ~isempty(vars_in_Aij)
                
                % Update FrictionS_Sym to FrictionS_Val*sign(qdot) if it exists in A(i, j)
                if ismember(FrictionS_Sym, vars_in_Aij)
                    A(i, j) = subs(A(i, j), FrictionS_Sym, FrictionS_Val*sign(qdot));
                    % remove FrictionS_Sym from vars_in_Aij
                    vars_in_Aij = setdiff(vars_in_Aij, FrictionS_Sym);
                end

                % Get the indices of the variables in x
                indices = arrayfun(@(v) str2double(regexp(char(v), '\d+', 'match')), vars_in_Aij);

                % Substitute the remaining variables with the corresponding values from x
                A(i, j) = subs(A(i, j), vars_in_Aij, transpose(x(indices)));
            end
        end
    end

    % Substitute the variables in B
    for i = 1:size(B, 1)
        for j = 1:size(B, 2)
            % Find the symbolic variables in B(i, j)
            vars_in_Bij = symvar(B(i, j));
            
            % If there are any symbolic variables, substitute them
            if ~isempty(vars_in_Bij)
                
                % Update FrictionS_Sym to FrictionS_Val*sign(qdot) if it exists in B(i, j)
                if ismember(FrictionS_Sym, vars_in_Bij)
                    B(i, j) = subs(B(i, j), FrictionS_Sym, FrictionS_Val*sign(qdot));
                    % remove FrictionS_Sym from vars_in_Bij
                    vars_in_Bij = setdiff(vars_in_Bij, FrictionS_Sym);
                end

                % Get the indices of the variables in x
                indices = arrayfun(@(v) str2double(regexp(char(v), '\d+', 'match')), vars_in_Bij);

                % Substitute the remaining variables with the corresponding values from x
                B(i, j) = subs(B(i, j), vars_in_Bij, transpose(x(indices)));
            end
        end
    end

    % convert A and B to double
    A = double(A);
    B = double(B);

    % check if lqr control law is true
    if lqr_control_law == 1
        % error from desired to current state
        % error = x - xd;
        % Get the gain matrix K from LQR
        K = lqr(A, B, Q, R);
        % calculate the control law
        u = -K*x; % for some reason it works only as a regulator not a tracker
    end

    % calculate the state space
    dxdt = A*x + B*u;

    %ode45 uses the state space model to calculate the next state
    %delta_t is not consistent
    %I need to calculate a new delta_t based on the current time
    
    % persistent accMatrix acc t_prev firstPass;
    % if isempty(firstPass)
    %     t_prev = 0;
    %     accMatrix = [];
    %     firstPass = 1;
    % end
    
    % % calculate the next time step
    % delta_t = t - t_prev;
    % t_prev = t;

    % % first two elements of dxdt are the velocity of the system
    % % calculate the acceleration of the system
    % acc = qdot/delta_t;
    % accMatrix = [accMatrix, acc];

    % round t to the nearest 0.1
    t_rounded = round(t, 1);
    
    % print dxdt when time is approximately a whole number
    if abs(mod(t_rounded, 0.25)) < 1e-6  % adjust the tolerance as needed
        disp('Running...')
        disp('dxdt at time t: ');
        disp(t);
        disp(transpose(dxdt));
    end

    % if t == tf
    %     disp('Time Matrix Updated in Workspace.')
    %     assignin('base', 'accMatrix', accMatrix);
    %     assignin('base', 't', t);
    % end

end