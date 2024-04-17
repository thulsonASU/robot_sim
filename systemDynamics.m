function dxdt = systemDynamics(t, x, u, A, B, FrictionV_Sym, FrictionV_Val, xd, Q, R, lqr_control_law)

    % FrictionV_Sym is the symbolic representation of the Fv friction from the governing equations the value is the user defined input
    % We need to update all FrictionV_Sym to FrictionV_Val*sign(qdot) in the A matrix if FricitonV_Sym exists in the A matrix
    % We need to update all FrictionV_Sym to FrictionV_Val*sign(qdot) in the B matrix if FricitonV_Sym exists in the B matrix
    % qdot is the x y z of the velocity of the system
    qdot = x(4:6);

    % Substitute the variables in A
    for i = 1:size(A, 1)
        for j = 1:size(A, 2)
            % Find the symbolic variables in A(i, j)
            vars_in_Aij = symvar(A(i, j));
            
            % If there are any symbolic variables, substitute them
            if ~isempty(vars_in_Aij)
                
                % Update FrictionV_Sym to FrictionV_Val*sign(qdot) if it exists in A(i, j)
                if ismember(FrictionV_Sym, vars_in_Aij)
                    A(i, j) = subs(A(i, j), FrictionV_Sym, FrictionV_Val*sign(qdot));
                    % remove FrictionV_Sym from vars_in_Aij
                    vars_in_Aij = setdiff(vars_in_Aij, FrictionV_Sym);
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
                
                % Update FrictionV_Sym to FrictionV_Val*sign(qdot) if it exists in B(i, j)
                if ismember(FrictionV_Sym, vars_in_Bij)
                    B(i, j) = subs(B(i, j), FrictionV_Sym, FrictionV_Val*sign(qdot));
                    % remove FrictionV_Sym from vars_in_Bij
                    vars_in_Bij = setdiff(vars_in_Bij, FrictionV_Sym);
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
        error = x - xd;
        % Get the gain matrix K from LQR
        K = lqr(A, B, Q, R);
        % calculate the control law
        u = -K*error;
    end

    % calculate the state space
    dxdt = A*x + B*u;

    % round t to the nearest 0.1
    t_rounded = round(t, 1);
    
    % print dxdt when time is approximately a whole number
    if abs(mod(t_rounded, 0.25)) < 1e-6  % adjust the tolerance as needed
        disp('Running...')
        disp('dxdt at time t: ');
        disp(t);
        disp(transpose(dxdt));
    end

end