function dxdt = systemDynamics_ss_test(t, x, u, A, B)

    % Substitute the variables in A
    for i = 1:size(A, 1)
        for j = 1:size(A, 2)
            % Find the symbolic variables in A(i, j)
            vars_in_Aij = symvar(A(i, j));
            
            % If there are any symbolic variables, substitute them
            if ~isempty(vars_in_Aij)
                % Get the indices of the variables in x
                indices = arrayfun(@(v) str2double(regexp(char(v), '\d+', 'match')), vars_in_Aij);
                
                % % display vars_in_Aij and indices
                % disp(vars_in_Aij);
                % disp(x(indices));

                % Substitute the variables with the corresponding values from x
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
                % Get the indices of the variables in x
                indices = arrayfun(@(v) str2double(regexp(char(v), '\d+', 'match')), vars_in_Bij);

                % % display vars_in_Aij and indices
                % disp(vars_in_Bij);
                % disp(x(indices));
                
                % Substitute the variables with the corresponding values from x
                B(i, j) = subs(B(i, j), vars_in_Bij, transpose(x(indices)));
            end
        end
    end

    % convert A and B to double
    A = double(A);
    B = double(B);

    disp(class(A))
    disp(class(B))
    disp(class(x))
    disp(class(u))

    % calculate the state space
    dxdt = A*x + B*u;

    % print dydt when time is approximately a whole number
    if abs(mod(t, 1)) < 1e-6  % adjust the tolerance as needed
        disp('dydt at time t: ');
        disp(t);
        disp(transpose(dxdt));
    end

end