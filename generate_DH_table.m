function dh_table = generate_DH_table(num_joints, joint_types)
    dh_table = sym(zeros(num_joints, 4));

    for i = 1:num_joints
        % Set joint type
        if strcmpi(joint_types(i), 'r') || strcmpi(joint_types(i), 'p')
            % Prompt user for input
            input_str = input(['Enter DH parameters for joint ' num2str(i) ' in the format [a, alpha, d, theta]: '], 's');
            
            % Remove brackets if they exist
            input_str = strrep(input_str, '[', '');
            input_str = strrep(input_str, ']', '');
            
            % Parse input string
            inputs = strsplit(input_str, ',');
            
            % Convert each input to either symbolic or numeric
            for j = 1:length(inputs)
                input_val = strtrim(inputs{j});
                if isempty(input_val)
                    error('Invalid input. Please enter all DH parameters.');
                end
                
                if ~isnan(str2double(input_val)) % Check if it's a number
                    inputs{j} = str2double(input_val); % Convert to number
                else
                    inputs{j} = str2sym(input_val); % Convert to symbolic variable using str2sym
                end
            end
            
            % Store DH parameters in the table
            dh_table(i, :) = [inputs{:}];
        else
            error('Unknown joint type. Please specify either "r" for revolute or "p" for prismatic.');
        end
    end
end
