function [Jacobian] = Jpj_Mi(i,num_joints,joint_types,Transformations)
for j = 1:num_joints
    if j < i
        %If it is less than or the same as i we actually calculate
        %stuff
        if joint_types(j) == 'p'
            %Zj-1
            if j == 1
                %Z0 is always this
                Jacobian(:,j) = sym([0 0 1]');
            else
                %Find Zj-1 from our trans matricies
                Jacobian(:,j) = Transformations(1:3,3,j-1);
            end
            
        else
            %REVOLTE
            if j == 1
                Jacobian(:,j) = cross([0,0,1]',...
                    (Transformations(1:3,4,i-1)...
                    - [0,0,0]'));
            else
                Jacobian(:,j) = cross(Transformations(1:3,3,j-1),...
                    (Transformations(1:3,4,i-1)...
                    - Transformations(1:3,3,j-1)));
            end
        end
    else
        %If j >= i then it is just all 0's
        Jacobian(:,j) = sym([0 0 0]');
    end
    
end
Jacobian = sym(Jacobian);
end

