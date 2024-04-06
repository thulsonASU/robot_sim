function [Jacobian] = Joj_Li(i,num_joints,joint_types,Transformations)
    for j = 1:num_joints
        if j <= i
            if joint_types(j) == 'p'
                %Prismatic
                Jacobian(:,j) = sym([0,0,0]');
            else
                %Revolute
                if j == 1
                    %Z0 is always this
                    Jacobian(:,j) = sym([0 0 1]');
                else
                    %Find Zj-1 from our trans matricies
                    Jacobian(:,j) = Transformations(1:3,3,j-1);
                end
            end
        else
            %If j > i then it is just all 0's
            Jacobian(:,j) = sym([0 0 0]');
        end
        
    end
    Jacobian = sym(Jacobian);
end

