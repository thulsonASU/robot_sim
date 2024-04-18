function [Jacobian] = Jpj_Li(i,num_joints,joint_types,Transformations,jointlen,jointhalves)
    for j = 1:num_joints
        if j <= i
            %If it is less than or the same as i we actually calculate
            %stuff
            if joint_types(j) == 'p'
                %Zj-1
                if j == 1
                    %Z0 is always this
                    Jacobian(:,j) = [0 0 1]';
                else
                    %Find Zj-1 from our trans matricies
                    Jacobian(:,j) = Transformations(1:3,3,j-1);
                end
                
            else
                if j == 1
                    %Because j = 1 is special we need this if function
                    Jacobian(:,j) = cross([0,0,1]',...
                        (subs(Transformations(1:3,4,i),str2sym("a"+i),str2sym("L"+i))...
                        - [0,0,0]'));
                else
                    Jacobian(:,j) = cross(Transformations(1:3,3,j-1),...
                        (subs(Transformations(1:3,4,i),str2sym("a"+i),str2sym("L"+i))...
                         - Transformations(1:3,4,j-1)));
                end
            end
        else
            %If j > i then it is just all 0's
            Jacobian(:,j) = [0 0 0]';
        end
        
    end
    Jacobian = sym(Jacobian);
end

