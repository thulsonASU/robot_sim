function x = Kinematics(q,l3,l1,l2)
    t1 = q(1);
    d2 = q(2);
    d3 = q(3);
    % H0_3 computed with the main.m
    H = [[ 0, -sin(t1), cos(t1), cos(t1)*(d3 + l3) + l1*cos(t1)]
        [ 0,  cos(t1), sin(t1), sin(t1)*(d3 + l3) + l1*sin(t1)]
        [-1,        0,       0,                        d2 + l2]
        [ 0,        0,       0,                              1]];

    H = double(H);
    
    % Return x with positions and orientations
    x = [ H(1:3,4); rotm2eul(H(1:3,1:3), 'ZYZ')'];
end