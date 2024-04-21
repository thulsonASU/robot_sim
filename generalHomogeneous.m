function H = generalHomogeneous(dhTable)
    
    % Input: DH table row
    a = dhTable(1,1);
    alpha = dhTable(1,2);
    d = dhTable(1,3);
    theta = dhTable(1,4);

    % all angles are in radians
    % Define the generalized homogeneous transformation matrix
    H = [cos(theta), -sin(theta)*cos(alpha), sin(theta)*sin(alpha), a*cos(theta);
                           sin(theta), cos(theta)*cos(alpha), -cos(theta)*sin(alpha), a*sin(theta);
                           0, sin(alpha), cos(alpha), d;
                           0, 0, 0, 1];
end