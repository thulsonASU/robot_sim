function q = InverseKinematics(p,l2,l1,l3)
    x = p(1);
    y = p(2);
    z = p(3);
    
    
    q2 = z-l2;
    q1 = atan2(y,x);
    q3 = -l1 + sqrt(x^2+y^2) - l3;
    
    q = [q1; q2; q3];
end