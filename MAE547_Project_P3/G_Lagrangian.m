function G = G_Lagrangian(q,m3,l1,m1,m2)
    t12 = q(1);
    d31 = q(2);
    g = 9.81;
    % Calculate G
    G = [ (g*cos(t12)*(2*d31*m3 + l1*m1 + 2*l1*m2 + 2*l1*m3))/2; 0; g*m3*sin(t12)];
end

