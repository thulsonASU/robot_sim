function J = Jacobian(q,l1,l3)
    t1 = q(1);
    d3 = q(3);


    J = [-sin(t1)*(l1 + l3 + d3) 0 cos(t1);
             cos(t1)*(l1 + l3 + d3) 0 sin(t1);
                                0 1       0;
                                0 0       0;
                                0 0       0;
                                1 0       0];
    J = double(J);
end