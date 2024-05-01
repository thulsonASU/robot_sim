function B = B_Lagrangian(q,b2,b3,m3,l1,m1,m2,l2,l3,r1)
    d3 = q(3);
    

    B = [[(b2^2*m2)/12 + (b3^2*m3)/12 + d3^2*m3 + l1^2*m1 + l1^2*m2 + l1^2*m3 + (l2^2*m2)/12 + (4*l3^2*m3)/3 + 6*m1*r1^2 + 2*d3*l1*m3,       0,  0]
        [                                                                                                                          0, m2 + m3,  0]
        [                                                                                                                          0,       0, m3]];

end