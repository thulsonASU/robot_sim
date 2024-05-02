function Ja = J_a(u)
%J_A    Jacobian of two-link planar arm.
%       Ja = J_A(u) returns 2-by-2 analytical Jacobian, where:

Ja = eye(2);
Ja(:,2) = [-u(4); u(2)];
Ja(:,1) = [-u(3); u(1)] + Ja(:,2);
