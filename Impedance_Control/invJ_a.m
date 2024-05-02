function y = invJ_a(w)
%INVJ_A Joint velocity with Jacobian inverse for two-link planar arm.
%       y = INVJ_A(w) returns joint velocity as:
%
%       y = inv(J_a(w(3:6))*w(1:2)

y = inv(J_a(w(3:6)))*w(1:2);
