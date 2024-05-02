function y = dJ_a(w)
%DJ_A   Derivative of Jacobian times joint velocity.
%       y = DJ_A(w) returns for a two-link planar arm:
%
%           dJ_a(w(3:6))
%       y = ------------*w(1:2)
%               dt

dJ(:,2) = [-w(4); -w(6)]*(w(1)+w(2));
dJ(:,1) = [-w(3); -w(5)]*w(1) + dJ(:,2);

y = dJ*w(1:2);
