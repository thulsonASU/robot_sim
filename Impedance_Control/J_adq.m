function x = J_adq(w)
%J_ADQ  Tip velocity for two-link planar arm.
%       x = J_ADQ(w) returns tip velocity as:
%
%       x = J_a(w(3:6))*w(1:2)


x = J_a(w(3:6))*w(1:2);
