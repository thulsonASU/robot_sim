function x = J_adq(w)
%J_ADQ  Tip velocity for two-link planar arm.
%       x = J_ADQ(w) returns tip velocity as:
%
%       x = J_a(w(3:6))*w(1:2)
%
%       where:
%
%       w(3:6)=[a(1)*c_1;a(2)*c_12;a(1)*s_1;a(2)*s_12]

x = J_a(w(3:6))*w(1:2);
