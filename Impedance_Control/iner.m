function  B = iner(c_2)
%INER   Inertia matrix of two-link planar arm.


global pi_m a k_r2

B(1,1) = a(1)*pi_m(1) + pi_m(2) + (a(2) + 2*a(1)*c_2)*pi_m(3) + pi_m(4);

B(1,2) = (a(2) + a(1)*c_2)*pi_m(3) + pi_m(4) + k_r2*pi_m(5);

B(2,1) = B(1,2);

B(2,2) = a(2)*pi_m(3) + pi_m(4) + k_r2*k_r2*pi_m(5);
