function dx = samplesys(t,x, extra_inputs) 
a0 = 1; a1 = 1; a2 = 1;
dx = zeros(3,1);
dx(1) = x(2);
dx(2) = 0.1*sin(t) -a1*x(2) - a0*(x(1)+x(3)); dx(3) = 5 -a2*(x(3)-x(1));