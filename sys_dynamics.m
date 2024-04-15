function dxdt = sys_dynamics(t, x, u)
    % x = [x1, x2, ..., xn]
    % u = [u1, u2, ..., um]
    % dxdt = [f1(x, u); f2(x, u); ...; fn(x, u)];
    n = length(x);
    dxdt = zeros(n, 1);
    
    for i = 1:n
        dxdt(i) = fi(x, u);
    end
    
end