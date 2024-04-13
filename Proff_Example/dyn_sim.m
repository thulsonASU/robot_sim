clear all 
close all
timeSpan = [0,50]; % Solve from t=0 to t=50
IC = [1.5;0;-0.5]; % initial conditions. Note that ic of y1dot (i.e. x2) is not specified, and hence assumed zero.
options = odeset('RelTol',1e-4,'AbsTol',[1e-5 1e-5 1e-5]);
[tsim,y] = ode45(@samplesys,timeSpan,IC); 
plot(tsim,y)
legend('y1','ydot1','y2')

% [tsim,y] = ode45(@samplesys,timeSpan,IC, options, extra_inputs); plot(tsim,y)


