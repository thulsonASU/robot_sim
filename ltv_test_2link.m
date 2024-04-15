m1 = 3;         % Mass of link 1, kg
m2 = 2;         % Mass of link 2, kg
L1 = 0.3;       % Length of link 1, meters
L2 = 0.3;       % Length of link 2, meters

T0 = 0;
Tf = 5;
tref = linspace(T0,Tf,400)';

io = [...
   linio('TwoLinkRobotSM/tau1',1,'input');...
   linio('TwoLinkRobotSM/tau2',1,'input');...
   linio('TwoLinkRobotSM/Mux',1,'output')];

tlin = linspace(T0,Tf,50);
op = findop('TwoLinkRobotSM',tlin);

linOpt = linearizeOptions('StoreOffsets',true);
[G,~,info] = linearize('TwoLinkRobotSM',io,op,linOpt);
G.u = 'tau';
G.y = 'x';
G.SamplingGrid = struct('Time',tlin);

Gltv = ssInterpolant(G,info.Offsets);
Gltv.StateName

Gltv = xperm(Gltv,[1 3 2 4]);
Gltv.StateName

xltv = lsim(Gltv,[tau1 tau2],tref,xinit);

figure
plot(tref,xltv(:,1),tref,xltv(:,2),tref,xltv(:,3),tref,xltv(:,4))
legend('q1','q2','q1dot','q2dot')
xlabel('Time (s)')
ylabel('State Variables')
title('Linear Time-Varying Simulation')