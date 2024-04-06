function G = FindG(GDir, M_L,M_M,JP_L,JP_M,Num,g)
G = sym(zeros(Num,1));
switch GDir
    case 'x'
        G_not = sym([g,0,0]');
    case 'y'
        G_not = sym([0,g,0]');
    case 'z'
        G_not = sym([0,0,g]');
    case '-x'
        G_not = sym([-g,0,0]');
    case '-y'
        G_not = sym([0,-g,0]');
    case '-z'
        G_not = sym([0,0,-g]');
end

for i = 1:Num
    tempG = sym([0]);
    for j = 1:Num
        LinkGj = M_L(j)*(G_not')*JP_L(:,i,j);
        MotorGj = M_M(j)*(G_not')*JP_M(:,i,j);
        tempG = tempG - (LinkGj+MotorGj);
    end
    G(i,:) = tempG;
end
G = simplify(G);
end

