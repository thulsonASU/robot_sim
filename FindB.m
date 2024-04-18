function [B] = FindB(M_LI,M_MI,JP_LI,JO_LI,I_LI,I_MI,Trans,JP_MI,JO_MI,Num)
B = zeros(Num);
for i = 1:Num
    TempB1 = simplify((M_LI(i)*(JP_LI(:,:,i)')*JP_LI(:,:,i)));
    TempB2 = simplify((JO_LI(:,:,i)')*Trans(1:3,1:3,i)*I_LI(i)*(Trans(1:3,1:3,i)')*JO_LI(:,:,i));
    TempB3 =  simplify(M_MI(i)*(JP_MI(:,:,i)')*JP_MI(:,:,i));
    TempB4 = simplify((JO_MI(:,:,i)')*Trans(1:3,1:3,i)*I_MI(i)*(Trans(1:3,1:3,i)')*JO_MI(:,:,i));
    TempB = simplify(TempB1+TempB2+TempB3+TempB4);
    B = B + TempB;
end
B = simplify(B);
end

